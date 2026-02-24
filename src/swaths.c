#define TOPOTOOLBOX_BUILD

#include <float.h>
#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "topotoolbox.h"

// 8-connected neighbor offsets (row, col) — shared by all D8 operations.
static const int k_di8[8] = {-1, -1, -1, 0, 0, 1, 1, 1};
static const int k_dj8[8] = {-1, 0, 1, -1, 1, -1, 0, 1};

/*
  Swath profile analysis — libtopotoolbox
  ========================================

  Three analysis modes, all operating in pixel-space float coordinates:

  Transverse  — collapses the swath into a single cross-sectional profile
                binned by signed perpendicular distance from the track.

  Longitudinal — per-track-point statistics using either:
    Case 1 (binning_distance <= 0): Bresenham cross-section at each point.
    Case 2 (binning_distance  > 0): Meijster EDT nearest-track-point
                assignment + sliding-window accumulation along the track.

  Windowed longitudinal — per-track-point statistics gathered from an
    oriented rectangle (PCA tangent, ±half_width × ±binning_distance).
    No pre-computed distance map required.

  Distance maps — two public variants:
    swath_compute_distance_map  : clipped to half_width, optional centre-line.
    swath_compute_full_distance_map : unbounded, optional mask/dem filter.

  Line simplification — swath_simplify_line (IEF engine, 7 stopping criteria).
*/

// ============================================================================
// Statistics accumulator — running mean, variance, min, max (online, O(1))
// ============================================================================

typedef struct {
  ptrdiff_t count;
  double sum;
  double sum_sq;
  float min_val;
  float max_val;
} swath_stats_accumulator;

// ============================================================================
// Percentile accumulator — collects values, sorts on demand, queries by rank
// ============================================================================

typedef struct {
  ptrdiff_t capacity;
  ptrdiff_t count;
  float *values;
} percentile_accumulator;

static inline void percentile_accumulator_init(percentile_accumulator *acc,
                                               ptrdiff_t initial_capacity) {
  acc->capacity = initial_capacity > 0 ? initial_capacity : 16;
  acc->count = 0;
  acc->values = (float *)malloc(acc->capacity * sizeof(float));
}

static inline void percentile_accumulator_add(percentile_accumulator *acc,
                                              float value) {
  if (isnan(value)) return;
  if (acc->count >= acc->capacity) {
    acc->capacity *= 2;
    acc->values = (float *)realloc(acc->values, acc->capacity * sizeof(float));
  }
  acc->values[acc->count++] = value;
}

static int compare_floats(const void *a, const void *b) {
  float fa = *(const float *)a;
  float fb = *(const float *)b;
  return (fa > fb) - (fa < fb);
}

static inline void percentile_accumulator_sort(percentile_accumulator *acc) {
  if (acc->count > 0) {
    qsort(acc->values, acc->count, sizeof(float), compare_floats);
  }
}

static inline float percentile_accumulator_get(
    const percentile_accumulator *acc, float p) {
  if (acc->count == 0) return NAN;
  if (p < 0.0f) p = 0.0f;
  if (p > 100.0f) p = 100.0f;
  float index = p / 100.0f * (acc->count - 1);
  ptrdiff_t lower = (ptrdiff_t)index;
  ptrdiff_t upper = lower + 1;
  if (upper >= acc->count) return acc->values[acc->count - 1];
  float fraction = index - lower;
  return acc->values[lower] * (1.0f - fraction) + acc->values[upper] * fraction;
}

static inline void percentile_accumulator_free(percentile_accumulator *acc) {
  free(acc->values);
  acc->values = NULL;
  acc->count = 0;
  acc->capacity = 0;
}

static inline void accumulator_init(swath_stats_accumulator *acc) {
  acc->count = 0;
  acc->sum = 0.0;
  acc->sum_sq = 0.0;
  acc->min_val = FLT_MAX;
  acc->max_val = -FLT_MAX;
}

static inline void accumulator_add(swath_stats_accumulator *acc, float value) {
  if (isnan(value)) return;
  acc->count++;
  acc->sum += value;
  acc->sum_sq += value * value;
  if (value < acc->min_val) acc->min_val = value;
  if (value > acc->max_val) acc->max_val = value;
}

static inline float accumulator_mean(const swath_stats_accumulator *acc) {
  return acc->count > 0 ? (float)(acc->sum / acc->count) : NAN;
}

static inline float accumulator_stddev(const swath_stats_accumulator *acc) {
  if (acc->count <= 1) return 0.0f;
  double variance =
      (acc->sum_sq - acc->sum * acc->sum / acc->count) / (acc->count - 1);
  return variance > 0.0 ? (float)sqrt(variance) : 0.0f;
}

// ============================================================================
// Min-heap — used by frontier Dijkstra and the IEF simplification engine
// ============================================================================

typedef struct {
  float abs_dist;  // absolute perpendicular distance in pixels (heap key)
  ptrdiff_t idx;   // linear pixel index
} heap_entry;

typedef struct {
  heap_entry *data;
  ptrdiff_t size;
  ptrdiff_t capacity;
} min_heap;

static void heap_init(min_heap *h, ptrdiff_t initial_capacity) {
  h->capacity = initial_capacity > 0 ? initial_capacity : 256;
  h->size = 0;
  h->data = (heap_entry *)malloc(h->capacity * sizeof(heap_entry));
}

static void heap_push(min_heap *h, float abs_dist, ptrdiff_t idx) {
  if (h->size >= h->capacity) {
    h->capacity *= 2;
    h->data = (heap_entry *)realloc(h->data, h->capacity * sizeof(heap_entry));
  }
  ptrdiff_t i = h->size++;
  h->data[i].abs_dist = abs_dist;
  h->data[i].idx = idx;
  // Sift up
  while (i > 0) {
    ptrdiff_t parent = (i - 1) / 2;
    if (h->data[parent].abs_dist <= h->data[i].abs_dist) break;
    heap_entry tmp = h->data[parent];
    h->data[parent] = h->data[i];
    h->data[i] = tmp;
    i = parent;
  }
}

static heap_entry heap_pop(min_heap *h) {
  heap_entry top = h->data[0];
  h->data[0] = h->data[--h->size];
  // Sift down
  ptrdiff_t i = 0;
  for (;;) {
    ptrdiff_t left = 2 * i + 1;
    ptrdiff_t right = 2 * i + 2;
    ptrdiff_t smallest = i;
    if (left < h->size && h->data[left].abs_dist < h->data[smallest].abs_dist)
      smallest = left;
    if (right < h->size && h->data[right].abs_dist < h->data[smallest].abs_dist)
      smallest = right;
    if (smallest == i) break;
    heap_entry tmp = h->data[smallest];
    h->data[smallest] = h->data[i];
    h->data[i] = tmp;
    i = smallest;
  }
  return top;
}

static void heap_free(min_heap *h) {
  free(h->data);
  h->data = NULL;
  h->size = 0;
  h->capacity = 0;
}

// ============================================================================
// Signed perpendicular distance from a point to a segment (pixel space)
// ============================================================================
// Returns positive if the point is left of a→b, negative if right.
// proj_x/proj_y: closest point on segment; lambda: clamped parameter in [0,1].
static float point_to_segment_distance(float px, float py, float ax, float ay,
                                       float bx, float by, float *proj_x,
                                       float *proj_y, float *lambda) {
  float dx = bx - ax;
  float dy = by - ay;
  float dpx = px - ax;
  float dpy = py - ay;
  float seg_length_sq = dx * dx + dy * dy;

  if (seg_length_sq < 1e-10f) {
    *lambda = 0.0f;
    *proj_x = ax;
    *proj_y = ay;
    return sqrtf(dpx * dpx + dpy * dpy);
  }

  float t = (dpx * dx + dpy * dy) / seg_length_sq;
  if (t < 0.0f)
    t = 0.0f;
  else if (t > 1.0f)
    t = 1.0f;

  *lambda = t;
  *proj_x = ax + t * dx;
  *proj_y = ay + t * dy;

  float to_p_x = px - *proj_x;
  float to_p_y = py - *proj_y;
  float dist = sqrtf(to_p_x * to_p_x + to_p_y * to_p_y);
  float cross = dpx * dy - dpy * dx;

  return cross >= 0.0f ? dist : -dist;
}

// ============================================================================
// Boundary Dijkstra — inward D8 wavefront from swath-edge pixels
// ============================================================================
//
// Used internally by swath_compute_distance_map to build the distance-from-
// boundary field needed for centre-line extraction.  Seeds are the boundary
// pixels on one side of the swath (positive or negative); expansion is
// restricted to pixels within the swath mask (best_abs <= hw_px).
static void boundary_dijkstra(float *restrict dist_out,
                              const float *restrict best_abs,
                              const ptrdiff_t *seeds, ptrdiff_t n_seeds,
                              ptrdiff_t dims[2], float hw_px) {
  ptrdiff_t total = dims[0] * dims[1];

  for (ptrdiff_t i = 0; i < total; i++) dist_out[i] = FLT_MAX;

  min_heap heap;
  heap_init(&heap, n_seeds * 4);

  for (ptrdiff_t s = 0; s < n_seeds; s++) {
    ptrdiff_t idx = seeds[s];
    dist_out[idx] = 0.0f;
    heap_push(&heap, 0.0f, idx);
  }

  static const float dd8[8] = {1.41421356f, 1.0f, 1.41421356f, 1.0f, 1.0f,
                               1.41421356f, 1.0f, 1.41421356f};

  while (heap.size > 0) {
    heap_entry top = heap_pop(&heap);
    ptrdiff_t idx = top.idx;

    if (top.abs_dist > dist_out[idx]) continue;  // stale

    ptrdiff_t ci = idx % dims[0];
    ptrdiff_t cj = idx / dims[0];

    for (int n = 0; n < 8; n++) {
      ptrdiff_t ni = ci + k_di8[n];
      ptrdiff_t nj = cj + k_dj8[n];
      if (ni < 0 || ni >= dims[0] || nj < 0 || nj >= dims[1]) continue;

      ptrdiff_t nidx = nj * dims[0] + ni;
      if (best_abs[nidx] == FLT_MAX || best_abs[nidx] > hw_px) continue;

      float new_dist = dist_out[idx] + dd8[n];
      if (new_dist < dist_out[nidx]) {
        dist_out[nidx] = new_dist;
        heap_push(&heap, new_dist, nidx);
      }
    }
  }

  heap_free(&heap);
}

// ============================================================================
// Frontier distance map — Dijkstra outward from track, bounded at max_dist_px
// ============================================================================
//
// Seeds: track segments rasterized at ~0.5 px spacing.
// Expansion: 8-connected, ordered by absolute perpendicular distance.
// Each pixel is settled once; stale heap entries are skipped.
//
// Outputs (caller-allocated, size dims[0]*dims[1]; NULL to skip):
//   best_abs    — absolute perpendicular distance in pixels (FLT_MAX =
//   unvisited) signed_dist — signed perpendicular distance in pixels
//   nearest_seg — index of nearest track segment (-1 = unvisited)

static void frontier_distance_map(
    float *restrict best_abs, float *restrict signed_dist,
    ptrdiff_t *restrict nearest_seg, const float *restrict track_i,
    const float *restrict track_j, ptrdiff_t n_track_points, ptrdiff_t dims[2],
    float max_dist_px, const float *dem, const int8_t *mask) {
  ptrdiff_t total = dims[0] * dims[1];

  // We always need best_abs for Dijkstra tracking
  int free_best = 0;
  if (best_abs == NULL) {
    best_abs = (float *)malloc(total * sizeof(float));
    free_best = 1;
  }

  // Initialize
  for (ptrdiff_t i = 0; i < total; i++) {
    best_abs[i] = FLT_MAX;
  }
  if (signed_dist) {
    for (ptrdiff_t i = 0; i < total; i++) signed_dist[i] = 0.0f;
  }
  if (nearest_seg) {
    for (ptrdiff_t i = 0; i < total; i++) nearest_seg[i] = -1;
  }

  min_heap heap;
  heap_init(&heap, n_track_points * 8);

// Helper: test a pixel against candidate segments, update if better
// clang-format off
#define TRY_PIXEL(pi, pj, seg_lo, seg_hi)                             \
  do {                                                                \
    if ((pi) >= 0 && (pi) < dims[0] && (pj) >= 0 && (pj) < dims[1]) { \
      ptrdiff_t _idx = (pj)*dims[0] + (pi);                           \
      if (dem && isnan(dem[_idx])) break;                             \
      if (mask && !mask[_idx]) break;                                 \
      float _px = (float)(pi);                                        \
      float _py = (float)(pj);                                        \
      float _best_sd = 0.0f, _best_ad = FLT_MAX;                      \
      ptrdiff_t _best_seg = -1;                                       \
      for (ptrdiff_t _sk = (seg_lo); _sk <= (seg_hi); _sk++) {        \
        float _proj_x, _proj_y, _lam;                                 \
        float _d = point_to_segment_distance(                         \
            _px, _py, track_i[_sk], track_j[_sk], track_i[_sk + 1],   \
            track_j[_sk + 1], &_proj_x, &_proj_y, &_lam);             \
        if (fabsf(_d) < _best_ad) {                                   \
          _best_ad = fabsf(_d);                                       \
          _best_sd = _d;                                              \
          _best_seg = _sk;                                            \
        }                                                             \
      }                                                               \
      if (_best_ad < best_abs[_idx] && _best_ad <= max_dist_px) {     \
        best_abs[_idx] = _best_ad;                                    \
        if (signed_dist) signed_dist[_idx] = _best_sd;                \
        if (nearest_seg) nearest_seg[_idx] = _best_seg;               \
        heap_push(&heap, _best_ad, _idx);                             \
      }                                                               \
    }                                                                 \
  } while (0)

  // Seed: rasterize track segments at ~0.5px spacing
  for (ptrdiff_t k = 0; k < n_track_points - 1; k++) {
    float di = track_i[k + 1] - track_i[k];
    float dj = track_j[k + 1] - track_j[k];
    float seg_len = sqrtf(di * di + dj * dj);
    ptrdiff_t n_steps = (ptrdiff_t)(seg_len * 2.0f) + 1;

    ptrdiff_t seg_lo = k > 0 ? k - 1 : 0;
    ptrdiff_t seg_hi = k + 1 < n_track_points - 1 ? k + 1 : k;

    for (ptrdiff_t s = 0; s <= n_steps; s++) {
      float t = n_steps > 0 ? (float)s / (float)n_steps : 0.0f;
      ptrdiff_t pi = (ptrdiff_t)(track_i[k] + t * di + 0.5f);
      ptrdiff_t pj = (ptrdiff_t)(track_j[k] + t * dj + 0.5f);
      TRY_PIXEL(pi, pj, seg_lo, seg_hi);
    }
  }

  // Dijkstra expansion
  while (heap.size > 0) {
    heap_entry top = heap_pop(&heap);
    ptrdiff_t idx = top.idx;

    // Skip stale entries
    if (top.abs_dist > best_abs[idx]) continue;

    ptrdiff_t ci = idx % dims[0];
    ptrdiff_t cj = idx / dims[0];
    ptrdiff_t seg = nearest_seg ? nearest_seg[idx] : 0;

    // Search range: parent's segment ±2
    ptrdiff_t seg_lo = seg > 1 ? seg - 2 : 0;
    ptrdiff_t seg_hi =
        seg + 2 < n_track_points - 1 ? seg + 2 : n_track_points - 2;

    for (int n = 0; n < 8; n++) {
      ptrdiff_t ni = ci + k_di8[n];
      ptrdiff_t nj = cj + k_dj8[n];
      TRY_PIXEL(ni, nj, seg_lo, seg_hi);
    }
  }

#undef TRY_PIXEL

  heap_free(&heap);
  if (free_best) free(best_abs);
}
// clang-format on

// ============================================================================
// Public API
// ============================================================================

// Returns the number of transverse bins covering [-half_width, +half_width]
// at the given bin_resolution.  Convenience function for pre-allocation.
TOPOTOOLBOX_API
ptrdiff_t swath_compute_nbins(float half_width, float bin_resolution) {
  if (bin_resolution <= 0.0f || half_width <= 0.0f) return 0;
  return (ptrdiff_t)(2.0f * half_width / bin_resolution) + 1;
}

// ============================================================================
// Clipped distance map — signed/unsigned, optional centre-line extraction
// ============================================================================
//
// Computes a perpendicular distance map clipped to half_width metres.
// Pixels beyond half_width receive NaN.
//
// If centre_line_i/j/width are non-NULL, also extracts the Voronoi ridge
// between the positive-side and negative-side swath boundaries (the
// geometric centre line) and the local swath width at each centre pixel.
// The ridge is thinned to 1-pixel width via elbow removal on the ordered path.
//
// Returns the number of centre-line pixels found (0 if not requested).

TOPOTOOLBOX_API
ptrdiff_t swath_compute_distance_map(
    float *restrict distance_from_track, ptrdiff_t *restrict nearest_segment,
    float *restrict dist_from_boundary, float *restrict centre_line_i,
    float *restrict centre_line_j, float *restrict centre_width,
    const float *restrict track_i, const float *restrict track_j,
    ptrdiff_t n_track_points, ptrdiff_t dims[2], float cellsize,
    float half_width, int compute_signed) {
  if (n_track_points < 2) return 0;

  ptrdiff_t total = dims[0] * dims[1];
  float hw_px = half_width / cellsize;

  // Allocate internal pixel-unit arrays
  float *best_abs = (float *)malloc(total * sizeof(float));
  float *signed_dist_px = (float *)malloc(total * sizeof(float));

  // Need nearest_seg for Dijkstra; share with caller if requested
  ptrdiff_t *nseg = NULL;
  int free_nseg = 0;
  if (nearest_segment != NULL) {
    nseg = nearest_segment;
  } else {
    nseg = (ptrdiff_t *)malloc(total * sizeof(ptrdiff_t));
    free_nseg = 1;
  }

  frontier_distance_map(best_abs, signed_dist_px, nseg, track_i, track_j,
                        n_track_points, dims, hw_px, NULL, NULL);

  // Fill distance_from_track output (meters)
  for (ptrdiff_t i = 0; i < total; i++) {
    if (best_abs[i] == FLT_MAX || best_abs[i] > hw_px) {
      distance_from_track[i] = NAN;
    } else if (compute_signed) {
      distance_from_track[i] = signed_dist_px[i] * cellsize;
    } else {
      distance_from_track[i] = best_abs[i] * cellsize;
    }
  }

  // Optional centre-line inversion
  ptrdiff_t n_centre = 0;
  int want_centre =
      (centre_line_i != NULL && centre_line_j != NULL && centre_width != NULL);
  int want_dfb = (dist_from_boundary != NULL);

  if (want_centre || want_dfb) {
    // Find boundary pixels, split by side
    static const int dn_i[4] = {-1, 1, 0, 0};
    static const int dn_j[4] = {0, 0, -1, 1};

    ptrdiff_t pos_cap = 4096, neg_cap = 4096;
    ptrdiff_t n_pos = 0, n_neg = 0;
    ptrdiff_t *pos_seeds = (ptrdiff_t *)calloc(pos_cap, sizeof(ptrdiff_t));
    ptrdiff_t *neg_seeds = (ptrdiff_t *)calloc(neg_cap, sizeof(ptrdiff_t));

    for (ptrdiff_t idx = 0; idx < total; idx++) {
      if (best_abs[idx] == FLT_MAX || best_abs[idx] > hw_px) continue;

      ptrdiff_t ci = idx % dims[0];
      ptrdiff_t cj = idx / dims[0];
      int is_boundary = 0;

      for (int n = 0; n < 4; n++) {
        ptrdiff_t ni = ci + dn_i[n];
        ptrdiff_t nj = cj + dn_j[n];
        if (ni < 0 || ni >= dims[0] || nj < 0 || nj >= dims[1]) {
          is_boundary = 1;
          break;
        }
        ptrdiff_t nidx = nj * dims[0] + ni;
        if (best_abs[nidx] == FLT_MAX || best_abs[nidx] > hw_px) {
          is_boundary = 1;
          break;
        }
      }
      if (!is_boundary) continue;

      if (signed_dist_px[idx] >= 0.0f) {
        if (n_pos >= pos_cap) {
          pos_cap *= 2;
          pos_seeds =
              (ptrdiff_t *)realloc(pos_seeds, pos_cap * sizeof(ptrdiff_t));
        }
        pos_seeds[n_pos++] = idx;
      }
      if (signed_dist_px[idx] <= 0.0f) {
        if (n_neg >= neg_cap) {
          neg_cap *= 2;
          neg_seeds =
              (ptrdiff_t *)realloc(neg_seeds, neg_cap * sizeof(ptrdiff_t));
        }
        neg_seeds[n_neg++] = idx;
      }
    }

    float *dist_pos = (float *)malloc(total * sizeof(float));
    float *dist_neg = (float *)malloc(total * sizeof(float));
    boundary_dijkstra(dist_pos, best_abs, pos_seeds, n_pos, dims, hw_px);
    boundary_dijkstra(dist_neg, best_abs, neg_seeds, n_neg, dims, hw_px);
    free(pos_seeds);
    free(neg_seeds);

    // Compute dfb = min(dist_pos, dist_neg)
    float *dfb = dist_from_boundary;
    int free_dfb = 0;
    if (dfb == NULL) {
      dfb = (float *)malloc(total * sizeof(float));
      free_dfb = 1;
    }
    for (ptrdiff_t i = 0; i < total; i++) {
      float dp = dist_pos[i], dn = dist_neg[i];
      float mn = dp < dn ? dp : dn;
      dfb[i] = (mn == FLT_MAX) ? NAN : mn * cellsize;
    }

    // Centre line: Voronoi boundary between positive-side and negative-side
    // boundary wavefronts, 1 pixel wide.
    //
    // Only keep pixels with my_sign >= 0 (equidistant or closer to positive
    // boundary) that have a -1 neighbour — this gives exactly one side of
    // the boundary so the line is 1 pixel wide.
    //
    // End-cap clipping: pixels whose nearest-segment projection falls outside
    // [0,1] on the first or last segment are beyond the track endpoints
    // (inside the semicircular cap) and are excluded.
    if (want_centre) {
      for (ptrdiff_t idx = 0; idx < total; idx++) {
        if (best_abs[idx] == FLT_MAX || best_abs[idx] > hw_px) continue;
        float dp = dist_pos[idx], dn = dist_neg[idx];
        if (dp == FLT_MAX || dn == FLT_MAX) continue;

        int my_sign = (dp < dn) ? 1 : (dp > dn) ? -1 : 0;

        // Keep only one side of the boundary (my_sign >= 0) → 1 pixel wide
        if (my_sign < 0) continue;

        ptrdiff_t ci = idx % dims[0];
        ptrdiff_t cj = idx / dims[0];

        // Equidistant pixels are always on centreline; my_sign==1 needs a -1
        // neighbour
        int on_cl = (my_sign == 0);
        if (!on_cl) {
          for (int n = 0; n < 8; n++) {
            ptrdiff_t ni = ci + k_di8[n];
            ptrdiff_t nj = cj + k_dj8[n];
            if (ni < 0 || ni >= dims[0] || nj < 0 || nj >= dims[1]) continue;
            ptrdiff_t nidx = nj * dims[0] + ni;
            if (best_abs[nidx] == FLT_MAX || best_abs[nidx] > hw_px) continue;
            float ndp = dist_pos[nidx], ndn = dist_neg[nidx];
            if (ndp == FLT_MAX || ndn == FLT_MAX) continue;
            if (ndp > ndn) {
              on_cl = 1;
              break;
            }  // nb_sign == -1
          }
        }
        if (!on_cl) continue;

        // End-cap clipping: recompute unclamped projection t onto nearest
        // segment
        ptrdiff_t seg = nseg[idx];
        float si = track_i[seg + 1] - track_i[seg];
        float sj = track_j[seg + 1] - track_j[seg];
        float len2 = si * si + sj * sj;
        if (len2 > 0.0f) {
          float t =
              ((ci - track_i[seg]) * si + (cj - track_j[seg]) * sj) / len2;
          if (seg == 0 && t < 0.0f) continue;
          if (seg == n_track_points - 2 && t > 1.0f) continue;
        }

        centre_line_i[n_centre] = (float)ci;
        centre_line_j[n_centre] = (float)cj;
        centre_width[n_centre] = (dp + dn) * cellsize;
        n_centre++;
      }
    }

    // Post-process: remove staircase elbows by working on the ORDERED path.
    // Previous approach checked D8-neighbour count in the pixel grid — this
    // fails when the Voronoi boundary has nc>2 at corners (other on_cl pixels
    // happen to be geometrically adjacent but are not path-neighbours).
    //
    // Instead: traverse the path in order, then in the ordered sequence a
    // pixel path[i] is an elbow iff path[i-1] and path[i+1] are D8-adjacent.
    // After removing path[i], skip path[i+1] to prevent adjacent removals
    // (which would disconnect the path).  Repeat until stable.
    if (want_centre && n_centre > 2) {
      // Cardinal and diagonal traversal offsets for ordered path extraction.
      // MUST check cardinal neighbours before diagonal: at a D4 corner the
      // diagonal next-next pixel is geometrically adjacent to the corner
      // pixel, so a diagonal-first ordering would skip the cardinal next
      // pixel and produce a wrong traversal sequence.
      static const int c_di[4] = {-1, 0, 0, 1};
      static const int c_dj[4] = {0, -1, 1, 0};
      static const int d_di[4] = {-1, -1, 1, 1};
      static const int d_dj[4] = {-1, 1, -1, 1};

      uint8_t *on_cl = (uint8_t *)calloc(total, 1);
      uint8_t *vis = (uint8_t *)calloc(total, 1);
      ptrdiff_t *path = (ptrdiff_t *)malloc(n_centre * sizeof(ptrdiff_t));

      if (on_cl && vis && path) {
        for (ptrdiff_t k = 0; k < n_centre; k++)
          on_cl[(ptrdiff_t)centre_line_j[k] * dims[0] +
                (ptrdiff_t)centre_line_i[k]] = 1;

        int changed = 1;
        while (changed) {
          changed = 0;

          // Find an endpoint (nc==1); fall back to first on_cl pixel
          ptrdiff_t start = -1;
          for (ptrdiff_t k = 0; k < n_centre && start < 0; k++) {
            ptrdiff_t ci = (ptrdiff_t)centre_line_i[k];
            ptrdiff_t cj = (ptrdiff_t)centre_line_j[k];
            if (!on_cl[cj * dims[0] + ci]) continue;
            int nc = 0;
            for (int n = 0; n < 8; n++) {
              ptrdiff_t ni = ci + k_di8[n], nj = cj + k_dj8[n];
              if (ni < 0 || ni >= dims[0] || nj < 0 || nj >= dims[1]) continue;
              if (on_cl[nj * dims[0] + ni]) nc++;
            }
            if (nc == 1) start = cj * dims[0] + ci;
          }
          if (start < 0)
            start = (ptrdiff_t)centre_line_j[0] * dims[0] +
                    (ptrdiff_t)centre_line_i[0];

          // Traverse path in order.
          memset(vis, 0, (size_t)total);
          ptrdiff_t n_path = 0;
          ptrdiff_t cur = start;
          while (cur >= 0 && n_path < n_centre) {
            path[n_path++] = cur;
            vis[cur] = 1;
            ptrdiff_t ci = cur % dims[0], cj = cur / dims[0];
            ptrdiff_t nxt = -1;
            for (int n = 0; n < 4 && nxt < 0; n++) {
              ptrdiff_t ni = ci + c_di[n], nj = cj + c_dj[n];
              if (ni < 0 || ni >= dims[0] || nj < 0 || nj >= dims[1]) continue;
              ptrdiff_t nidx = nj * dims[0] + ni;
              if (on_cl[nidx] && !vis[nidx]) nxt = nidx;
            }
            for (int n = 0; n < 4 && nxt < 0; n++) {
              ptrdiff_t ni = ci + d_di[n], nj = cj + d_dj[n];
              if (ni < 0 || ni >= dims[0] || nj < 0 || nj >= dims[1]) continue;
              ptrdiff_t nidx = nj * dims[0] + ni;
              if (on_cl[nidx] && !vis[nidx]) nxt = nidx;
            }
            cur = nxt;
          }

          // Remove elbows in the ordered sequence.
          // path[i] is an elbow if its ordered predecessor path[i-1] and
          // successor path[i+1] are D8-adjacent.  Skip path[i+1] after
          // removing path[i] to prevent removing two adjacent pixels.
          for (ptrdiff_t i = 1; i < n_path - 1; i++) {
            ptrdiff_t ai = path[i - 1] % dims[0], aj = path[i - 1] / dims[0];
            ptrdiff_t bi = path[i + 1] % dims[0], bj = path[i + 1] / dims[0];
            ptrdiff_t di = ai - bi, dj = aj - bj;
            if (di >= -1 && di <= 1 && dj >= -1 && dj <= 1) {
              on_cl[path[i]] = 0;
              changed = 1;
              i++;  // skip successor — it now has a different predecessor
            }
          }
        }
      }

      free(vis);
      free(path);

      if (on_cl) {
        ptrdiff_t new_n = 0;
        for (ptrdiff_t k = 0; k < n_centre; k++) {
          ptrdiff_t ci = (ptrdiff_t)centre_line_i[k];
          ptrdiff_t cj = (ptrdiff_t)centre_line_j[k];
          if (on_cl[cj * dims[0] + ci]) {
            centre_line_i[new_n] = centre_line_i[k];
            centre_line_j[new_n] = centre_line_j[k];
            centre_width[new_n] = centre_width[k];
            new_n++;
          }
        }
        n_centre = new_n;
        free(on_cl);
      }
    }

    free(dist_pos);
    free(dist_neg);
    if (free_dfb) free(dfb);
  }

  free(best_abs);
  free(signed_dist_px);
  if (free_nseg) free(nseg);

  return n_centre;
}

// ============================================================================
// Full (unbounded) distance map — expands to all reachable pixels
// ============================================================================
// Optional dem/mask restrict expansion to non-NaN / non-zero pixels.

TOPOTOOLBOX_API
void swath_compute_full_distance_map(
    float *restrict distance, ptrdiff_t *restrict nearest_segment,
    const float *restrict track_i, const float *restrict track_j,
    ptrdiff_t n_track_points, ptrdiff_t dims[2], float cellsize,
    const float *dem, const int8_t *mask, int compute_signed) {
  if (n_track_points < 2) return;

  ptrdiff_t total = dims[0] * dims[1];

  ptrdiff_t *nseg = NULL;
  int free_nseg = 0;
  if (nearest_segment != NULL) {
    nseg = nearest_segment;
  } else {
    nseg = (ptrdiff_t *)malloc(total * sizeof(ptrdiff_t));
    free_nseg = 1;
  }

  float *best_abs = (float *)malloc(total * sizeof(float));
  float *signed_dist_px = NULL;
  if (compute_signed) {
    signed_dist_px = (float *)malloc(total * sizeof(float));
  }

  frontier_distance_map(best_abs, signed_dist_px, nseg, track_i, track_j,
                        n_track_points, dims, FLT_MAX, dem, mask);

  for (ptrdiff_t i = 0; i < total; i++) {
    if (best_abs[i] == FLT_MAX) {
      distance[i] = NAN;
    } else if (compute_signed && signed_dist_px) {
      distance[i] = signed_dist_px[i] * cellsize;
    } else {
      distance[i] = best_abs[i] * cellsize;
    }
  }

  free(best_abs);
  if (signed_dist_px) free(signed_dist_px);
  if (free_nseg) free(nseg);
}

// ============================================================================
// Transverse swath — DEM pixels binned by signed perpendicular distance
// ============================================================================
// Requires a pre-computed signed distance map (swath_compute_distance_map).
// If normalize != 0, elevations are shifted so the track-centre bin is at zero.

TOPOTOOLBOX_API
void swath_transverse(float *restrict bin_distances, float *restrict bin_means,
                      float *restrict bin_stddevs, float *restrict bin_mins,
                      float *restrict bin_maxs, ptrdiff_t *restrict bin_counts,
                      float *restrict bin_medians, float *restrict bin_q1,
                      float *restrict bin_q3,
                      const int *restrict percentile_list,
                      ptrdiff_t n_percentiles, float *restrict bin_percentiles,
                      const float *restrict dem,
                      const float *restrict distance_from_track,
                      ptrdiff_t dims[2], float half_width, float bin_resolution,
                      ptrdiff_t n_bins, int normalize) {
  if (n_bins <= 0) return;

  ptrdiff_t total_pixels = dims[0] * dims[1];

  int compute_percentiles =
      (bin_medians != NULL || bin_q1 != NULL || bin_q3 != NULL ||
       (percentile_list != NULL && n_percentiles > 0 &&
        bin_percentiles != NULL));

  // Allocate accumulators
  swath_stats_accumulator *accumulators = (swath_stats_accumulator *)calloc(
      n_bins, sizeof(swath_stats_accumulator));
  for (ptrdiff_t b = 0; b < n_bins; b++) accumulator_init(&accumulators[b]);

  percentile_accumulator *p_accumulators = NULL;
  if (compute_percentiles) {
    p_accumulators = (percentile_accumulator *)calloc(
        n_bins, sizeof(percentile_accumulator));
    for (ptrdiff_t b = 0; b < n_bins; b++)
      percentile_accumulator_init(&p_accumulators[b], 64);
  }

  // Normalization: mean elevation near track center
  float reference_elevation = 0.0f;
  if (normalize) {
    swath_stats_accumulator track_acc;
    accumulator_init(&track_acc);
    for (ptrdiff_t idx = 0; idx < total_pixels; idx++) {
      if (isnan(distance_from_track[idx])) continue;
      if (fabsf(distance_from_track[idx]) <= bin_resolution) {
        accumulator_add(&track_acc, dem[idx]);
      }
    }
    reference_elevation = accumulator_mean(&track_acc);
  }

  // Bin pixels by signed perpendicular distance (input is in meters)
  for (ptrdiff_t idx = 0; idx < total_pixels; idx++) {
    float dist_m = distance_from_track[idx];
    if (isnan(dist_m)) continue;
    if (fabsf(dist_m) > half_width) continue;
    if (isnan(dem[idx])) continue;

    ptrdiff_t bin_idx =
        (ptrdiff_t)((dist_m + half_width) / bin_resolution + 0.5f);
    if (bin_idx < 0) bin_idx = 0;
    if (bin_idx >= n_bins) bin_idx = n_bins - 1;

    float value = dem[idx];
    if (normalize) value -= reference_elevation;

    accumulator_add(&accumulators[bin_idx], value);
    if (compute_percentiles) {
      percentile_accumulator_add(&p_accumulators[bin_idx], value);
    }
  }

  // Percentiles
  if (compute_percentiles) {
    for (ptrdiff_t b = 0; b < n_bins; b++) {
      percentile_accumulator_sort(&p_accumulators[b]);

      if (bin_medians != NULL) {
        bin_medians[b] = percentile_accumulator_get(&p_accumulators[b], 50.0f);
        if (normalize && !isnan(bin_medians[b]))
          bin_medians[b] += reference_elevation;
      }
      if (bin_q1 != NULL) {
        bin_q1[b] = percentile_accumulator_get(&p_accumulators[b], 25.0f);
        if (normalize && !isnan(bin_q1[b])) bin_q1[b] += reference_elevation;
      }
      if (bin_q3 != NULL) {
        bin_q3[b] = percentile_accumulator_get(&p_accumulators[b], 75.0f);
        if (normalize && !isnan(bin_q3[b])) bin_q3[b] += reference_elevation;
      }
      if (percentile_list != NULL && n_percentiles > 0 &&
          bin_percentiles != NULL) {
        for (ptrdiff_t p = 0; p < n_percentiles; p++) {
          float pval = percentile_accumulator_get(&p_accumulators[b],
                                                  (float)percentile_list[p]);
          if (normalize && !isnan(pval)) pval += reference_elevation;
          bin_percentiles[b * n_percentiles + p] = pval;
        }
      }
    }
  }

  // Finalize bin statistics
  for (ptrdiff_t b = 0; b < n_bins; b++) {
    bin_distances[b] = -half_width + b * bin_resolution;
    float mean = accumulator_mean(&accumulators[b]);
    float stddev = accumulator_stddev(&accumulators[b]);

    if (normalize && !isnan(mean)) mean += reference_elevation;

    bin_means[b] = mean;
    bin_stddevs[b] = stddev;
    bin_mins[b] = accumulators[b].count > 0 ? accumulators[b].min_val : NAN;
    bin_maxs[b] = accumulators[b].count > 0 ? accumulators[b].max_val : NAN;
    bin_counts[b] = accumulators[b].count;

    if (normalize && bin_counts[b] > 0) {
      bin_mins[b] += reference_elevation;
      bin_maxs[b] += reference_elevation;
    }
  }

  // Cleanup
  free(accumulators);
  if (compute_percentiles) {
    for (ptrdiff_t b = 0; b < n_bins; b++)
      percentile_accumulator_free(&p_accumulators[b]);
    free(p_accumulators);
  }
}

// ============================================================================
// Local tangent via PCA — used by longitudinal and windowed swath functions
// ============================================================================
//
// Fits a 2×2 covariance matrix over a window of n_points_regression
// neighbouring track points centred on pt.  The principal eigenvector gives
// the local tangent direction; orientation is fixed to agree with the
// forward-difference at pt.  Falls back to the endpoint direction when
// the covariance is degenerate (e.g. integer-quantised track coordinates).

static void compute_local_tangent(const float *track_i, const float *track_j,
                                  ptrdiff_t n_track_points, ptrdiff_t pt,
                                  ptrdiff_t n_points_regression, float *ti,
                                  float *tj) {
  // Window centred on pt, clamped to track bounds, at least 2 points
  ptrdiff_t half_n = n_points_regression / 2;
  if (half_n < 1) half_n = 1;

  ptrdiff_t win_lo = pt - half_n;
  ptrdiff_t win_hi = pt + half_n;
  if (win_lo < 0) win_lo = 0;
  if (win_hi >= n_track_points) win_hi = n_track_points - 1;
  if (win_hi - win_lo < 1) {
    if (pt < n_track_points - 1) {
      win_lo = pt;
      win_hi = pt + 1;
    } else {
      win_lo = pt - 1;
      win_hi = pt;
    }
  }

  ptrdiff_t n = win_hi - win_lo + 1;

  // Mean
  float mean_i = 0.0f, mean_j = 0.0f;
  for (ptrdiff_t k = win_lo; k <= win_hi; k++) {
    mean_i += track_i[k];
    mean_j += track_j[k];
  }
  mean_i /= n;
  mean_j /= n;

  // 2x2 covariance matrix
  float cov_ii = 0.0f, cov_ij = 0.0f, cov_jj = 0.0f;
  for (ptrdiff_t k = win_lo; k <= win_hi; k++) {
    float di = track_i[k] - mean_i;
    float dj = track_j[k] - mean_j;
    cov_ii += di * di;
    cov_ij += di * dj;
    cov_jj += dj * dj;
  }

  float tang_i, tang_j;

  // Principal eigenvector of the 2x2 covariance matrix [[a, b], [b, c]]
  // via direct algebraic formula (no atan2/trig, avoids quadrant issues).
  //
  // λ_max = (a + c + D) / 2  where D = sqrt((a-c)² + 4b²)
  // Eigenvector for λ_max: proportional to (2b, c - a + D)
  //   or equivalently (a - c + D, 2b) when c < a.
  //
  // Falls back to endpoint direction when the covariance is degenerate
  // (quantised/integer track coordinates with < 1 pixel cross-axis variation).

  float diff = cov_ii - cov_jj;
  float D = sqrtf(diff * diff + 4.0f * cov_ij * cov_ij);
  float vi, vj;

  if (cov_ii >= cov_jj) {
    // Use form (a - c + D, 2b) — stable when cov_ii >= cov_jj
    vi = diff + D;
    vj = 2.0f * cov_ij;
  } else {
    // Use form (2b, c - a + D) — stable when cov_jj > cov_ii
    vi = 2.0f * cov_ij;
    vj = -diff + D;
  }

  float vlen = sqrtf(vi * vi + vj * vj);

  if (vlen > 1e-10f) {
    tang_i = vi / vlen;
    tang_j = vj / vlen;
  } else {
    // Degenerate: cov_ij ≈ 0 and cov_ii ≈ cov_jj (or both ≈ 0).
    // Fall back to endpoint direction, expanding window until resolved.
    tang_i = 0.0f;
    tang_j = 0.0f;
    ptrdiff_t lo = win_lo, hi = win_hi;
    for (;;) {
      float di = track_i[hi] - track_i[lo];
      float dj = track_j[hi] - track_j[lo];
      float len = sqrtf(di * di + dj * dj);

      if (len > 0.5f) {
        float minor = (fabsf(di) < fabsf(dj)) ? fabsf(di) : fabsf(dj);
        if (minor >= 2.0f || (lo == 0 && hi == n_track_points - 1)) {
          tang_i = di / len;
          tang_j = dj / len;
          break;
        }
      }

      if (lo == 0 && hi == n_track_points - 1) {
        if (len > 0.0f) {
          tang_i = di / len;
          tang_j = dj / len;
        } else {
          tang_i = 1.0f;
          tang_j = 0.0f;
        }
        break;
      }

      if (lo > 0) lo--;
      if (hi < n_track_points - 1) hi++;
    }
  }

  // Orient consistently: dot with forward difference at pt
  float fwd_i, fwd_j;
  if (pt > 0 && pt < n_track_points - 1) {
    fwd_i = track_i[pt + 1] - track_i[pt - 1];
    fwd_j = track_j[pt + 1] - track_j[pt - 1];
  } else if (pt == 0) {
    fwd_i = track_i[1] - track_i[0];
    fwd_j = track_j[1] - track_j[0];
  } else {
    fwd_i = track_i[pt] - track_i[pt - 1];
    fwd_j = track_j[pt] - track_j[pt - 1];
  }
  if (tang_i * fwd_i + tang_j * fwd_j < 0.0f) {
    tang_i = -tang_i;
    tang_j = -tang_j;
  }

  *ti = tang_i;
  *tj = tang_j;
}

// ============================================================================
// Longitudinal swath — per-track-point statistics (two modes)
// ============================================================================
//
// Case 1 (binning_distance <= 0):
//   For each track point, walk a Bresenham line along the orthogonal direction
//   (±half_width pixels) and accumulate pixels whose distance_from_track is
//   within half_width metres.
//
// Case 2 (binning_distance > 0):
//   Assigns every swath pixel to its nearest track point via Meijster EDT
//   (compute_track_assignment), then applies a sliding window of width
//   ±binning_distance metres along the cumulative track distance, using
//   Lemire monotonic deques for O(1) amortised min/max per output point.

// ============================================================================
// Bresenham rasterization — D8 (diagonal allowed) and D4 (cardinal only)
// ============================================================================
//
// Both variants produce pixel sequences from (i0,j0) to (i1,j1).
// D4 splits diagonal steps into two cardinal steps, used by
// sample_points_between_refs when strict 4-connectivity is needed.

// D8: each step moves exactly 1 pixel (cardinal or diagonal).
static ptrdiff_t bresenham_d8(ptrdiff_t *out_i, ptrdiff_t *out_j, ptrdiff_t i0,
                              ptrdiff_t j0, ptrdiff_t i1, ptrdiff_t j1,
                              int skip_first) {
  ptrdiff_t di = i1 - i0;
  ptrdiff_t dj = j1 - j0;
  ptrdiff_t si = di > 0 ? 1 : (di < 0 ? -1 : 0);
  ptrdiff_t sj = dj > 0 ? 1 : (dj < 0 ? -1 : 0);
  ptrdiff_t adi = di < 0 ? -di : di;
  ptrdiff_t adj = dj < 0 ? -dj : dj;

  ptrdiff_t count = 0;
  ptrdiff_t ci = i0, cj = j0;

  if (!skip_first) {
    out_i[count] = ci;
    out_j[count] = cj;
    count++;
  }

  if (adi >= adj) {
    // i is the driving axis
    ptrdiff_t err = adi / 2;
    for (ptrdiff_t step = 0; step < adi; step++) {
      err -= adj;
      if (err < 0) {
        cj += sj;
        err += adi;
      }
      ci += si;
      out_i[count] = ci;
      out_j[count] = cj;
      count++;
    }
  } else {
    // j is the driving axis
    ptrdiff_t err = adj / 2;
    for (ptrdiff_t step = 0; step < adj; step++) {
      err -= adi;
      if (err < 0) {
        ci += si;
        err += adj;
      }
      cj += sj;
      out_i[count] = ci;
      out_j[count] = cj;
      count++;
    }
  }

  return count;
}

// D4 Bresenham: only cardinal moves. When the standard algorithm would step
// diagonally, two cardinal steps are emitted instead.
static ptrdiff_t bresenham_d4(ptrdiff_t *out_i, ptrdiff_t *out_j, ptrdiff_t i0,
                              ptrdiff_t j0, ptrdiff_t i1, ptrdiff_t j1,
                              int skip_first) {
  ptrdiff_t di = i1 - i0;
  ptrdiff_t dj = j1 - j0;
  ptrdiff_t si = di > 0 ? 1 : (di < 0 ? -1 : 0);
  ptrdiff_t sj = dj > 0 ? 1 : (dj < 0 ? -1 : 0);
  ptrdiff_t adi = di < 0 ? -di : di;
  ptrdiff_t adj = dj < 0 ? -dj : dj;

  ptrdiff_t count = 0;
  ptrdiff_t ci = i0, cj = j0;

  if (!skip_first) {
    out_i[count] = ci;
    out_j[count] = cj;
    count++;
  }

  if (adi >= adj) {
    // i is the driving axis
    ptrdiff_t err = adi / 2;
    for (ptrdiff_t step = 0; step < adi; step++) {
      err -= adj;
      if (err < 0) {
        // Would be diagonal: emit two cardinal steps
        // Order based on error: step in j first if error is more negative
        cj += sj;
        out_i[count] = ci;
        out_j[count] = cj;
        count++;
        err += adi;
      }
      ci += si;
      out_i[count] = ci;
      out_j[count] = cj;
      count++;
    }
  } else {
    // j is the driving axis
    ptrdiff_t err = adj / 2;
    for (ptrdiff_t step = 0; step < adj; step++) {
      err -= adi;
      if (err < 0) {
        // Would be diagonal: emit two cardinal steps
        ci += si;
        out_i[count] = ci;
        out_j[count] = cj;
        count++;
        err += adj;
      }
      cj += sj;
      out_i[count] = ci;
      out_j[count] = cj;
      count++;
    }
  }

  return count;
}

// ============================================================================
// Meijster EDT nearest-track-point labeling (internal, Case 2 only)
// ============================================================================
//
// Two-pass algorithm: column scan builds row distances, row parabola envelope
// resolves the 2-D nearest-seed per pixel.
//
// Seeds: one per track point (rounded position) or one per rasterized segment
// pixel (use_segment_seeds != 0, more accurate for coarsely sampled tracks).
//
// Returns a malloc'd array [dims[0]*dims[1]]; caller frees.
// Value = nearest track-point index, or -1 for pixels outside the mask.
static ptrdiff_t *compute_track_assignment(const float *track_i,
                                           const float *track_j,
                                           ptrdiff_t n_track_points,
                                           const float *distance_from_track,
                                           ptrdiff_t dims[2], float half_width,
                                           int use_segment_seeds) {
  ptrdiff_t nrows = dims[0], ncols = dims[1];
  ptrdiff_t total = nrows * ncols;
  ptrdiff_t INF_VAL = nrows + ncols + 2;

  // g[j*nrows+i]  = row distance from row i to nearest seed in column j
  // lg[j*nrows+i] = track point index of that nearest seed (-1 = none)
  ptrdiff_t *g = (ptrdiff_t *)malloc(total * sizeof(ptrdiff_t));
  ptrdiff_t *lg = (ptrdiff_t *)malloc(total * sizeof(ptrdiff_t));
  for (ptrdiff_t idx = 0; idx < total; idx++) {
    g[idx] = INF_VAL;
    lg[idx] = -1;
  }

  if (use_segment_seeds) {
    // Rasterize each segment at ~0.5px spacing. For each rasterized pixel,
    // assign the nearest track point by Euclidean distance.
    for (ptrdiff_t k = 0; k < n_track_points - 1; k++) {
      float di = track_i[k + 1] - track_i[k];
      float dj = track_j[k + 1] - track_j[k];
      float seg_len = sqrtf(di * di + dj * dj);
      ptrdiff_t n_steps = (ptrdiff_t)(seg_len * 2.0f) + 1;

      for (ptrdiff_t s = 0; s <= n_steps; s++) {
        float t = n_steps > 0 ? (float)s / (float)n_steps : 0.0f;
        ptrdiff_t si = (ptrdiff_t)(track_i[k] + t * di + 0.5f);
        ptrdiff_t sj = (ptrdiff_t)(track_j[k] + t * dj + 0.5f);
        if (si < 0 || si >= nrows || sj < 0 || sj >= ncols) continue;
        ptrdiff_t sidx = sj * nrows + si;

        float min_d2 = FLT_MAX;
        ptrdiff_t best_pt = -1;
        for (ptrdiff_t p = 0; p < n_track_points; p++) {
          float tdi = (float)si - track_i[p];
          float tdj = (float)sj - track_j[p];
          float d2 = tdi * tdi + tdj * tdj;
          if (d2 < min_d2) {
            min_d2 = d2;
            best_pt = p;
          }
        }

        if (lg[sidx] < 0) {
          g[sidx] = 0;
          lg[sidx] = best_pt;
        } else {
          float tdi = (float)si - track_i[lg[sidx]];
          float tdj = (float)sj - track_j[lg[sidx]];
          if (min_d2 < tdi * tdi + tdj * tdj) {
            g[sidx] = 0;
            lg[sidx] = best_pt;
          }
        }
      }
    }
  } else {
    // One seed per track point at its rounded grid position.
    for (ptrdiff_t k = 0; k < n_track_points; k++) {
      ptrdiff_t si = (ptrdiff_t)(track_i[k] + 0.5f);
      ptrdiff_t sj = (ptrdiff_t)(track_j[k] + 0.5f);
      if (si < 0 || si >= nrows || sj < 0 || sj >= ncols) continue;
      ptrdiff_t sidx = sj * nrows + si;
      if (lg[sidx] < 0) {
        g[sidx] = 0;
        lg[sidx] = k;
      } else {
        float aki = (float)si - track_i[k], akj = (float)sj - track_j[k];
        float api = (float)si - track_i[lg[sidx]],
              apj = (float)sj - track_j[lg[sidx]];
        if (aki * aki + akj * akj < api * api + apj * apj) lg[sidx] = k;
      }
    }
  }

  // Phase 1: column scan — for each column j, propagate nearest-seed row
  // distance up and down so g[j*nrows+i] = distance to nearest seed in col j.
  for (ptrdiff_t j = 0; j < ncols; j++) {
    for (ptrdiff_t i = 1; i < nrows; i++) {
      ptrdiff_t cur = j * nrows + i, prev = j * nrows + (i - 1);
      if (g[prev] < INF_VAL && g[prev] + 1 < g[cur]) {
        g[cur] = g[prev] + 1;
        lg[cur] = lg[prev];
      }
    }
    for (ptrdiff_t i = nrows - 2; i >= 0; i--) {
      ptrdiff_t cur = j * nrows + i, next = j * nrows + (i + 1);
      if (g[next] < INF_VAL && g[next] + 1 < g[cur]) {
        g[cur] = g[next] + 1;
        lg[cur] = lg[next];
      }
    }
  }

  // Phase 2: row parabola envelope (Meijster).
  // For each row i, find which column j's seed is the globally nearest
  // using the parabola Sep function: first col where u beats s is
  //   w = floor((u²+gu² - s²-gs²) / (2*(u-s))) + 1
  ptrdiff_t *assignment = (ptrdiff_t *)malloc(total * sizeof(ptrdiff_t));
  for (ptrdiff_t idx = 0; idx < total; idx++) assignment[idx] = -1;

  ptrdiff_t *stk_s = (ptrdiff_t *)malloc(ncols * sizeof(ptrdiff_t));
  ptrdiff_t *stk_t = (ptrdiff_t *)malloc(ncols * sizeof(ptrdiff_t));
  ptrdiff_t *stk_l = (ptrdiff_t *)malloc(ncols * sizeof(ptrdiff_t));

  for (ptrdiff_t i = 0; i < nrows; i++) {
    ptrdiff_t q = -1;

    // Forward scan: build parabola stack
    for (ptrdiff_t u = 0; u < ncols; u++) {
      ptrdiff_t gu = g[u * nrows + i];
      if (gu >= INF_VAL) continue;

      // Pop while u dominates top at top's start column
      while (q >= 0) {
        ptrdiff_t s = stk_s[q], gs = g[s * nrows + i];
        double sep = ((double)(u * u) - (double)(s * s) + (double)(gu * gu) -
                      (double)(gs * gs)) /
                     (2.0 * (double)(u - s));
        ptrdiff_t w = (ptrdiff_t)floor(sep) + 1;
        if (w > stk_t[q])
          break;  // u only takes over after top's range: keep top
        q--;
      }

      // Compute where u starts being the closest
      ptrdiff_t w_start;
      if (q < 0) {
        w_start = 0;
      } else {
        ptrdiff_t s = stk_s[q], gs = g[s * nrows + i];
        double sep = ((double)(u * u) - (double)(s * s) + (double)(gu * gu) -
                      (double)(gs * gs)) /
                     (2.0 * (double)(u - s));
        w_start = (ptrdiff_t)floor(sep) + 1;
      }

      if (w_start < ncols) {
        q++;
        stk_s[q] = u;
        stk_t[q] = w_start;
        stk_l[q] = lg[u * nrows + i];
      }
    }

    if (q < 0) continue;

    // Backward scan: assign each masked pixel to its nearest track point
    for (ptrdiff_t u = ncols - 1; u >= 0; u--) {
      while (q > 0 && u < stk_t[q]) q--;
      ptrdiff_t idx = u * nrows + i;
      float d = distance_from_track[idx];
      if (!isnan(d) && fabsf(d) <= half_width) assignment[idx] = stk_l[q];
    }
  }

  free(g);
  free(lg);
  free(stk_s);
  free(stk_t);
  free(stk_l);
  return assignment;
}

TOPOTOOLBOX_API
ptrdiff_t swath_longitudinal(
    float *restrict point_means, float *restrict point_stddevs,
    float *restrict point_mins, float *restrict point_maxs,
    ptrdiff_t *restrict point_counts, float *restrict point_medians,
    float *restrict point_q1, float *restrict point_q3,
    const int *restrict percentile_list, ptrdiff_t n_percentiles,
    float *restrict point_percentiles, const float *restrict dem,
    const float *restrict track_i, const float *restrict track_j,
    ptrdiff_t n_track_points, const float *restrict distance_from_track,
    ptrdiff_t dims[2], float cellsize, float half_width, float binning_distance,
    ptrdiff_t n_points_regression, ptrdiff_t use_segment_seeds, ptrdiff_t skip,
    float *result_track_i, float *result_track_j) {
  if (n_track_points < 2) return 0;
  ptrdiff_t effective_skip = (skip < 1) ? 1 : skip;

  float hw_px = half_width / cellsize;

  int compute_percentiles =
      (point_medians != NULL || point_q1 != NULL || point_q3 != NULL ||
       (percentile_list != NULL && n_percentiles > 0 &&
        point_percentiles != NULL));

  // Cumulative along-track distance (meters) — needed for binning_distance > 0
  float *cum_dist = NULL;
  if (binning_distance > 0.0f) {
    cum_dist = (float *)malloc(n_track_points * sizeof(float));
    cum_dist[0] = 0.0f;
    for (ptrdiff_t k = 0; k < n_track_points - 1; k++) {
      float di = track_i[k + 1] - track_i[k];
      float dj = track_j[k + 1] - track_j[k];
      cum_dist[k + 1] = cum_dist[k] + sqrtf(di * di + dj * dj) * cellsize;
    }
  }

  // Per-point accumulators — only needed for Case 1 (Bresenham cross-sections)
  swath_stats_accumulator *accumulators = NULL;
  if (binning_distance <= 0.0f) {
    accumulators = (swath_stats_accumulator *)calloc(
        n_track_points, sizeof(swath_stats_accumulator));
    for (ptrdiff_t k = 0; k < n_track_points; k++)
      accumulator_init(&accumulators[k]);
  }

  percentile_accumulator *p_accumulators = NULL;
  if (compute_percentiles) {
    p_accumulators = (percentile_accumulator *)calloc(
        n_track_points, sizeof(percentile_accumulator));
    for (ptrdiff_t k = 0; k < n_track_points; k++)
      percentile_accumulator_init(&p_accumulators[k], 64);
  }

  if (binning_distance <= 0.0f) {
    // --- Case 1: single cross-section via Bresenham (per track point) ---
    ptrdiff_t bres_cap = (ptrdiff_t)(2.0f * hw_px) + 16;
    ptrdiff_t *bres_i = (ptrdiff_t *)malloc(bres_cap * sizeof(ptrdiff_t));
    ptrdiff_t *bres_j = (ptrdiff_t *)malloc(bres_cap * sizeof(ptrdiff_t));

    for (ptrdiff_t pt = 0; pt < n_track_points; pt++) {
      if (effective_skip > 1 && pt % effective_skip != 0) continue;
      float tang_i, tang_j;
      compute_local_tangent(track_i, track_j, n_track_points, pt,
                            n_points_regression, &tang_i, &tang_j);
      float orth_i = -tang_j;
      float orth_j = tang_i;

      ptrdiff_t ai = (ptrdiff_t)(track_i[pt] - orth_i * hw_px + 0.5f);
      ptrdiff_t aj = (ptrdiff_t)(track_j[pt] - orth_j * hw_px + 0.5f);
      ptrdiff_t bi = (ptrdiff_t)(track_i[pt] + orth_i * hw_px + 0.5f);
      ptrdiff_t bj = (ptrdiff_t)(track_j[pt] + orth_j * hw_px + 0.5f);

      ptrdiff_t n_bres = bresenham_d8(bres_i, bres_j, ai, aj, bi, bj, 0);

      for (ptrdiff_t p = 0; p < n_bres; p++) {
        ptrdiff_t pi = bres_i[p], pj = bres_j[p];
        if (pi < 0 || pi >= dims[0] || pj < 0 || pj >= dims[1]) continue;
        ptrdiff_t idx = pj * dims[0] + pi;
        float dist_m = distance_from_track[idx];
        if (isnan(dist_m)) continue;
        if (fabsf(dist_m) > half_width) continue;
        if (isnan(dem[idx])) continue;

        accumulator_add(&accumulators[pt], dem[idx]);
        if (compute_percentiles)
          percentile_accumulator_add(&p_accumulators[pt], dem[idx]);
      }
    }
    free(bres_i);
    free(bres_j);
  } else {
    // --- Case 2: sliding window over pre-computed per-block stats ---
    ptrdiff_t total = dims[0] * dims[1];
    ptrdiff_t *assignment = compute_track_assignment(
        track_i, track_j, n_track_points, distance_from_track, dims, half_width,
        (int)use_segment_seeds);

    // Per-block precomputed stats
    double *blk_sum = (double *)calloc(n_track_points, sizeof(double));
    double *blk_sum2 = (double *)calloc(n_track_points, sizeof(double));
    ptrdiff_t *blk_count =
        (ptrdiff_t *)calloc(n_track_points, sizeof(ptrdiff_t));
    float *blk_min = (float *)malloc(n_track_points * sizeof(float));
    float *blk_max = (float *)malloc(n_track_points * sizeof(float));
    for (ptrdiff_t k = 0; k < n_track_points; k++) {
      blk_min[k] = FLT_MAX;
      blk_max[k] = -FLT_MAX;
    }

    // CSR pixel lists (only needed for percentiles)
    ptrdiff_t *pt_offset = NULL, *pt_pixels = NULL, *pt_fill_arr = NULL;
    if (compute_percentiles) {
      ptrdiff_t *cnt = (ptrdiff_t *)calloc(n_track_points, sizeof(ptrdiff_t));
      ptrdiff_t n_assigned = 0;
      for (ptrdiff_t idx = 0; idx < total; idx++) {
        if (assignment[idx] >= 0 && !isnan(dem[idx])) {
          cnt[assignment[idx]]++;
          n_assigned++;
        }
      }
      pt_offset = (ptrdiff_t *)malloc((n_track_points + 1) * sizeof(ptrdiff_t));
      pt_offset[0] = 0;
      for (ptrdiff_t k = 0; k < n_track_points; k++)
        pt_offset[k + 1] = pt_offset[k] + cnt[k];
      free(cnt);
      pt_pixels = (ptrdiff_t *)malloc(n_assigned * sizeof(ptrdiff_t));
      pt_fill_arr = (ptrdiff_t *)calloc(n_track_points, sizeof(ptrdiff_t));
      // Fill block stats and CSR in one pass
      for (ptrdiff_t idx = 0; idx < total; idx++) {
        ptrdiff_t a = assignment[idx];
        if (a < 0 || isnan(dem[idx])) continue;
        double v = (double)dem[idx];
        blk_sum[a] += v;
        blk_sum2[a] += v * v;
        blk_count[a]++;
        if (dem[idx] < blk_min[a]) blk_min[a] = dem[idx];
        if (dem[idx] > blk_max[a]) blk_max[a] = dem[idx];
        pt_pixels[pt_offset[a] + pt_fill_arr[a]] = idx;
        pt_fill_arr[a]++;
      }
    } else {
      // Block stats only
      for (ptrdiff_t idx = 0; idx < total; idx++) {
        ptrdiff_t a = assignment[idx];
        if (a < 0 || isnan(dem[idx])) continue;
        double v = (double)dem[idx];
        blk_sum[a] += v;
        blk_sum2[a] += v * v;
        blk_count[a]++;
        if (dem[idx] < blk_min[a]) blk_min[a] = dem[idx];
        if (dem[idx] > blk_max[a]) blk_max[a] = dem[idx];
      }
    }
    free(assignment);

    // Monotonic deques for sliding window min/max
    ptrdiff_t *min_dq = (ptrdiff_t *)malloc(n_track_points * sizeof(ptrdiff_t));
    ptrdiff_t *max_dq = (ptrdiff_t *)malloc(n_track_points * sizeof(ptrdiff_t));
    ptrdiff_t min_head = 0, min_tail = -1;
    ptrdiff_t max_head = 0, max_tail = -1;
    double win_sum = 0.0, win_sum2 = 0.0;
    ptrdiff_t win_count = 0;
    ptrdiff_t lo = 0, hi = -1;

    for (ptrdiff_t pt = 0; pt < n_track_points; pt++) {
      // Expand right: add blocks within binning_distance ahead of pt
      while (hi < n_track_points - 1 &&
             cum_dist[hi + 1] - cum_dist[pt] <= binning_distance) {
        hi++;
        win_sum += blk_sum[hi];
        win_sum2 += blk_sum2[hi];
        win_count += blk_count[hi];
        if (blk_count[hi] > 0) {
          while (min_head <= min_tail &&
                 blk_min[min_dq[min_tail]] >= blk_min[hi])
            min_tail--;
          min_dq[++min_tail] = hi;
          while (max_head <= max_tail &&
                 blk_max[max_dq[max_tail]] <= blk_max[hi])
            max_tail--;
          max_dq[++max_tail] = hi;
        }
      }
      // Shrink left: remove blocks outside binning_distance behind pt
      while (lo <= hi && cum_dist[pt] - cum_dist[lo] > binning_distance) {
        win_sum -= blk_sum[lo];
        win_sum2 -= blk_sum2[lo];
        win_count -= blk_count[lo];
        if (min_head <= min_tail && min_dq[min_head] == lo) min_head++;
        if (max_head <= max_tail && max_dq[max_head] == lo) max_head++;
        lo++;
      }
      // Write output only for kept track points
      if (effective_skip <= 1 || pt % effective_skip == 0) {
        ptrdiff_t out = pt / effective_skip;
        if (win_count > 0) {
          double mean = win_sum / (double)win_count;
          double var = win_sum2 / (double)win_count - mean * mean;
          point_counts[out] = win_count;
          point_means[out] = (float)mean;
          point_stddevs[out] = (float)sqrt(var > 0.0 ? var : 0.0);
          point_mins[out] =
              (min_head <= min_tail) ? blk_min[min_dq[min_head]] : NAN;
          point_maxs[out] =
              (max_head <= max_tail) ? blk_max[max_dq[max_head]] : NAN;
        } else {
          point_counts[out] = 0;
          point_means[out] = NAN;
          point_stddevs[out] = NAN;
          point_mins[out] = NAN;
          point_maxs[out] = NAN;
        }
      }
      // Percentile accumulation for this window
      if (compute_percentiles &&
          (effective_skip <= 1 || pt % effective_skip == 0)) {
        for (ptrdiff_t k = lo; k <= hi; k++)
          for (ptrdiff_t q = pt_offset[k]; q < pt_offset[k + 1]; q++)
            percentile_accumulator_add(&p_accumulators[pt], dem[pt_pixels[q]]);
      }
    }

    free(blk_sum);
    free(blk_sum2);
    free(blk_count);
    free(blk_min);
    free(blk_max);
    free(min_dq);
    free(max_dq);
    if (compute_percentiles) {
      free(pt_offset);
      free(pt_pixels);
      free(pt_fill_arr);
    }
  }

  // Finalize statistics
  ptrdiff_t n_result = 0;
  for (ptrdiff_t k = 0; k < n_track_points; k++) {
    if (effective_skip > 1 && k % effective_skip != 0) continue;
    ptrdiff_t out = k / effective_skip;
    n_result = out + 1;

    if (binning_distance <= 0.0f) {
      point_means[out] = accumulator_mean(&accumulators[k]);
      point_stddevs[out] = accumulator_stddev(&accumulators[k]);
      point_mins[out] =
          accumulators[k].count > 0 ? accumulators[k].min_val : NAN;
      point_maxs[out] =
          accumulators[k].count > 0 ? accumulators[k].max_val : NAN;
      point_counts[out] = accumulators[k].count;
    }

    if (compute_percentiles) {
      percentile_accumulator_sort(&p_accumulators[k]);

      if (point_medians != NULL)
        point_medians[out] =
            percentile_accumulator_get(&p_accumulators[k], 50.0f);
      if (point_q1 != NULL)
        point_q1[out] = percentile_accumulator_get(&p_accumulators[k], 25.0f);
      if (point_q3 != NULL)
        point_q3[out] = percentile_accumulator_get(&p_accumulators[k], 75.0f);

      if (percentile_list != NULL && n_percentiles > 0 &&
          point_percentiles != NULL) {
        for (ptrdiff_t p = 0; p < n_percentiles; p++) {
          point_percentiles[out * n_percentiles + p] =
              percentile_accumulator_get(&p_accumulators[k],
                                         (float)percentile_list[p]);
        }
      }
    }

    if (result_track_i) result_track_i[out] = track_i[k];
    if (result_track_j) result_track_j[out] = track_j[k];
  }

  // Cleanup
  if (cum_dist) free(cum_dist);
  if (accumulators) free(accumulators);
  if (compute_percentiles) {
    for (ptrdiff_t k = 0; k < n_track_points; k++)
      percentile_accumulator_free(&p_accumulators[k]);
    free(p_accumulators);
  }
  return n_result;
}

// ============================================================================
// Windowed longitudinal swath — oriented rectangle, all interior pixels
// ============================================================================
//
// For each (kept) track point:
//   1. Compute local PCA tangent (ti, tj); orthogonal = (-tj, ti).
//   2. Define rectangle: ±binning_distance along track, ±half_width orthogonal.
//   3. Scan the axis-aligned bounding box; accept pixels satisfying both
//      |Δ·tang| ≤ bd_px  AND  |Δ·orth| ≤ hw_px.
//   4. Pass 1: accumulate mean/std/min/max/count.
//   5. Pass 2 (if percentiles requested): 2048-bin histogram over [vmin, vmax],
//      queries resolved by rank scan.
//
// No pre-computed distance map required.  Independent of swath_longitudinal.

TOPOTOOLBOX_API
ptrdiff_t swath_longitudinal_windowed(
    float *point_means, float *point_stddevs, float *point_mins,
    float *point_maxs, ptrdiff_t *point_counts, float *point_medians,
    float *point_q1, float *point_q3, const int *percentile_list,
    ptrdiff_t n_percentiles, float *point_percentiles, const float *dem,
    const float *track_i, const float *track_j, ptrdiff_t n_track_points,
    ptrdiff_t dims[2], float cellsize, float half_width, float binning_distance,
    ptrdiff_t n_points_regression, ptrdiff_t skip, float *result_track_i,
    float *result_track_j) {
  if (n_track_points < 2) return 0;
  ptrdiff_t effective_skip = (skip < 1) ? 1 : skip;

  ptrdiff_t nrows = dims[0], ncols = dims[1];
  float hw_px = half_width / cellsize;
  float bd_px = binning_distance / cellsize;
  static const int slw_n_bins = 2048;

  int compute_percentiles =
      (point_medians != NULL || point_q1 != NULL || point_q3 != NULL ||
       (percentile_list != NULL && n_percentiles > 0 &&
        point_percentiles != NULL));

  ptrdiff_t *hist = compute_percentiles
                        ? (ptrdiff_t *)calloc(slw_n_bins, sizeof(ptrdiff_t))
                        : NULL;

  ptrdiff_t n_result = 0;
  for (ptrdiff_t pt = 0; pt < n_track_points; pt++) {
    if (effective_skip > 1 && pt % effective_skip != 0) continue;
    ptrdiff_t out = pt / effective_skip;
    float ti, tj;
    compute_local_tangent(track_i, track_j, n_track_points, pt,
                          n_points_regression, &ti, &tj);
    float oi = -tj, oj = ti;  // orthogonal direction
    float ci = track_i[pt], cj = track_j[pt];

    // Axis-aligned bounding box of the oriented rectangle
    float R = hw_px + bd_px;
    ptrdiff_t pi_lo = (ptrdiff_t)(ci - R);
    if (pi_lo < 0) pi_lo = 0;
    ptrdiff_t pi_hi = (ptrdiff_t)(ci + R) + 1;
    if (pi_hi >= nrows) pi_hi = nrows - 1;
    ptrdiff_t pj_lo = (ptrdiff_t)(cj - R);
    if (pj_lo < 0) pj_lo = 0;
    ptrdiff_t pj_hi = (ptrdiff_t)(cj + R) + 1;
    if (pj_hi >= ncols) pj_hi = ncols - 1;

    // Pass 1: mean / std / min / max / count
    double sum = 0.0, sum2 = 0.0;
    ptrdiff_t count = 0;
    float vmin = FLT_MAX, vmax = -FLT_MAX;

    for (ptrdiff_t pj = pj_lo; pj <= pj_hi; pj++) {
      for (ptrdiff_t pi = pi_lo; pi <= pi_hi; pi++) {
        float di = (float)pi - ci, dj = (float)pj - cj;
        if (fabsf(di * ti + dj * tj) > bd_px) continue;  // along-track
        if (fabsf(di * oi + dj * oj) > hw_px) continue;  // orthogonal
        float v = dem[pj * nrows + pi];
        if (isnan(v)) continue;
        sum += (double)v;
        sum2 += (double)v * (double)v;
        count++;
        if (v < vmin) vmin = v;
        if (v > vmax) vmax = v;
      }
    }

    if (count > 0) {
      double mean = sum / (double)count;
      double var = sum2 / (double)count - mean * mean;
      point_counts[out] = count;
      point_means[out] = (float)mean;
      point_stddevs[out] = (float)sqrt(var > 0.0 ? var : 0.0);
      point_mins[out] = vmin;
      point_maxs[out] = vmax;
    } else {
      point_counts[out] = 0;
      point_means[out] = NAN;
      point_stddevs[out] = NAN;
      point_mins[out] = NAN;
      point_maxs[out] = NAN;
    }

    // Pass 2 (percentiles): histogram over window pixels
    if (compute_percentiles) {
      if (count <= 0) {
        if (point_medians) point_medians[out] = NAN;
        if (point_q1) point_q1[out] = NAN;
        if (point_q3) point_q3[out] = NAN;
        if (percentile_list && n_percentiles > 0 && point_percentiles)
          for (ptrdiff_t p = 0; p < n_percentiles; p++)
            point_percentiles[out * n_percentiles + p] = NAN;
      } else if (vmax <= vmin) {
        if (point_medians) point_medians[out] = vmin;
        if (point_q1) point_q1[out] = vmin;
        if (point_q3) point_q3[out] = vmin;
        if (percentile_list && n_percentiles > 0 && point_percentiles)
          for (ptrdiff_t p = 0; p < n_percentiles; p++)
            point_percentiles[out * n_percentiles + p] = vmin;
      } else {
        float bw = (vmax - vmin) / (float)slw_n_bins;

        memset(hist, 0, slw_n_bins * sizeof(ptrdiff_t));

        for (ptrdiff_t pj = pj_lo; pj <= pj_hi; pj++) {
          for (ptrdiff_t pi = pi_lo; pi <= pi_hi; pi++) {
            float di = (float)pi - ci, dj = (float)pj - cj;
            if (fabsf(di * ti + dj * tj) > bd_px) continue;
            if (fabsf(di * oi + dj * oj) > hw_px) continue;
            float v = dem[pj * nrows + pi];
            if (isnan(v)) continue;
            int b = (int)((v - vmin) / bw);
            if (b < 0) b = 0;
            if (b >= slw_n_bins) b = slw_n_bins - 1;
            hist[b]++;
          }
        }

#define HPCT(pval, dst)                                              \
  do {                                                               \
    ptrdiff_t _t = (ptrdiff_t)ceilf((pval) / 100.0f * (float)count); \
    _t = (_t < 1) ? 1 : (_t > count ? count : _t);                   \
    ptrdiff_t _c = 0;                                                \
    (dst) = vmax;                                                    \
    for (ptrdiff_t _b = 0; _b < slw_n_bins; _b++) {                  \
      _c += hist[_b];                                                \
      if (_c >= _t) {                                                \
        (dst) = vmin + ((float)_b + 0.5f) * bw;                      \
        break;                                                       \
      }                                                              \
    }                                                                \
  } while (0)

        if (point_medians) HPCT(50.0f, point_medians[out]);
        if (point_q1) HPCT(25.0f, point_q1[out]);
        if (point_q3) HPCT(75.0f, point_q3[out]);
        if (percentile_list && n_percentiles > 0 && point_percentiles)
          for (ptrdiff_t p = 0; p < n_percentiles; p++)
            HPCT((float)percentile_list[p],
                 point_percentiles[out * n_percentiles + p]);

#undef HPCT
      }
    }

    if (result_track_i) result_track_i[out] = track_i[pt];
    if (result_track_j) result_track_j[out] = track_j[pt];
    n_result = out + 1;
  }

  if (hist) free(hist);
  return n_result;
}

// Returns integer pixel coordinates of all pixels inside the oriented rectangle
// for a single track point.  Mirrors swath_longitudinal_windowed exactly.
// Safe upper bound for output arrays: dims[0] * dims[1].
TOPOTOOLBOX_API
ptrdiff_t swath_windowed_get_point_samples(
    ptrdiff_t *pixels_i, ptrdiff_t *pixels_j, const float *track_i,
    const float *track_j, ptrdiff_t n_track_points, ptrdiff_t point_index,
    ptrdiff_t dims[2], float cellsize, float half_width, float binning_distance,
    ptrdiff_t n_points_regression) {
  if (n_track_points < 2) return 0;
  if (point_index < 0 || point_index >= n_track_points) return 0;

  ptrdiff_t nrows = dims[0], ncols = dims[1];
  float hw_px = half_width / cellsize;
  float bd_px = binning_distance / cellsize;

  float ti, tj;
  compute_local_tangent(track_i, track_j, n_track_points, point_index,
                        n_points_regression, &ti, &tj);
  float oi = -tj, oj = ti;
  float ci = track_i[point_index], cj = track_j[point_index];

  float R = hw_px + bd_px;
  ptrdiff_t pi_lo = (ptrdiff_t)(ci - R);
  if (pi_lo < 0) pi_lo = 0;
  ptrdiff_t pi_hi = (ptrdiff_t)(ci + R) + 1;
  if (pi_hi >= nrows) pi_hi = nrows - 1;
  ptrdiff_t pj_lo = (ptrdiff_t)(cj - R);
  if (pj_lo < 0) pj_lo = 0;
  ptrdiff_t pj_hi = (ptrdiff_t)(cj + R) + 1;
  if (pj_hi >= ncols) pj_hi = ncols - 1;

  ptrdiff_t n_out = 0;
  for (ptrdiff_t pj = pj_lo; pj <= pj_hi; pj++) {
    for (ptrdiff_t pi = pi_lo; pi <= pi_hi; pi++) {
      float di = (float)pi - ci, dj = (float)pj - cj;
      if (fabsf(di * ti + dj * tj) > bd_px) continue;
      if (fabsf(di * oi + dj * oj) > hw_px) continue;
      pixels_i[n_out] = pi;
      pixels_j[n_out] = pj;
      n_out++;
    }
  }
  return n_out;
}

// ============================================================================
// Per-point pixel retrieval — mirrors swath_longitudinal (vanilla)
// ============================================================================
//
// Case 1 (binning_distance <= 0): Bresenham cross-section filtered by
//   distance_from_track.
// Case 2 (binning_distance  > 0): Meijster EDT assignment; gathers all pixels
//   whose nearest track point falls within [pt_lo, pt_hi] (the
//   ±binning_distance window around point_index along the cumulative distance).

TOPOTOOLBOX_API
ptrdiff_t swath_get_point_pixels(
    ptrdiff_t *restrict pixels_i, ptrdiff_t *restrict pixels_j,
    const float *restrict track_i, const float *restrict track_j,
    ptrdiff_t n_track_points, ptrdiff_t point_index,
    const float *restrict distance_from_track, ptrdiff_t dims[2],
    float cellsize, float half_width, float binning_distance,
    ptrdiff_t n_points_regression, ptrdiff_t use_segment_seeds) {
  if (n_track_points < 2) return 0;
  if (point_index < 0 || point_index >= n_track_points) return 0;

  float hw_px = half_width / cellsize;
  ptrdiff_t n_pixels = 0;

  float tang_i, tang_j;
  compute_local_tangent(track_i, track_j, n_track_points, point_index,
                        n_points_regression, &tang_i, &tang_j);
  float orth_i = -tang_j;
  float orth_j = tang_i;

  if (binning_distance <= 0.0f) {
    // --- Case 1: single cross-section via Bresenham ---
    ptrdiff_t ai = (ptrdiff_t)(track_i[point_index] - orth_i * hw_px + 0.5f);
    ptrdiff_t aj = (ptrdiff_t)(track_j[point_index] - orth_j * hw_px + 0.5f);
    ptrdiff_t bi = (ptrdiff_t)(track_i[point_index] + orth_i * hw_px + 0.5f);
    ptrdiff_t bj = (ptrdiff_t)(track_j[point_index] + orth_j * hw_px + 0.5f);

    ptrdiff_t bres_cap = (ptrdiff_t)(2.0f * hw_px) + 16;
    ptrdiff_t *bres_i = (ptrdiff_t *)malloc(bres_cap * sizeof(ptrdiff_t));
    ptrdiff_t *bres_j = (ptrdiff_t *)malloc(bres_cap * sizeof(ptrdiff_t));

    ptrdiff_t n_bres = bresenham_d8(bres_i, bres_j, ai, aj, bi, bj, 0);

    for (ptrdiff_t p = 0; p < n_bres; p++) {
      ptrdiff_t pi = bres_i[p], pj = bres_j[p];
      if (pi < 0 || pi >= dims[0] || pj < 0 || pj >= dims[1]) continue;
      ptrdiff_t idx = pj * dims[0] + pi;
      float dist_m = distance_from_track[idx];
      if (isnan(dist_m)) continue;
      if (fabsf(dist_m) > half_width) continue;

      pixels_i[n_pixels] = pi;
      pixels_j[n_pixels] = pj;
      n_pixels++;
    }

    free(bres_i);
    free(bres_j);
  } else {
    // --- Case 2: Meijster EDT nearest-track-point assignment ---

    float *cum_dist = (float *)malloc(n_track_points * sizeof(float));
    cum_dist[0] = 0.0f;
    for (ptrdiff_t k = 0; k < n_track_points - 1; k++) {
      float di = track_i[k + 1] - track_i[k];
      float dj = track_j[k + 1] - track_j[k];
      cum_dist[k + 1] = cum_dist[k] + sqrtf(di * di + dj * dj) * cellsize;
    }

    ptrdiff_t pt_lo = point_index, pt_hi = point_index;
    while (pt_lo > 0 &&
           cum_dist[point_index] - cum_dist[pt_lo - 1] <= binning_distance)
      pt_lo--;
    while (pt_hi < n_track_points - 1 &&
           cum_dist[pt_hi + 1] - cum_dist[point_index] <= binning_distance)
      pt_hi++;

    ptrdiff_t *assignment = compute_track_assignment(
        track_i, track_j, n_track_points, distance_from_track, dims, half_width,
        (int)use_segment_seeds);

    // Gather pixels assigned to track points in [pt_lo..pt_hi]
    for (ptrdiff_t pj = 0; pj < dims[1]; pj++) {
      for (ptrdiff_t pi = 0; pi < dims[0]; pi++) {
        ptrdiff_t idx = pj * dims[0] + pi;
        ptrdiff_t a = assignment[idx];
        if (a < pt_lo || a > pt_hi) continue;

        pixels_i[n_pixels] = pi;
        pixels_j[n_pixels] = pj;
        n_pixels++;
      }
    }

    free(assignment);
    free(cum_dist);
  }

  return n_pixels;
}

// Rasterizes a polyline through ref_i/ref_j reference pixels.
// close_loop != 0 adds a segment from the last ref back to the first.
// use_d4 != 0 uses D4 (cardinal-only) steps; otherwise D8.
// Returns total number of pixels written to out_i/out_j.
TOPOTOOLBOX_API
ptrdiff_t sample_points_between_refs(ptrdiff_t *out_i, ptrdiff_t *out_j,
                                     const ptrdiff_t *ref_i,
                                     const ptrdiff_t *ref_j, ptrdiff_t n_refs,
                                     int close_loop, int use_d4) {
  if (n_refs < 2) return 0;

  ptrdiff_t total = 0;
  ptrdiff_t n_segments = close_loop ? n_refs : n_refs - 1;

  for (ptrdiff_t k = 0; k < n_segments; k++) {
    ptrdiff_t i0 = ref_i[k];
    ptrdiff_t j0 = ref_j[k];
    ptrdiff_t i1 = ref_i[(k + 1) % n_refs];
    ptrdiff_t j1 = ref_j[(k + 1) % n_refs];

    // First segment emits start point; subsequent segments skip it
    // (it was the end point of the previous segment)
    int skip_first = (k > 0) ? 1 : 0;

    ptrdiff_t n;
    if (use_d4) {
      n = bresenham_d4(out_i + total, out_j + total, i0, j0, i1, j1,
                       skip_first);
    } else {
      n = bresenham_d8(out_i + total, out_j + total, i0, j0, i1, j1,
                       skip_first);
    }
    total += n;
  }

  return total;
}

// ============================================================================
// Line simplification — Iterative End-Point Fit (IEF) engine
// ============================================================================
//
// ief_build constructs the full ordered insertion sequence in O(n log n):
//   seq[k]        = original index of the k-th inserted interior point
//   rmse_curve[k] = RMSE of the k+2-point approximation (k=0: endpoints only)
//
// simplify_line selects a stopping point on that curve via one of 7 criteria
// (methods 0–6; see SIMPLIFY_* defines in topotoolbox.h), then reconstructs
// the kept-point set.  Method 5 (VW-area) uses a separate doubly-linked-list
// insertion path and does not call ief_build.

// Perpendicular distance from point p to infinite line through a→b.
static float seg_perp_dist(float ax, float ay, float bx, float by, float px,
                           float py) {
  float dx = bx - ax, dy = by - ay;
  float len2 = dx * dx + dy * dy;
  if (len2 < 1e-12f)
    return sqrtf((px - ax) * (px - ax) + (py - ay) * (py - ay));
  float cross = (px - ax) * dy - (py - ay) * dx;
  return fabsf(cross) / sqrtf(len2);
}

// Area of triangle (a, b, c).
static float seg_tri_area(float ax, float ay, float bx, float by, float cx,
                          float cy) {
  return 0.5f * fabsf((bx - ax) * (cy - ay) - (by - ay) * (cx - ax));
}

// ---- IEF heap (max-deviation via min-heap on neg_dev) ----

typedef struct {
  float neg_dev;  // -max_perp_dist in segment (heap key)
  float seg_rss;  // sum of squared perpendicular distances in segment
  ptrdiff_t lo, hi, max_k;
} ief_seg;

typedef struct {
  ief_seg *data;
  ptrdiff_t size, cap;
} ief_heap;

static void ief_heap_push(ief_heap *h, ief_seg e) {
  if (h->size >= h->cap) {
    h->cap = h->cap ? h->cap * 2 : 64;
    h->data = (ief_seg *)realloc(h->data, (size_t)h->cap * sizeof(ief_seg));
  }
  ptrdiff_t i = h->size++;
  h->data[i] = e;
  while (i > 0) {
    ptrdiff_t p = (i - 1) / 2;
    if (h->data[p].neg_dev <= h->data[i].neg_dev) break;
    ief_seg tmp = h->data[p];
    h->data[p] = h->data[i];
    h->data[i] = tmp;
    i = p;
  }
}

static ief_seg ief_heap_pop(ief_heap *h) {
  ief_seg top = h->data[0];
  h->data[0] = h->data[--h->size];
  ptrdiff_t i = 0;
  for (;;) {
    ptrdiff_t l = 2 * i + 1, r = 2 * i + 2, s = i;
    if (l < h->size && h->data[l].neg_dev < h->data[s].neg_dev) s = l;
    if (r < h->size && h->data[r].neg_dev < h->data[s].neg_dev) s = r;
    if (s == i) break;
    ief_seg tmp = h->data[s];
    h->data[s] = h->data[i];
    h->data[i] = tmp;
    i = s;
  }
  return top;
}

// Scan interior of [lo, hi] for max perp-dist point and segment RSS.
static ief_seg ief_make_seg(const float *in_i, const float *in_j, ptrdiff_t lo,
                            ptrdiff_t hi) {
  ief_seg s;
  s.lo = lo;
  s.hi = hi;
  s.max_k = lo;
  s.seg_rss = 0.0f;
  s.neg_dev = 0.0f;
  float max_d = -1.0f;
  for (ptrdiff_t k = lo + 1; k < hi; k++) {
    float d =
        seg_perp_dist(in_i[lo], in_j[lo], in_i[hi], in_j[hi], in_i[k], in_j[k]);
    s.seg_rss += d * d;
    if (d > max_d) {
      max_d = d;
      s.max_k = k;
    }
  }
  if (hi - lo > 1) s.neg_dev = -max_d;
  return s;
}

// Build IEF insertion sequence.
// seq[k]        = k-th inserted interior point index, k=0..n-3
// rmse_curve[k] = RMSE with k+2 points (k=0..n-2):
//   rmse_curve[0] = endpoints only; rmse_curve[k+1] = after inserting seq[k]
// Returns n-2 (number of interior points).
static ptrdiff_t ief_build(const float *in_i, const float *in_j, ptrdiff_t n,
                           ptrdiff_t *seq, float *rmse_curve) {
  if (n <= 2) return 0;

  ief_heap heap = {NULL, 0, 0};
  ief_seg init = ief_make_seg(in_i, in_j, 0, n - 1);
  double total_rss = (double)init.seg_rss;
  ief_heap_push(&heap, init);

  rmse_curve[0] = (float)sqrt(total_rss / (double)n);

  ptrdiff_t n_steps = n - 2;
  for (ptrdiff_t t = 0; t < n_steps; t++) {
    if (heap.size == 0) {
      seq[t] = 0;
      rmse_curve[t + 1] = 0.0f;
      continue;
    }
    ief_seg top = ief_heap_pop(&heap);
    seq[t] = top.max_k;
    total_rss -= (double)top.seg_rss;

    ief_seg left = ief_make_seg(in_i, in_j, top.lo, top.max_k);
    ief_seg right = ief_make_seg(in_i, in_j, top.max_k, top.hi);
    total_rss += (double)left.seg_rss + (double)right.seg_rss;
    if (total_rss < 0.0) total_rss = 0.0;

    if (left.hi - left.lo > 1) ief_heap_push(&heap, left);
    if (right.hi - right.lo > 1) ief_heap_push(&heap, right);

    rmse_curve[t + 1] = (float)sqrt(total_rss / (double)n);
  }
  free(heap.data);
  return n_steps;
}

// OLS residual sum of squares for y[from..to] vs index x=0..to-from.
static float ols_residual(const float *y, ptrdiff_t from, ptrdiff_t to) {
  ptrdiff_t m = to - from + 1;
  if (m < 2) return 0.0f;
  double sx = 0, sy = 0, sxy = 0, sxx = 0;
  for (ptrdiff_t k = 0; k < m; k++) {
    double x = (double)k;
    double yi = (double)y[from + k];
    sx += x;
    sy += yi;
    sxy += x * yi;
    sxx += x * x;
  }
  double denom = (double)m * sxx - sx * sx;
  if (denom < 1e-15) return 0.0f;
  double slope = ((double)m * sxy - sx * sy) / denom;
  double intercept = (sy - slope * sx) / (double)m;
  float res = 0.0f;
  for (ptrdiff_t k = 0; k < m; k++) {
    double pred = slope * (double)k + intercept;
    double r = (double)y[from + k] - pred;
    res += (float)(r * r);
  }
  return res;
}

TOPOTOOLBOX_API
ptrdiff_t simplify_line(float *out_i, float *out_j, const float *track_i,
                        const float *track_j, ptrdiff_t n_points,
                        float tolerance, int method) {
  if (n_points <= 2) {
    for (ptrdiff_t k = 0; k < n_points; k++) {
      out_i[k] = track_i[k];
      out_j[k] = track_j[k];
    }
    return n_points;
  }

  ptrdiff_t n = n_points;

  // ---- Method 5: VW-area IEF (separate code path) ----
  if (method == 5) {
    ptrdiff_t *prv_kept = (ptrdiff_t *)malloc((size_t)n * sizeof(ptrdiff_t));
    ptrdiff_t *nxt_kept = (ptrdiff_t *)malloc((size_t)n * sizeof(ptrdiff_t));
    float *eff_area = (float *)malloc((size_t)n * sizeof(float));
    uint8_t *inserted = (uint8_t *)calloc((size_t)n, sizeof(uint8_t));

    inserted[0] = 1;
    inserted[n - 1] = 1;
    for (ptrdiff_t k = 1; k < n - 1; k++) {
      prv_kept[k] = 0;
      nxt_kept[k] = n - 1;
      eff_area[k] = seg_tri_area(track_i[0], track_j[0], track_i[k], track_j[k],
                                 track_i[n - 1], track_j[n - 1]);
    }
    min_heap heap;
    heap_init(&heap, n);
    for (ptrdiff_t k = 1; k < n - 1; k++) heap_push(&heap, -eff_area[k], k);

    while (heap.size > 0) {
      heap_entry top = heap_pop(&heap);
      ptrdiff_t k = top.idx;
      if (inserted[k]) continue;
      float area = -top.abs_dist;
      if (area < tolerance) break;
      if (top.abs_dist != -eff_area[k]) continue;  // stale

      inserted[k] = 1;
      ptrdiff_t p = prv_kept[k], nx = nxt_kept[k];

      // Update nxt_kept for uninserted points in (p, k)
      for (ptrdiff_t q = p + 1; q < k; q++) {
        if (!inserted[q]) {
          nxt_kept[q] = k;
          float new_a =
              seg_tri_area(track_i[prv_kept[q]], track_j[prv_kept[q]],
                           track_i[q], track_j[q], track_i[k], track_j[k]);
          eff_area[q] = new_a;
          heap_push(&heap, -new_a, q);
        }
      }
      // Update prv_kept for uninserted points in (k, nx)
      for (ptrdiff_t q = k + 1; q < nx; q++) {
        if (!inserted[q]) {
          prv_kept[q] = k;
          float new_a =
              seg_tri_area(track_i[k], track_j[k], track_i[q], track_j[q],
                           track_i[nxt_kept[q]], track_j[nxt_kept[q]]);
          eff_area[q] = new_a;
          heap_push(&heap, -new_a, q);
        }
      }
    }
    heap_free(&heap);

    ptrdiff_t cnt = 0;
    for (ptrdiff_t k = 0; k < n; k++)
      if (inserted[k]) {
        out_i[cnt] = track_i[k];
        out_j[cnt] = track_j[k];
        cnt++;
      }
    free(prv_kept);
    free(nxt_kept);
    free(eff_area);
    free(inserted);
    return cnt;
  }

  // ---- Methods 0-4, 6: IEF engine + stopping criterion ----

  // rmse_curve[k] = RMSE with k+2 points (k=0..n-2), size n-1
  // seq[k] = k-th inserted interior point (k=0..n-3), size n-2
  ptrdiff_t *seq = (ptrdiff_t *)malloc((size_t)(n - 2) * sizeof(ptrdiff_t));
  float *rmse_curve = (float *)malloc((size_t)(n - 1) * sizeof(float));
  ief_build(track_i, track_j, n, seq, rmse_curve);

  // n_target = number of points to keep (2..n)
  ptrdiff_t n_target = 2;

  if (method == 0) {
    // Fixed N
    n_target = (ptrdiff_t)tolerance;
    if (n_target < 2) n_target = 2;
    if (n_target > n) n_target = n;

  } else if (method == 1) {
    // Kneedle: minimize d[k] = y[k] + x[k] - 1 on normalised RMSE curve
    float y0 = rmse_curve[0];
    ptrdiff_t best_k = 0;
    float best_d = FLT_MAX;
    for (ptrdiff_t k = 0; k < n - 1; k++) {
      float xk = (n > 2) ? (float)k / (float)(n - 2) : 0.0f;
      float yk = (y0 > 1e-12f) ? rmse_curve[k] / y0 : 0.0f;
      float d = yk + xk - 1.0f;
      if (d < best_d) {
        best_d = d;
        best_k = k;
      }
    }
    n_target = best_k + 2;

  } else if (method == 2) {
    // AIC: 2*(k+2) + n*ln(rmse²)
    // tolerance = RMSE noise floor: log term stops improving below it,
    // preventing the criterion from always preferring all n points.
    float min_rmse2 = tolerance * tolerance;
    if (min_rmse2 < 1e-30f) min_rmse2 = 1e-30f;
    ptrdiff_t best_k = 0;
    float best_aic = FLT_MAX;
    for (ptrdiff_t k = 0; k < n - 1; k++) {
      float rmse2 = rmse_curve[k] * rmse_curve[k];
      if (rmse2 < min_rmse2) rmse2 = min_rmse2;
      float aic = 2.0f * (float)(k + 2) + (float)n * logf(rmse2);
      if (aic < best_aic) {
        best_aic = aic;
        best_k = k;
      }
    }
    n_target = best_k + 2;

  } else if (method == 3) {
    // BIC: (k+2)*ln(n) + n*ln(rmse²)
    // tolerance = RMSE noise floor (same reasoning as AIC).
    float min_rmse2 = tolerance * tolerance;
    if (min_rmse2 < 1e-30f) min_rmse2 = 1e-30f;
    float logn = (n > 1) ? logf((float)n) : 0.0f;
    ptrdiff_t best_k = 0;
    float best_bic = FLT_MAX;
    for (ptrdiff_t k = 0; k < n - 1; k++) {
      float rmse2 = rmse_curve[k] * rmse_curve[k];
      if (rmse2 < min_rmse2) rmse2 = min_rmse2;
      float bic = (float)(k + 2) * logn + (float)n * logf(rmse2);
      if (bic < best_bic) {
        best_bic = bic;
        best_k = k;
      }
    }
    n_target = best_k + 2;

  } else if (method == 4) {
    // MDL: 2*(k+2)*log2(coord_range) + n*log2(rmse)
    // tolerance = RMSE noise floor.
    float min_i = track_i[0], max_i = track_i[0];
    float min_j = track_j[0], max_j = track_j[0];
    for (ptrdiff_t k = 1; k < n; k++) {
      if (track_i[k] < min_i) min_i = track_i[k];
      if (track_i[k] > max_i) max_i = track_i[k];
      if (track_j[k] < min_j) min_j = track_j[k];
      if (track_j[k] > max_j) max_j = track_j[k];
    }
    float coord_range =
        (max_i - min_i > max_j - min_j) ? max_i - min_i : max_j - min_j;
    if (coord_range < 1.0f) coord_range = 1.0f;
    float log2_range = log2f(coord_range);
    float log2e = 1.0f / logf(2.0f);
    float min_rmse = (tolerance > 1e-15f) ? tolerance : 1e-15f;

    ptrdiff_t best_k = 0;
    float best_mdl = FLT_MAX;
    for (ptrdiff_t k = 0; k < n - 1; k++) {
      float rmse_k = rmse_curve[k];
      if (rmse_k < min_rmse) rmse_k = min_rmse;
      float mdl =
          2.0f * (float)(k + 2) * log2_range + (float)n * log2e * logf(rmse_k);
      if (mdl < best_mdl) {
        best_mdl = mdl;
        best_k = k;
      }
    }
    n_target = best_k + 2;

  } else {
    // Method 6: L-method — find split m minimising combined OLS residuals
    ptrdiff_t best_m = 0;
    float best_score = FLT_MAX;
    for (ptrdiff_t m = 1; m <= n - 3; m++) {
      float score =
          ols_residual(rmse_curve, 0, m) + ols_residual(rmse_curve, m, n - 2);
      if (score < best_score) {
        best_score = score;
        best_m = m;
      }
    }
    n_target = best_m + 2;
  }

  if (n_target < 2) n_target = 2;
  if (n_target > n) n_target = n;

  // Reconstruct: keep endpoints + seq[0..n_target-3]
  uint8_t *keep = (uint8_t *)calloc((size_t)n, 1);
  keep[0] = 1;
  keep[n - 1] = 1;
  for (ptrdiff_t k = 0; k < n_target - 2; k++) keep[seq[k]] = 1;

  ptrdiff_t cnt = 0;
  for (ptrdiff_t k = 0; k < n; k++)
    if (keep[k]) {
      out_i[cnt] = track_i[k];
      out_j[cnt] = track_j[k];
      cnt++;
    }

  free(seq);
  free(rmse_curve);
  free(keep);
  return cnt;
}
