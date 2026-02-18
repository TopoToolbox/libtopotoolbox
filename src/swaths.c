#define TOPOTOOLBOX_BUILD

#include <float.h>
#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "topotoolbox.h"

/*
  Swath profile analysis functions for libtopotoolbox.

  This module provides two complementary approaches for swath profile analysis:

  1. Transverse: Aggregates DEM pixels by perpendicular distance to track,
     creating a single averaged cross-sectional profile.

  2. Longitudinal: Computes individual swath statistics for each track point,
     showing how the profile changes along the track.

  Core algorithm: frontier-based Dijkstra expansion from the track outward.
  Only pixels within the specified distance are visited, making cost
  proportional to swath area rather than total grid size.
*/

// ============================================================================
// Statistics Accumulator - Welford's online algorithm for numerical stability
// ============================================================================

typedef struct {
  ptrdiff_t count;
  double sum;
  double sum_sq;
  float min_val;
  float max_val;
} swath_stats_accumulator;

// ============================================================================
// Percentile Accumulator - Stores values for percentile computation
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
  if (isnan(value))
    return;
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

static inline float percentile_accumulator_get(const percentile_accumulator *acc,
                                                float p) {
  if (acc->count == 0)
    return NAN;
  if (p < 0.0f) p = 0.0f;
  if (p > 100.0f) p = 100.0f;
  float index = p / 100.0f * (acc->count - 1);
  ptrdiff_t lower = (ptrdiff_t)index;
  ptrdiff_t upper = lower + 1;
  if (upper >= acc->count)
    return acc->values[acc->count - 1];
  float fraction = index - lower;
  return acc->values[lower] * (1.0f - fraction) +
         acc->values[upper] * fraction;
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
  if (isnan(value))
    return;
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
  if (acc->count <= 1)
    return 0.0f;
  double variance = (acc->sum_sq - acc->sum * acc->sum / acc->count) /
                    (acc->count - 1);
  return variance > 0.0 ? (float)sqrt(variance) : 0.0f;
}

// ============================================================================
// Min-Heap for Frontier Expansion (Dijkstra)
// ============================================================================

typedef struct {
  float abs_dist;   // absolute perpendicular distance in pixels (heap key)
  ptrdiff_t idx;    // linear pixel index
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
// Geometric Distance Computation
// ============================================================================

// Signed perpendicular distance from point to line segment (pixel space).
// Negative = left of segment direction, positive = right.
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
  if (t < 0.0f) t = 0.0f;
  else if (t > 1.0f) t = 1.0f;

  *lambda = t;
  *proj_x = ax + t * dx;
  *proj_y = ay + t * dy;

  float to_p_x = px - *proj_x;
  float to_p_y = py - *proj_y;
  float dist = sqrtf(to_p_x * to_p_x + to_p_y * to_p_y);
  float cross = dpx * dy - dpy * dx;

  return cross >= 0.0f ? dist : -dist;
}

// Forward declaration
static void boundary_dijkstra(
    float *restrict dist_out,
    const float *restrict best_abs,
    const ptrdiff_t *seeds, ptrdiff_t n_seeds,
    ptrdiff_t dims[2], float hw_px);

// ============================================================================
// Frontier-based Distance Map (Dijkstra from track, bounded by max_dist_px)
// ============================================================================
//
// Expands outward from track pixels using a min-heap ordered by absolute
// perpendicular distance. Each pixel is processed once. Only pixels within
// max_dist_px (pixel units) are visited.
//
// Outputs (caller-allocated, size total_pixels):
//   best_abs:     absolute perpendicular distance in pixels (FLT_MAX if unvisited)
//   signed_dist:  signed perpendicular distance in pixels (0 if unvisited)
//   nearest_seg:  nearest segment index (-1 if unvisited)
//
// Any output pointer can be NULL to skip that output, but best_abs is always
// needed internally and will be allocated temporarily if NULL.

static void frontier_distance_map(
    float *restrict best_abs,
    float *restrict signed_dist,
    ptrdiff_t *restrict nearest_seg,
    const float *restrict track_i,
    const float *restrict track_j,
    ptrdiff_t n_track_points, ptrdiff_t dims[2],
    float max_dist_px,
    const float *dem,
    const int8_t *mask) {
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
  #define TRY_PIXEL(pi, pj, seg_lo, seg_hi) do {                             \
    if ((pi) >= 0 && (pi) < dims[0] && (pj) >= 0 && (pj) < dims[1]) {       \
      ptrdiff_t _idx = (pj) * dims[0] + (pi);                                \
      if (dem && isnan(dem[_idx])) break;                                      \
      if (mask && !mask[_idx]) break;                                          \
      float _px = (float)(pi);                                                \
      float _py = (float)(pj);                                                \
      float _best_sd = 0.0f, _best_ad = FLT_MAX;                             \
      ptrdiff_t _best_seg = -1;                                               \
      for (ptrdiff_t _sk = (seg_lo); _sk <= (seg_hi); _sk++) {               \
        float _proj_x, _proj_y, _lam;                                        \
        float _d = point_to_segment_distance(_px, _py,                        \
            track_i[_sk], track_j[_sk], track_i[_sk+1], track_j[_sk+1],      \
            &_proj_x, &_proj_y, &_lam);                                      \
        if (fabsf(_d) < _best_ad) {                                           \
          _best_ad = fabsf(_d);                                                \
          _best_sd = _d;                                                       \
          _best_seg = _sk;                                                     \
        }                                                                      \
      }                                                                        \
      if (_best_ad < best_abs[_idx] && _best_ad <= max_dist_px) {             \
        best_abs[_idx] = _best_ad;                                             \
        if (signed_dist) signed_dist[_idx] = _best_sd;                        \
        if (nearest_seg) nearest_seg[_idx] = _best_seg;                       \
        heap_push(&heap, _best_ad, _idx);                                     \
      }                                                                        \
    }                                                                          \
  } while(0)

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

  // 8-connected neighbor offsets
  static const int di8[8] = {-1, -1, -1, 0, 0, 1, 1, 1};
  static const int dj8[8] = {-1, 0, 1, -1, 1, -1, 0, 1};

  // Dijkstra expansion
  while (heap.size > 0) {
    heap_entry top = heap_pop(&heap);
    ptrdiff_t idx = top.idx;

    // Skip stale entries
    if (top.abs_dist > best_abs[idx])
      continue;

    ptrdiff_t ci = idx % dims[0];
    ptrdiff_t cj = idx / dims[0];
    ptrdiff_t seg = nearest_seg ? nearest_seg[idx] : 0;

    // Search range: parent's segment ±2
    ptrdiff_t seg_lo = seg > 1 ? seg - 2 : 0;
    ptrdiff_t seg_hi = seg + 2 < n_track_points - 1 ? seg + 2
                                                      : n_track_points - 2;

    for (int n = 0; n < 8; n++) {
      ptrdiff_t ni = ci + di8[n];
      ptrdiff_t nj = cj + dj8[n];
      TRY_PIXEL(ni, nj, seg_lo, seg_hi);
    }
  }

  #undef TRY_PIXEL

  heap_free(&heap);
  if (free_best) free(best_abs);
}

// ============================================================================
// Public API Functions
// ============================================================================

TOPOTOOLBOX_API
ptrdiff_t swath_compute_nbins(float half_width, float bin_resolution) {
  if (bin_resolution <= 0.0f || half_width <= 0.0f)
    return 0;
  return (ptrdiff_t)(2.0f * half_width / bin_resolution) + 1;
}

// ============================================================================
// Clipped Distance Map with Optional Centre-Line Inversion (Public API)
// ============================================================================

TOPOTOOLBOX_API
ptrdiff_t swath_compute_distance_map(
    float *restrict distance_from_track,
    ptrdiff_t *restrict nearest_segment,
    float *restrict dist_from_boundary,
    float *restrict centre_line_i,
    float *restrict centre_line_j,
    float *restrict centre_width,
    const float *restrict track_i,
    const float *restrict track_j,
    ptrdiff_t n_track_points,
    ptrdiff_t dims[2],
    float cellsize,
    float half_width,
    int compute_signed)
{
  if (n_track_points < 2)
    return 0;

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

  frontier_distance_map(best_abs, signed_dist_px, nseg,
                        track_i, track_j, n_track_points, dims, hw_px,
                        NULL, NULL);

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
  int want_centre = (centre_line_i != NULL && centre_line_j != NULL &&
                     centre_width != NULL);
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
      if (best_abs[idx] == FLT_MAX || best_abs[idx] > hw_px)
        continue;

      ptrdiff_t ci = idx % dims[0];
      ptrdiff_t cj = idx / dims[0];
      int is_boundary = 0;

      for (int n = 0; n < 4; n++) {
        ptrdiff_t ni = ci + dn_i[n];
        ptrdiff_t nj = cj + dn_j[n];
        if (ni < 0 || ni >= dims[0] || nj < 0 || nj >= dims[1]) {
          is_boundary = 1; break;
        }
        ptrdiff_t nidx = nj * dims[0] + ni;
        if (best_abs[nidx] == FLT_MAX || best_abs[nidx] > hw_px) {
          is_boundary = 1; break;
        }
      }
      if (!is_boundary) continue;

      if (signed_dist_px[idx] >= 0.0f) {
        if (n_pos >= pos_cap) {
          pos_cap *= 2;
          pos_seeds = (ptrdiff_t *)realloc(pos_seeds, pos_cap * sizeof(ptrdiff_t));
        }
        pos_seeds[n_pos++] = idx;
      }
      if (signed_dist_px[idx] <= 0.0f) {
        if (n_neg >= neg_cap) {
          neg_cap *= 2;
          neg_seeds = (ptrdiff_t *)realloc(neg_seeds, neg_cap * sizeof(ptrdiff_t));
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

    // Centre line: ridge detection on dfb (pixel-unit min for comparison)
    if (want_centre) {
      for (ptrdiff_t idx = 0; idx < total; idx++) {
        if (best_abs[idx] == FLT_MAX || best_abs[idx] > hw_px)
          continue;
        float dp = dist_pos[idx], dn = dist_neg[idx];
        float val = dp < dn ? dp : dn;  // pixel units for comparison
        if (val == FLT_MAX || val == 0.0f)
          continue;

        ptrdiff_t ci = idx % dims[0];
        ptrdiff_t cj = idx / dims[0];

        int ridge_i = 0;
        if (ci > 0 && ci < dims[0] - 1) {
          float l = dist_pos[idx-1] < dist_neg[idx-1] ? dist_pos[idx-1] : dist_neg[idx-1];
          float r = dist_pos[idx+1] < dist_neg[idx+1] ? dist_pos[idx+1] : dist_neg[idx+1];
          if (val >= l && val >= r) ridge_i = 1;
        }
        int ridge_j = 0;
        if (cj > 0 && cj < dims[1] - 1) {
          ptrdiff_t u = idx - dims[0], d = idx + dims[0];
          float uv = dist_pos[u] < dist_neg[u] ? dist_pos[u] : dist_neg[u];
          float dv = dist_pos[d] < dist_neg[d] ? dist_pos[d] : dist_neg[d];
          if (val >= uv && val >= dv) ridge_j = 1;
        }
        if (!ridge_i && !ridge_j) continue;

        centre_line_i[n_centre] = (float)ci;
        centre_line_j[n_centre] = (float)cj;
        centre_width[n_centre] = (dp + dn) * cellsize;
        n_centre++;
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
// Full (Unclipped) Distance Map (Public API)
// ============================================================================

TOPOTOOLBOX_API
void swath_compute_full_distance_map(
    float *restrict distance,
    ptrdiff_t *restrict nearest_segment,
    const float *restrict track_i,
    const float *restrict track_j,
    ptrdiff_t n_track_points,
    ptrdiff_t dims[2],
    float cellsize,
    const float *dem,
    const int8_t *mask,
    int compute_signed)
{
  if (n_track_points < 2)
    return;

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

  frontier_distance_map(best_abs, signed_dist_px, nseg,
                        track_i, track_j, n_track_points, dims, FLT_MAX,
                        dem, mask);

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
// Transverse Swath Profile - Averaged Cross-Section
// ============================================================================

TOPOTOOLBOX_API
void swath_transverse(
    float *restrict bin_distances, float *restrict bin_means,
    float *restrict bin_stddevs, float *restrict bin_mins,
    float *restrict bin_maxs, ptrdiff_t *restrict bin_counts,
    float *restrict bin_medians, float *restrict bin_q1,
    float *restrict bin_q3, const int *restrict percentile_list,
    ptrdiff_t n_percentiles, float *restrict bin_percentiles,
    const float *restrict dem,
    const float *restrict distance_from_track,
    ptrdiff_t dims[2],
    float cellsize, float half_width, float bin_resolution, ptrdiff_t n_bins,
    int normalize) {
  if (n_bins <= 0)
    return;

  ptrdiff_t total_pixels = dims[0] * dims[1];

  int compute_percentiles =
      (bin_medians != NULL || bin_q1 != NULL || bin_q3 != NULL ||
       (percentile_list != NULL && n_percentiles > 0 && bin_percentiles != NULL));

  // Allocate accumulators
  swath_stats_accumulator *accumulators =
      (swath_stats_accumulator *)calloc(n_bins, sizeof(swath_stats_accumulator));
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
        if (normalize && !isnan(bin_q1[b]))
          bin_q1[b] += reference_elevation;
      }
      if (bin_q3 != NULL) {
        bin_q3[b] = percentile_accumulator_get(&p_accumulators[b], 75.0f);
        if (normalize && !isnan(bin_q3[b]))
          bin_q3[b] += reference_elevation;
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
// Local Tangent via PCA (for orthogonal cross-sections)
// ============================================================================
//
// Computes the local tangent direction at a track point by PCA on a
// window of N neighbouring track points. Returns a unit tangent (ti, tj)
// in pixel space, oriented consistently along the track direction.

static void compute_local_tangent(
    const float *track_i, const float *track_j,
    ptrdiff_t n_track_points,
    ptrdiff_t pt,
    ptrdiff_t n_points_regression,
    float *ti, float *tj)
{
  // Window centred on pt, clamped to track bounds, at least 2 points
  ptrdiff_t half_n = n_points_regression / 2;
  if (half_n < 1) half_n = 1;

  ptrdiff_t win_lo = pt - half_n;
  ptrdiff_t win_hi = pt + half_n;
  if (win_lo < 0) win_lo = 0;
  if (win_hi >= n_track_points) win_hi = n_track_points - 1;
  if (win_hi - win_lo < 1) {
    if (pt < n_track_points - 1) {
      win_lo = pt; win_hi = pt + 1;
    } else {
      win_lo = pt - 1; win_hi = pt;
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
// Longitudinal Swath Profile - Orthogonal Cross-Section Approach
// ============================================================================
//
// Per-point approach: for each track point, compute a local tangent via
// PCA regression, then gather pixels along the orthogonal cross-section
// (binning_distance <= 0) or within the bounding box defined by edge
// orthogonals (binning_distance > 0).

// Forward declaration for Bresenham (defined later in file)
static ptrdiff_t bresenham_d8(ptrdiff_t *out_i, ptrdiff_t *out_j,
                              ptrdiff_t i0, ptrdiff_t j0,
                              ptrdiff_t i1, ptrdiff_t j1,
                              int skip_first);

TOPOTOOLBOX_API
void swath_longitudinal(
    float *restrict point_means, float *restrict point_stddevs,
    float *restrict point_mins, float *restrict point_maxs,
    ptrdiff_t *restrict point_counts, float *restrict point_medians,
    float *restrict point_q1, float *restrict point_q3,
    const int *restrict percentile_list, ptrdiff_t n_percentiles,
    float *restrict point_percentiles, const float *restrict dem,
    const float *restrict track_i, const float *restrict track_j,
    ptrdiff_t n_track_points,
    const float *restrict distance_from_track,
    ptrdiff_t dims[2], float cellsize,
    float half_width, float binning_distance,
    ptrdiff_t n_points_regression,
    ptrdiff_t n_internal_cuts) {
  if (n_track_points < 2)
    return;

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

  // Initialize per-point accumulators
  swath_stats_accumulator *accumulators = (swath_stats_accumulator *)calloc(
      n_track_points, sizeof(swath_stats_accumulator));
  for (ptrdiff_t k = 0; k < n_track_points; k++)
    accumulator_init(&accumulators[k]);

  percentile_accumulator *p_accumulators = NULL;
  if (compute_percentiles) {
    p_accumulators = (percentile_accumulator *)calloc(
        n_track_points, sizeof(percentile_accumulator));
    for (ptrdiff_t k = 0; k < n_track_points; k++)
      percentile_accumulator_init(&p_accumulators[k], 64);
  }

  // Bresenham work buffers (reused across track points)
  // Max Bresenham length for hw_px on each side: 2*hw_px + 1 is sufficient
  ptrdiff_t bres_cap = (ptrdiff_t)(2.0f * hw_px) + 16;
  ptrdiff_t *bres_i = (ptrdiff_t *)malloc(bres_cap * sizeof(ptrdiff_t));
  ptrdiff_t *bres_j = (ptrdiff_t *)malloc(bres_cap * sizeof(ptrdiff_t));

  // Per-point processing
  for (ptrdiff_t pt = 0; pt < n_track_points; pt++) {
    float tang_i, tang_j;
    compute_local_tangent(track_i, track_j, n_track_points, pt,
                          n_points_regression, &tang_i, &tang_j);
    // Orthogonal direction
    float orth_i = -tang_j;
    float orth_j = tang_i;

    if (binning_distance <= 0.0f) {
      // --- Case 1: single cross-section via Bresenham ---
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
    } else {
      // --- Case 2: regression rectangle on binning window ---
      //
      // Linear regression (PCA) on all track points in [pt_lo..pt_hi]
      // gives one tangent direction. Build one oriented rectangle along
      // that direction. Pixel is selected if inside the rectangle AND
      // within the distance map band.

      // Find edge track points of the binning window
      ptrdiff_t pt_lo = pt, pt_hi = pt;
      while (pt_lo > 0 &&
             cum_dist[pt] - cum_dist[pt_lo - 1] <= binning_distance)
        pt_lo--;
      while (pt_hi < n_track_points - 1 &&
             cum_dist[pt_hi + 1] - cum_dist[pt] <= binning_distance)
        pt_hi++;

      // PCA on track points [pt_lo..pt_hi]
      ptrdiff_t nw = pt_hi - pt_lo + 1;
      float mean_i = 0, mean_j = 0;
      for (ptrdiff_t k = pt_lo; k <= pt_hi; k++) {
        mean_i += track_i[k];
        mean_j += track_j[k];
      }
      mean_i /= (float)nw;
      mean_j /= (float)nw;

      float cov_ii = 0, cov_ij = 0, cov_jj = 0;
      for (ptrdiff_t k = pt_lo; k <= pt_hi; k++) {
        float di = track_i[k] - mean_i;
        float dj = track_j[k] - mean_j;
        cov_ii += di * di;
        cov_ij += di * dj;
        cov_jj += dj * dj;
      }

      // Principal eigenvector (tangent direction)
      float diff = cov_ii - cov_jj;
      float D = sqrtf(diff * diff + 4.0f * cov_ij * cov_ij);
      float reg_ti, reg_tj;
      if (cov_ii >= cov_jj) {
        reg_ti = diff + D;
        reg_tj = 2.0f * cov_ij;
      } else {
        reg_ti = 2.0f * cov_ij;
        reg_tj = -diff + D;
      }
      float vlen = sqrtf(reg_ti * reg_ti + reg_tj * reg_tj);
      if (vlen > 1e-10f) {
        reg_ti /= vlen;
        reg_tj /= vlen;
      } else {
        // Degenerate: fall back to endpoint direction
        reg_ti = track_i[pt_hi] - track_i[pt_lo];
        reg_tj = track_j[pt_hi] - track_j[pt_lo];
        vlen = sqrtf(reg_ti * reg_ti + reg_tj * reg_tj);
        if (vlen > 0) { reg_ti /= vlen; reg_tj /= vlen; }
        else { reg_ti = 1; reg_tj = 0; }
      }

      // Orient consistently with track direction
      float fwd_i = track_i[pt_hi] - track_i[pt_lo];
      float fwd_j = track_j[pt_hi] - track_j[pt_lo];
      if (fwd_i * reg_ti + fwd_j * reg_tj < 0) {
        reg_ti = -reg_ti;
        reg_tj = -reg_tj;
      }

      float reg_oi = -reg_tj, reg_oj = reg_ti;  // orthogonal

      // Project all window track points onto tangent to get along-track extent
      float along_min = 0, along_max = 0;
      for (ptrdiff_t k = pt_lo; k <= pt_hi; k++) {
        float d = (track_i[k] - mean_i) * reg_ti +
                  (track_j[k] - mean_j) * reg_tj;
        if (d < along_min) along_min = d;
        if (d > along_max) along_max = d;
      }
      // Small buffer for edge pixels
      along_min -= 1.0f;
      along_max += 1.0f;

      // Rectangle corners (centered on mean, along tangent, ±hw_px orthogonal)
      float c0_i = mean_i + reg_ti * along_min + reg_oi * hw_px;
      float c0_j = mean_j + reg_tj * along_min + reg_oj * hw_px;
      float c1_i = mean_i + reg_ti * along_max + reg_oi * hw_px;
      float c1_j = mean_j + reg_tj * along_max + reg_oj * hw_px;
      float c2_i = mean_i + reg_ti * along_max - reg_oi * hw_px;
      float c2_j = mean_j + reg_tj * along_max - reg_oj * hw_px;
      float c3_i = mean_i + reg_ti * along_min - reg_oi * hw_px;
      float c3_j = mean_j + reg_tj * along_min - reg_oj * hw_px;

      // Axis-aligned bounding box of the 4 corners
      float corners_i[4] = {c0_i, c1_i, c2_i, c3_i};
      float corners_j[4] = {c0_j, c1_j, c2_j, c3_j};
      float bb_i0 = corners_i[0], bb_i1 = corners_i[0];
      float bb_j0 = corners_j[0], bb_j1 = corners_j[0];
      for (int c = 1; c < 4; c++) {
        if (corners_i[c] < bb_i0) bb_i0 = corners_i[c];
        if (corners_i[c] > bb_i1) bb_i1 = corners_i[c];
        if (corners_j[c] < bb_j0) bb_j0 = corners_j[c];
        if (corners_j[c] > bb_j1) bb_j1 = corners_j[c];
      }

      ptrdiff_t i0 = (ptrdiff_t)bb_i0 - 1;
      ptrdiff_t i1 = (ptrdiff_t)bb_i1 + 2;
      ptrdiff_t j0 = (ptrdiff_t)bb_j0 - 1;
      ptrdiff_t j1 = (ptrdiff_t)bb_j1 + 2;
      if (i0 < 0) i0 = 0;
      if (j0 < 0) j0 = 0;
      if (i1 >= dims[0]) i1 = dims[0] - 1;
      if (j1 >= dims[1]) j1 = dims[1] - 1;

      for (ptrdiff_t pj = j0; pj <= j1; pj++) {
        for (ptrdiff_t pi = i0; pi <= i1; pi++) {
          ptrdiff_t idx = pj * dims[0] + pi;
          float dist_m = distance_from_track[idx];
          if (isnan(dist_m)) continue;
          if (fabsf(dist_m) > half_width) continue;
          if (isnan(dem[idx])) continue;

          // Oriented rectangle check: project onto tangent and orthogonal
          float di = (float)pi - mean_i;
          float dj = (float)pj - mean_j;
          float along = di * reg_ti + dj * reg_tj;
          float across = di * reg_oi + dj * reg_oj;
          if (along < along_min || along > along_max) continue;
          if (across < -hw_px || across > hw_px) continue;

          accumulator_add(&accumulators[pt], dem[idx]);
          if (compute_percentiles)
            percentile_accumulator_add(&p_accumulators[pt], dem[idx]);
        }
      }
    }
  }

  // Finalize statistics
  for (ptrdiff_t k = 0; k < n_track_points; k++) {
    point_means[k] = accumulator_mean(&accumulators[k]);
    point_stddevs[k] = accumulator_stddev(&accumulators[k]);
    point_mins[k] = accumulators[k].count > 0 ? accumulators[k].min_val : NAN;
    point_maxs[k] = accumulators[k].count > 0 ? accumulators[k].max_val : NAN;
    point_counts[k] = accumulators[k].count;

    if (compute_percentiles) {
      percentile_accumulator_sort(&p_accumulators[k]);

      if (point_medians != NULL)
        point_medians[k] =
            percentile_accumulator_get(&p_accumulators[k], 50.0f);
      if (point_q1 != NULL)
        point_q1[k] = percentile_accumulator_get(&p_accumulators[k], 25.0f);
      if (point_q3 != NULL)
        point_q3[k] = percentile_accumulator_get(&p_accumulators[k], 75.0f);

      if (percentile_list != NULL && n_percentiles > 0 &&
          point_percentiles != NULL) {
        for (ptrdiff_t p = 0; p < n_percentiles; p++) {
          point_percentiles[k * n_percentiles + p] =
              percentile_accumulator_get(&p_accumulators[k],
                                         (float)percentile_list[p]);
        }
      }
    }
  }

  // Cleanup
  free(bres_i);
  free(bres_j);
  if (cum_dist) free(cum_dist);
  free(accumulators);
  if (compute_percentiles) {
    for (ptrdiff_t k = 0; k < n_track_points; k++)
      percentile_accumulator_free(&p_accumulators[k]);
    free(p_accumulators);
  }
}

// ============================================================================
// Per-Point Pixel Retrieval (Orthogonal Cross-Section)
// ============================================================================
//
// Returns pixel coordinates associated with a single track point using
// the same regression + orthogonal approach as swath_longitudinal.

TOPOTOOLBOX_API
ptrdiff_t swath_get_point_pixels(
    ptrdiff_t *restrict pixels_i, ptrdiff_t *restrict pixels_j,
    const float *restrict track_i, const float *restrict track_j,
    ptrdiff_t n_track_points, ptrdiff_t point_index,
    const float *restrict distance_from_track,
    ptrdiff_t dims[2],
    float cellsize, float half_width, float binning_distance,
    ptrdiff_t n_points_regression,
    ptrdiff_t n_internal_cuts) {
  if (n_track_points < 2)
    return 0;
  if (point_index < 0 || point_index >= n_track_points)
    return 0;

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
    // --- Case 2: regression rectangle on binning window ---

    // Cumulative along-track distance
    float *cum_dist = (float *)malloc(n_track_points * sizeof(float));
    cum_dist[0] = 0.0f;
    for (ptrdiff_t k = 0; k < n_track_points - 1; k++) {
      float di = track_i[k + 1] - track_i[k];
      float dj = track_j[k + 1] - track_j[k];
      cum_dist[k + 1] = cum_dist[k] + sqrtf(di * di + dj * dj) * cellsize;
    }

    // Find edge track points of the binning window
    ptrdiff_t pt_lo = point_index, pt_hi = point_index;
    while (pt_lo > 0 &&
           cum_dist[point_index] - cum_dist[pt_lo - 1] <= binning_distance)
      pt_lo--;
    while (pt_hi < n_track_points - 1 &&
           cum_dist[pt_hi + 1] - cum_dist[point_index] <= binning_distance)
      pt_hi++;

    // PCA on track points [pt_lo..pt_hi]
    ptrdiff_t nw = pt_hi - pt_lo + 1;
    float mean_i = 0, mean_j = 0;
    for (ptrdiff_t k = pt_lo; k <= pt_hi; k++) {
      mean_i += track_i[k];
      mean_j += track_j[k];
    }
    mean_i /= (float)nw;
    mean_j /= (float)nw;

    float cov_ii = 0, cov_ij = 0, cov_jj = 0;
    for (ptrdiff_t k = pt_lo; k <= pt_hi; k++) {
      float di = track_i[k] - mean_i;
      float dj = track_j[k] - mean_j;
      cov_ii += di * di;
      cov_ij += di * dj;
      cov_jj += dj * dj;
    }

    float diff = cov_ii - cov_jj;
    float D = sqrtf(diff * diff + 4.0f * cov_ij * cov_ij);
    float reg_ti, reg_tj;
    if (cov_ii >= cov_jj) {
      reg_ti = diff + D;
      reg_tj = 2.0f * cov_ij;
    } else {
      reg_ti = 2.0f * cov_ij;
      reg_tj = -diff + D;
    }
    float vlen = sqrtf(reg_ti * reg_ti + reg_tj * reg_tj);
    if (vlen > 1e-10f) {
      reg_ti /= vlen;
      reg_tj /= vlen;
    } else {
      reg_ti = track_i[pt_hi] - track_i[pt_lo];
      reg_tj = track_j[pt_hi] - track_j[pt_lo];
      vlen = sqrtf(reg_ti * reg_ti + reg_tj * reg_tj);
      if (vlen > 0) { reg_ti /= vlen; reg_tj /= vlen; }
      else { reg_ti = 1; reg_tj = 0; }
    }

    float fwd_i = track_i[pt_hi] - track_i[pt_lo];
    float fwd_j = track_j[pt_hi] - track_j[pt_lo];
    if (fwd_i * reg_ti + fwd_j * reg_tj < 0) {
      reg_ti = -reg_ti;
      reg_tj = -reg_tj;
    }

    float reg_oi = -reg_tj, reg_oj = reg_ti;

    // Project all window track points onto tangent to get along-track extent
    float along_min = 0, along_max = 0;
    for (ptrdiff_t k = pt_lo; k <= pt_hi; k++) {
      float d = (track_i[k] - mean_i) * reg_ti +
                (track_j[k] - mean_j) * reg_tj;
      if (d < along_min) along_min = d;
      if (d > along_max) along_max = d;
    }
    along_min -= 1.0f;
    along_max += 1.0f;

    // Rectangle corners
    float c0_i = mean_i + reg_ti * along_min + reg_oi * hw_px;
    float c0_j = mean_j + reg_tj * along_min + reg_oj * hw_px;
    float c1_i = mean_i + reg_ti * along_max + reg_oi * hw_px;
    float c1_j = mean_j + reg_tj * along_max + reg_oj * hw_px;
    float c2_i = mean_i + reg_ti * along_max - reg_oi * hw_px;
    float c2_j = mean_j + reg_tj * along_max - reg_oj * hw_px;
    float c3_i = mean_i + reg_ti * along_min - reg_oi * hw_px;
    float c3_j = mean_j + reg_tj * along_min - reg_oj * hw_px;

    float corners_i[4] = {c0_i, c1_i, c2_i, c3_i};
    float corners_j[4] = {c0_j, c1_j, c2_j, c3_j};
    float bb_i0 = corners_i[0], bb_i1 = corners_i[0];
    float bb_j0 = corners_j[0], bb_j1 = corners_j[0];
    for (int c = 1; c < 4; c++) {
      if (corners_i[c] < bb_i0) bb_i0 = corners_i[c];
      if (corners_i[c] > bb_i1) bb_i1 = corners_i[c];
      if (corners_j[c] < bb_j0) bb_j0 = corners_j[c];
      if (corners_j[c] > bb_j1) bb_j1 = corners_j[c];
    }

    ptrdiff_t i0 = (ptrdiff_t)bb_i0 - 1;
    ptrdiff_t i1 = (ptrdiff_t)bb_i1 + 2;
    ptrdiff_t j0 = (ptrdiff_t)bb_j0 - 1;
    ptrdiff_t j1 = (ptrdiff_t)bb_j1 + 2;
    if (i0 < 0) i0 = 0;
    if (j0 < 0) j0 = 0;
    if (i1 >= dims[0]) i1 = dims[0] - 1;
    if (j1 >= dims[1]) j1 = dims[1] - 1;

    for (ptrdiff_t pj = j0; pj <= j1; pj++) {
      for (ptrdiff_t pi = i0; pi <= i1; pi++) {
        ptrdiff_t idx = pj * dims[0] + pi;
        float dist_m = distance_from_track[idx];
        if (isnan(dist_m)) continue;
        if (fabsf(dist_m) > half_width) continue;

        float di = (float)pi - mean_i;
        float dj = (float)pj - mean_j;
        float along = di * reg_ti + dj * reg_tj;
        float across = di * reg_oi + dj * reg_oj;
        if (along < along_min || along > along_max) continue;
        if (across < -hw_px || across > hw_px) continue;

        pixels_i[n_pixels] = pi;
        pixels_j[n_pixels] = pj;
        n_pixels++;
      }
    }

    free(cum_dist);
  }

  return n_pixels;
}

// ============================================================================
// Inverted Distance Map — Centre Line + Width from Boundary-Inward Dijkstra
// ============================================================================
//
// Given a track and half_width, computes the swath mask via the global distance
// map, then grows inward from the mask boundary using two Dijkstra passes
// (one from each side's boundary). The centre line is the ridge of the
// distance-from-boundary field. Width at each centre-line pixel is the sum
// of distances from the positive-side and negative-side boundaries.

// Helper: single-source Dijkstra from a set of boundary seed pixels, expanding
// D8 inward within the mask. Fills dist_out[total] with Euclidean distance
// from nearest boundary seed.
static void boundary_dijkstra(
    float *restrict dist_out,
    const float *restrict best_abs,
    const ptrdiff_t *seeds, ptrdiff_t n_seeds,
    ptrdiff_t dims[2], float hw_px)
{
  ptrdiff_t total = dims[0] * dims[1];

  for (ptrdiff_t i = 0; i < total; i++)
    dist_out[i] = FLT_MAX;

  min_heap heap;
  heap_init(&heap, n_seeds * 4);

  // Seed boundary pixels with distance 0
  for (ptrdiff_t s = 0; s < n_seeds; s++) {
    ptrdiff_t idx = seeds[s];
    dist_out[idx] = 0.0f;
    heap_push(&heap, 0.0f, idx);
  }

  // D8 neighbor offsets and distances
  static const int di8[8] = {-1, -1, -1, 0, 0, 1, 1, 1};
  static const int dj8[8] = {-1, 0, 1, -1, 1, -1, 0, 1};
  static const float dd8[8] = {1.41421356f, 1.0f, 1.41421356f, 1.0f,
                                1.0f, 1.41421356f, 1.0f, 1.41421356f};

  while (heap.size > 0) {
    heap_entry top = heap_pop(&heap);
    ptrdiff_t idx = top.idx;

    if (top.abs_dist > dist_out[idx])
      continue;  // stale

    ptrdiff_t ci = idx % dims[0];
    ptrdiff_t cj = idx / dims[0];

    for (int n = 0; n < 8; n++) {
      ptrdiff_t ni = ci + di8[n];
      ptrdiff_t nj = cj + dj8[n];
      if (ni < 0 || ni >= dims[0] || nj < 0 || nj >= dims[1]) continue;

      ptrdiff_t nidx = nj * dims[0] + ni;
      // Only expand within the swath mask
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

// swath_invert_distance_map removed — logic now in swath_compute_distance_map

// ============================================================================
// Bresenham Path Rasterization Between Reference Points
// ============================================================================

// D8 Bresenham: standard line algorithm, each step moves 1 pixel
// (cardinal or diagonal).
static ptrdiff_t bresenham_d8(ptrdiff_t *out_i, ptrdiff_t *out_j,
                              ptrdiff_t i0, ptrdiff_t j0,
                              ptrdiff_t i1, ptrdiff_t j1,
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
static ptrdiff_t bresenham_d4(ptrdiff_t *out_i, ptrdiff_t *out_j,
                              ptrdiff_t i0, ptrdiff_t j0,
                              ptrdiff_t i1, ptrdiff_t j1,
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

TOPOTOOLBOX_API
ptrdiff_t sample_points_between_refs(
    ptrdiff_t *out_i,
    ptrdiff_t *out_j,
    const ptrdiff_t *ref_i,
    const ptrdiff_t *ref_j,
    ptrdiff_t n_refs,
    int close_loop,
    int use_d4)
{
  if (n_refs < 2)
    return 0;

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
