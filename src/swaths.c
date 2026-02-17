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

// Find nearest segment in a range. Returns signed distance in pixel space.
static ptrdiff_t find_nearest_segment(float px, float py,
                                      const float *restrict track_i,
                                      const float *restrict track_j,
                                      ptrdiff_t n_track_points,
                                      ptrdiff_t search_start,
                                      ptrdiff_t search_end, float *min_dist,
                                      float *proj_x, float *proj_y,
                                      float *best_lambda) {
  ptrdiff_t best_segment = -1;
  *min_dist = FLT_MAX;
  *best_lambda = 0.0f;

  if (search_start < 0) search_start = 0;
  if (search_end >= n_track_points) search_end = n_track_points - 1;

  for (ptrdiff_t k = search_start; k < search_end; k++) {
    float proj_i, proj_j, lambda;
    float dist = point_to_segment_distance(px, py,
        track_i[k], track_j[k], track_i[k + 1], track_j[k + 1],
        &proj_i, &proj_j, &lambda);

    if (fabsf(dist) < fabsf(*min_dist)) {
      *min_dist = dist;
      best_segment = k;
      *proj_x = proj_i;
      *proj_y = proj_j;
      *best_lambda = lambda;
    }
  }

  return best_segment;
}

// ============================================================================
// Tree March on Global Distance Map
// ============================================================================
//
// Starting from seed pixels (rasterized sub-track), marches D4 outward on the
// GLOBAL distance map following increasing abs_dist. Each side of the track
// (positive/negative signed_dist) is processed sequentially. The result is a
// tree of activated pixels (each visited once). Holes enclosed by the tree
// within its bounding box are filled via flood-fill from the boundary.
//
// pixels_ptr/pixels_cap: dynamic buffer for output linear pixel indices.
// best_abs, signed_dist: global distance map arrays [total_pixels].
// seed_indices/n_seeds: rasterized sub-track pixel linear indices.
// visited: char[total_pixels], must be 0 on entry, restored to 0 on exit.
//
// Returns count of activated pixels (seeds + march + fill).

static ptrdiff_t tree_march_fill(
    ptrdiff_t **pixels_ptr, ptrdiff_t *pixels_cap,
    const float *best_abs,
    const float *signed_dist,
    const ptrdiff_t *seed_indices, ptrdiff_t n_seeds,
    ptrdiff_t dims[2], float hw_px,
    char *visited)
{
  // D4 neighbor offsets
  static const int dn_i[4] = {-1, 1, 0, 0};
  static const int dn_j[4] = {0, 0, -1, 1};

  ptrdiff_t count = 0;
  ptrdiff_t *pixels = *pixels_ptr;
  ptrdiff_t cap = *pixels_cap;

  #define TM_ADD(idx_val) do { \
    if (count >= cap) { cap *= 2; \
      pixels = (ptrdiff_t *)realloc(pixels, cap * sizeof(ptrdiff_t)); } \
    pixels[count++] = (idx_val); \
  } while(0)

  // --- Add seeds (track pixels present in the global distance map) ---
  for (ptrdiff_t s = 0; s < n_seeds; s++) {
    ptrdiff_t idx = seed_indices[s];
    if (visited[idx]) continue;
    if (best_abs[idx] == FLT_MAX) continue;
    visited[idx] = 1;
    TM_ADD(idx);
  }
  ptrdiff_t n_seeds_actual = count;
  if (n_seeds_actual == 0) goto cleanup;

  // --- Positive side march (signed_dist >= 0): BFS uphill from seeds ---
  {
    ptrdiff_t q_head = 0;
    while (q_head < count) {
      ptrdiff_t idx = pixels[q_head++];
      ptrdiff_t ci = idx % dims[0];
      ptrdiff_t cj = idx / dims[0];
      float my_abs = best_abs[idx];

      for (int n = 0; n < 4; n++) {
        ptrdiff_t ni = ci + dn_i[n];
        ptrdiff_t nj = cj + dn_j[n];
        if (ni < 0 || ni >= dims[0] || nj < 0 || nj >= dims[1]) continue;

        ptrdiff_t nidx = nj * dims[0] + ni;
        if (visited[nidx]) continue;
        if (best_abs[nidx] == FLT_MAX) continue;
        if (best_abs[nidx] > hw_px) continue;
        if (best_abs[nidx] <= my_abs) continue;   // uphill only
        if (signed_dist[nidx] < 0.0f) continue;   // positive side

        visited[nidx] = 1;
        TM_ADD(nidx);
      }
    }
  }

  // --- Negative side march (signed_dist <= 0): re-expand from seeds ---
  {
    // Seed the negative side from original seed pixels
    ptrdiff_t neg_start = count;
    for (ptrdiff_t s = 0; s < n_seeds_actual; s++) {
      ptrdiff_t idx = pixels[s];
      ptrdiff_t ci = idx % dims[0];
      ptrdiff_t cj = idx / dims[0];
      float my_abs = best_abs[idx];

      for (int n = 0; n < 4; n++) {
        ptrdiff_t ni = ci + dn_i[n];
        ptrdiff_t nj = cj + dn_j[n];
        if (ni < 0 || ni >= dims[0] || nj < 0 || nj >= dims[1]) continue;

        ptrdiff_t nidx = nj * dims[0] + ni;
        if (visited[nidx]) continue;
        if (best_abs[nidx] == FLT_MAX) continue;
        if (best_abs[nidx] > hw_px) continue;
        if (best_abs[nidx] <= my_abs) continue;
        if (signed_dist[nidx] > 0.0f) continue;   // negative side

        visited[nidx] = 1;
        TM_ADD(nidx);
      }
    }

    ptrdiff_t q_head = neg_start;
    while (q_head < count) {
      ptrdiff_t idx = pixels[q_head++];
      ptrdiff_t ci = idx % dims[0];
      ptrdiff_t cj = idx / dims[0];
      float my_abs = best_abs[idx];

      for (int n = 0; n < 4; n++) {
        ptrdiff_t ni = ci + dn_i[n];
        ptrdiff_t nj = cj + dn_j[n];
        if (ni < 0 || ni >= dims[0] || nj < 0 || nj >= dims[1]) continue;

        ptrdiff_t nidx = nj * dims[0] + ni;
        if (visited[nidx]) continue;
        if (best_abs[nidx] == FLT_MAX) continue;
        if (best_abs[nidx] > hw_px) continue;
        if (best_abs[nidx] <= my_abs) continue;
        if (signed_dist[nidx] > 0.0f) continue;   // negative side

        visited[nidx] = 1;
        TM_ADD(nidx);
      }
    }
  }

  // --- Fill holes within bounding box of activated pixels ---
  // Flood-fill from the bounding-box boundary to mark exterior pixels,
  // then fill anything interior that is within hw_px.
  if (count > 0) {
    ptrdiff_t bb_i0 = pixels[0] % dims[0], bb_i1 = bb_i0;
    ptrdiff_t bb_j0 = pixels[0] / dims[0], bb_j1 = bb_j0;
    for (ptrdiff_t p = 1; p < count; p++) {
      ptrdiff_t pi = pixels[p] % dims[0];
      ptrdiff_t pj = pixels[p] / dims[0];
      if (pi < bb_i0) bb_i0 = pi;
      if (pi > bb_i1) bb_i1 = pi;
      if (pj < bb_j0) bb_j0 = pj;
      if (pj > bb_j1) bb_j1 = pj;
    }

    ptrdiff_t bb_w = bb_i1 - bb_i0 + 1;
    ptrdiff_t bb_h = bb_j1 - bb_j0 + 1;
    ptrdiff_t bb_total = bb_w * bb_h;

    char *exterior = (char *)calloc(bb_total, sizeof(char));
    ptrdiff_t *fq = (ptrdiff_t *)malloc(bb_total * sizeof(ptrdiff_t));
    ptrdiff_t fq_head = 0, fq_tail = 0;

    // Seed exterior flood-fill from bounding-box border pixels that
    // are NOT activated
    for (ptrdiff_t bj = 0; bj < bb_h; bj++) {
      for (ptrdiff_t bi = 0; bi < bb_w; bi++) {
        if (bj > 0 && bj < bb_h - 1 && bi > 0 && bi < bb_w - 1) continue;
        ptrdiff_t gidx = (bb_j0 + bj) * dims[0] + (bb_i0 + bi);
        ptrdiff_t bidx = bj * bb_w + bi;
        if (!visited[gidx]) {
          exterior[bidx] = 1;
          fq[fq_tail++] = bidx;
        }
      }
    }

    // Flood-fill exterior
    while (fq_head < fq_tail) {
      ptrdiff_t bidx = fq[fq_head++];
      ptrdiff_t bi = bidx % bb_w;
      ptrdiff_t bj = bidx / bb_w;

      for (int n = 0; n < 4; n++) {
        ptrdiff_t nbi = bi + dn_i[n];
        ptrdiff_t nbj = bj + dn_j[n];
        if (nbi < 0 || nbi >= bb_w || nbj < 0 || nbj >= bb_h) continue;

        ptrdiff_t nbidx = nbj * bb_w + nbi;
        if (exterior[nbidx]) continue;

        ptrdiff_t gidx = (bb_j0 + nbj) * dims[0] + (bb_i0 + nbi);
        if (!visited[gidx]) {
          exterior[nbidx] = 1;
          fq[fq_tail++] = nbidx;
        }
      }
    }

    // Fill interior holes: non-exterior, non-activated, within hw_px
    for (ptrdiff_t bj = 0; bj < bb_h; bj++) {
      for (ptrdiff_t bi = 0; bi < bb_w; bi++) {
        ptrdiff_t bidx = bj * bb_w + bi;
        if (exterior[bidx]) continue;

        ptrdiff_t gidx = (bb_j0 + bj) * dims[0] + (bb_i0 + bi);
        if (visited[gidx]) continue;
        if (best_abs[gidx] == FLT_MAX) continue;
        if (best_abs[gidx] > hw_px) continue;

        visited[gidx] = 1;
        TM_ADD(gidx);
      }
    }

    free(exterior);
    free(fq);
  }

cleanup:
  // Restore visited for all touched pixels
  for (ptrdiff_t p = 0; p < count; p++) {
    visited[pixels[p]] = 0;
  }

  #undef TM_ADD

  *pixels_ptr = pixels;
  *pixels_cap = cap;
  return count;
}

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
    float max_dist_px) {
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
// Distance Map (Public API)
// ============================================================================

TOPOTOOLBOX_API
void swath_distance_map(float *restrict distance,
                        ptrdiff_t *restrict nearest_segment,
                        const float *restrict track_i,
                        const float *restrict track_j,
                        ptrdiff_t n_track_points, ptrdiff_t dims[2],
                        float cellsize) {
  if (n_track_points < 2)
    return;

  ptrdiff_t total_pixels = dims[0] * dims[1];

  // Allocate temporary nearest_segment if caller doesn't want it
  ptrdiff_t *temp_nearest = NULL;
  int need_free = 0;
  if (nearest_segment == NULL) {
    temp_nearest = (ptrdiff_t *)malloc(total_pixels * sizeof(ptrdiff_t));
    need_free = 1;
  } else {
    temp_nearest = nearest_segment;
  }

  // Use frontier with no distance cutoff (process all pixels)
  float *best_abs = (float *)malloc(total_pixels * sizeof(float));
  frontier_distance_map(best_abs, distance, temp_nearest,
                        track_i, track_j, n_track_points, dims, FLT_MAX);

  // Convert signed distance from pixels to meters
  for (ptrdiff_t i = 0; i < total_pixels; i++) {
    distance[i] *= cellsize;
  }

  free(best_abs);
  if (need_free) free(temp_nearest);
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
    const float *restrict dem, const float *restrict track_i,
    const float *restrict track_j, ptrdiff_t n_track_points, ptrdiff_t dims[2],
    float cellsize, float half_width, float bin_resolution, ptrdiff_t n_bins,
    int normalize) {
  if (n_track_points < 2 || n_bins <= 0)
    return;

  ptrdiff_t total_pixels = dims[0] * dims[1];
  float hw_px = half_width / cellsize;

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

  // Frontier expansion bounded by half_width
  float *best_abs = (float *)malloc(total_pixels * sizeof(float));
  float *signed_dist_px = (float *)malloc(total_pixels * sizeof(float));
  frontier_distance_map(best_abs, signed_dist_px, NULL,
                        track_i, track_j, n_track_points, dims, hw_px);

  // Normalization: mean elevation near track center
  float reference_elevation = 0.0f;
  if (normalize) {
    swath_stats_accumulator track_acc;
    accumulator_init(&track_acc);
    float bin_res_px = bin_resolution / cellsize;
    for (ptrdiff_t idx = 0; idx < total_pixels; idx++) {
      if (best_abs[idx] <= bin_res_px) {
        accumulator_add(&track_acc, dem[idx]);
      }
    }
    reference_elevation = accumulator_mean(&track_acc);
  }

  // Bin pixels by signed perpendicular distance
  for (ptrdiff_t idx = 0; idx < total_pixels; idx++) {
    if (best_abs[idx] > hw_px) continue;  // unvisited or outside
    if (isnan(dem[idx])) continue;

    float dist_m = signed_dist_px[idx] * cellsize;

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

  free(best_abs);
  free(signed_dist_px);

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
// Longitudinal Swath Profile - Along-Track Variation
// ============================================================================
//
// Uses frontier-based distance map (bounded by half_width) to find all
// relevant pixels in one pass. Each pixel's along-track position is computed
// from its nearest segment + lambda, then binary search finds matching
// track points within ±binning_distance.

// Binary search helpers on sorted cum_dist array
static ptrdiff_t lower_bound(const float *cum_dist, ptrdiff_t n, float value) {
  ptrdiff_t lo = 0, hi = n;
  while (lo < hi) {
    ptrdiff_t mid = lo + (hi - lo) / 2;
    if (cum_dist[mid] < value) lo = mid + 1;
    else hi = mid;
  }
  return lo;
}

static ptrdiff_t upper_bound(const float *cum_dist, ptrdiff_t n, float value) {
  ptrdiff_t lo = 0, hi = n;
  while (lo < hi) {
    ptrdiff_t mid = lo + (hi - lo) / 2;
    if (cum_dist[mid] <= value) lo = mid + 1;
    else hi = mid;
  }
  return lo;
}

TOPOTOOLBOX_API
void swath_longitudinal(
    float *restrict point_means, float *restrict point_stddevs,
    float *restrict point_mins, float *restrict point_maxs,
    ptrdiff_t *restrict point_counts, float *restrict point_medians,
    float *restrict point_q1, float *restrict point_q3,
    const int *restrict percentile_list, ptrdiff_t n_percentiles,
    float *restrict point_percentiles, const float *restrict dem,
    const float *restrict track_i, const float *restrict track_j,
    ptrdiff_t n_track_points, ptrdiff_t dims[2], float cellsize,
    float half_width, float binning_distance,
    int exclude_extended_bin) {
  if (n_track_points < 2)
    return;

  ptrdiff_t total_pixels = dims[0] * dims[1];
  float hw_px = half_width / cellsize;

  int compute_percentiles =
      (point_medians != NULL || point_q1 != NULL || point_q3 != NULL ||
       (percentile_list != NULL && n_percentiles > 0 &&
        point_percentiles != NULL));

  // Step 1: Global frontier distance map bounded by half_width
  // Need best_abs + signed_dist for tree march, nearest_seg for exclude check
  float *best_abs = (float *)malloc(total_pixels * sizeof(float));
  float *signed_dist_map = (float *)malloc(total_pixels * sizeof(float));
  ptrdiff_t *nearest_seg = (ptrdiff_t *)malloc(total_pixels * sizeof(ptrdiff_t));
  frontier_distance_map(best_abs, signed_dist_map, nearest_seg,
                        track_i, track_j, n_track_points, dims, hw_px);

  // Step 2: Cumulative along-track distance (meters)
  float *cum_dist = (float *)malloc(n_track_points * sizeof(float));
  cum_dist[0] = 0.0f;
  for (ptrdiff_t k = 1; k < n_track_points; k++) {
    float di = track_i[k] - track_i[k - 1];
    float dj = track_j[k] - track_j[k - 1];
    cum_dist[k] = cum_dist[k - 1] + sqrtf(di * di + dj * dj) * cellsize;
  }

  // Step 3: Initialize per-point accumulators
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

  // Step 4: Per-bin tree march on the global distance map
  // Work arrays reused across track points
  char *visited = (char *)calloc(total_pixels, sizeof(char));
  ptrdiff_t pixel_cap = 4096;
  ptrdiff_t *pixel_list = (ptrdiff_t *)malloc(pixel_cap * sizeof(ptrdiff_t));
  ptrdiff_t seed_cap = 1024;
  ptrdiff_t *seed_buf = (ptrdiff_t *)malloc(seed_cap * sizeof(ptrdiff_t));

  for (ptrdiff_t pt = 0; pt < n_track_points; pt++) {
    // Determine sub-track segment range for this bin
    ptrdiff_t seg_start, seg_end;
    if (binning_distance <= 0.0f) {
      seg_start = pt > 0 ? pt - 1 : 0;
      seg_end = pt < n_track_points - 1 ? pt : n_track_points - 2;
    } else {
      float range_lo = cum_dist[pt] - binning_distance;
      float range_hi = cum_dist[pt] + binning_distance;
      ptrdiff_t lb = lower_bound(cum_dist, n_track_points, range_lo);
      seg_start = lb > 0 ? lb - 1 : 0;
      ptrdiff_t ub = upper_bound(cum_dist, n_track_points, range_hi);
      seg_end = ub - 1;
      if (seg_end > n_track_points - 2) seg_end = n_track_points - 2;
      if (seg_end < seg_start) seg_end = seg_start;
    }

    // Rasterize sub-track segments as seeds
    ptrdiff_t n_seeds = 0;
    for (ptrdiff_t k = seg_start; k <= seg_end; k++) {
      float si = track_i[k + 1] - track_i[k];
      float sj = track_j[k + 1] - track_j[k];
      float seg_len = sqrtf(si * si + sj * sj);
      ptrdiff_t n_steps = (ptrdiff_t)(seg_len * 2.0f) + 1;

      for (ptrdiff_t s = 0; s <= n_steps; s++) {
        float t = n_steps > 0 ? (float)s / (float)n_steps : 0.0f;

        // exclude_extended_bin: skip seeds at full-track endpoints
        if (exclude_extended_bin) {
          if (k == 0 && t < 1e-3f) continue;
          if (k == n_track_points - 2 && t > (1.0f - 1e-3f)) continue;
        }

        ptrdiff_t pi = (ptrdiff_t)(track_i[k] + t * si + 0.5f);
        ptrdiff_t pj = (ptrdiff_t)(track_j[k] + t * sj + 0.5f);
        if (pi < 0 || pi >= dims[0] || pj < 0 || pj >= dims[1]) continue;

        ptrdiff_t idx = pj * dims[0] + pi;
        // Avoid duplicate seeds
        int dup = 0;
        for (ptrdiff_t d = 0; d < n_seeds; d++) {
          if (seed_buf[d] == idx) { dup = 1; break; }
        }
        if (dup) continue;

        if (n_seeds >= seed_cap) {
          seed_cap *= 2;
          seed_buf = (ptrdiff_t *)realloc(seed_buf, seed_cap * sizeof(ptrdiff_t));
        }
        seed_buf[n_seeds++] = idx;
      }
    }

    if (exclude_extended_bin) {
      // Tree march from seeds on the global distance map (no end caps)
      ptrdiff_t n_pixels = tree_march_fill(
          &pixel_list, &pixel_cap,
          best_abs, signed_dist_map,
          seed_buf, n_seeds,
          dims, hw_px, visited);

      for (ptrdiff_t p = 0; p < n_pixels; p++) {
        ptrdiff_t idx = pixel_list[p];
        if (isnan(dem[idx])) continue;
        accumulator_add(&accumulators[pt], dem[idx]);
        if (compute_percentiles) {
          percentile_accumulator_add(&p_accumulators[pt], dem[idx]);
        }
      }
    } else {
      // Full bean: all pixels within hw_px of the sub-track segments
      float bb_i0 = track_i[seg_start], bb_i1 = track_i[seg_start];
      float bb_j0 = track_j[seg_start], bb_j1 = track_j[seg_start];
      for (ptrdiff_t k = seg_start; k <= seg_end + 1; k++) {
        if (track_i[k] < bb_i0) bb_i0 = track_i[k];
        if (track_i[k] > bb_i1) bb_i1 = track_i[k];
        if (track_j[k] < bb_j0) bb_j0 = track_j[k];
        if (track_j[k] > bb_j1) bb_j1 = track_j[k];
      }
      ptrdiff_t i0 = (ptrdiff_t)(bb_i0 - hw_px - 1.0f);
      ptrdiff_t i1 = (ptrdiff_t)(bb_i1 + hw_px + 2.0f);
      ptrdiff_t j0 = (ptrdiff_t)(bb_j0 - hw_px - 1.0f);
      ptrdiff_t j1 = (ptrdiff_t)(bb_j1 + hw_px + 2.0f);
      if (i0 < 0) i0 = 0;
      if (j0 < 0) j0 = 0;
      if (i1 >= dims[0]) i1 = dims[0] - 1;
      if (j1 >= dims[1]) j1 = dims[1] - 1;

      for (ptrdiff_t pj = j0; pj <= j1; pj++) {
        for (ptrdiff_t pi = i0; pi <= i1; pi++) {
          float px = (float)pi;
          float py = (float)pj;
          float min_d = FLT_MAX;

          for (ptrdiff_t k = seg_start; k <= seg_end; k++) {
            float proj_i, proj_j, lam;
            float d = point_to_segment_distance(
                px, py, track_i[k], track_j[k],
                track_i[k + 1], track_j[k + 1],
                &proj_i, &proj_j, &lam);
            float ad = fabsf(d);
            if (ad < min_d) min_d = ad;
          }

          if (min_d <= hw_px) {
            ptrdiff_t idx = pj * dims[0] + pi;
            if (isnan(dem[idx])) continue;
            accumulator_add(&accumulators[pt], dem[idx]);
            if (compute_percentiles) {
              percentile_accumulator_add(&p_accumulators[pt], dem[idx]);
            }
          }
        }
      }
    }
  }

  free(visited);
  free(pixel_list);
  free(seed_buf);

  // Step 5: Finalize statistics
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
  free(best_abs);
  free(signed_dist_map);
  free(nearest_seg);
  free(cum_dist);
  free(accumulators);
  if (compute_percentiles) {
    for (ptrdiff_t k = 0; k < n_track_points; k++)
      percentile_accumulator_free(&p_accumulators[k]);
    free(p_accumulators);
  }
}

// ============================================================================
// Per-Point Pixel Retrieval
// ============================================================================
//
// Returns pixel coordinates associated with a single track point using
// the sub-track + perpendicular distance approach with bounding box.
// This function is called per-point and does not use the frontier.

TOPOTOOLBOX_API
ptrdiff_t swath_get_point_pixels(
    ptrdiff_t *restrict pixels_i, ptrdiff_t *restrict pixels_j,
    const float *restrict track_i, const float *restrict track_j,
    ptrdiff_t n_track_points, ptrdiff_t point_index, ptrdiff_t dims[2],
    float cellsize, float half_width, float binning_distance,
    int exclude_extended_bin) {
  if (n_track_points < 2)
    return 0;
  if (point_index < 0 || point_index >= n_track_points)
    return 0;

  ptrdiff_t total_pixels = dims[0] * dims[1];
  float hw_px = half_width / cellsize;

  // Global distance map (nearest_seg is required for correct Dijkstra expansion)
  float *best_abs = (float *)malloc(total_pixels * sizeof(float));
  float *signed_dist_map = (float *)malloc(total_pixels * sizeof(float));
  ptrdiff_t *nearest_seg = (ptrdiff_t *)malloc(total_pixels * sizeof(ptrdiff_t));
  frontier_distance_map(best_abs, signed_dist_map, nearest_seg,
                        track_i, track_j, n_track_points, dims, hw_px);

  // Cumulative along-track distance for segment range
  float *cum_dist = (float *)malloc(n_track_points * sizeof(float));
  cum_dist[0] = 0.0f;
  for (ptrdiff_t k = 1; k < n_track_points; k++) {
    float di = track_i[k] - track_i[k - 1];
    float dj = track_j[k] - track_j[k - 1];
    cum_dist[k] = cum_dist[k - 1] + sqrtf(di * di + dj * dj) * cellsize;
  }

  // Determine sub-track segment range
  ptrdiff_t seg_start, seg_end;
  if (binning_distance <= 0.0f) {
    seg_start = point_index > 0 ? point_index - 1 : 0;
    seg_end = point_index < n_track_points - 1
                  ? point_index : n_track_points - 2;
  } else {
    float range_lo = cum_dist[point_index] - binning_distance;
    float range_hi = cum_dist[point_index] + binning_distance;
    ptrdiff_t lb = lower_bound(cum_dist, n_track_points, range_lo);
    seg_start = lb > 0 ? lb - 1 : 0;
    ptrdiff_t ub = upper_bound(cum_dist, n_track_points, range_hi);
    seg_end = ub - 1;
    if (seg_end > n_track_points - 2) seg_end = n_track_points - 2;
    if (seg_end < seg_start) seg_end = seg_start;
  }

  // Rasterize sub-track as seeds
  ptrdiff_t seed_cap = 1024;
  ptrdiff_t n_seeds = 0;
  ptrdiff_t *seed_buf = (ptrdiff_t *)malloc(seed_cap * sizeof(ptrdiff_t));

  for (ptrdiff_t k = seg_start; k <= seg_end; k++) {
    float si = track_i[k + 1] - track_i[k];
    float sj = track_j[k + 1] - track_j[k];
    float seg_len = sqrtf(si * si + sj * sj);
    ptrdiff_t n_steps = (ptrdiff_t)(seg_len * 2.0f) + 1;

    for (ptrdiff_t s = 0; s <= n_steps; s++) {
      float t = n_steps > 0 ? (float)s / (float)n_steps : 0.0f;

      if (exclude_extended_bin) {
        if (k == 0 && t < 1e-3f) continue;
        if (k == n_track_points - 2 && t > (1.0f - 1e-3f)) continue;
      }

      ptrdiff_t pi = (ptrdiff_t)(track_i[k] + t * si + 0.5f);
      ptrdiff_t pj = (ptrdiff_t)(track_j[k] + t * sj + 0.5f);
      if (pi < 0 || pi >= dims[0] || pj < 0 || pj >= dims[1]) continue;

      ptrdiff_t idx = pj * dims[0] + pi;
      int dup = 0;
      for (ptrdiff_t d = 0; d < n_seeds; d++) {
        if (seed_buf[d] == idx) { dup = 1; break; }
      }
      if (dup) continue;

      if (n_seeds >= seed_cap) {
        seed_cap *= 2;
        seed_buf = (ptrdiff_t *)realloc(seed_buf, seed_cap * sizeof(ptrdiff_t));
      }
      seed_buf[n_seeds++] = idx;
    }
  }

  ptrdiff_t n_pixels = 0;

  if (exclude_extended_bin) {
    // Tree march on global distance map (no end caps)
    char *visited = (char *)calloc(total_pixels, sizeof(char));
    ptrdiff_t pixel_cap = 4096;
    ptrdiff_t *pixel_list = (ptrdiff_t *)malloc(pixel_cap * sizeof(ptrdiff_t));

    n_pixels = tree_march_fill(
        &pixel_list, &pixel_cap,
        best_abs, signed_dist_map,
        seed_buf, n_seeds,
        dims, hw_px, visited);

    for (ptrdiff_t p = 0; p < n_pixels; p++) {
      ptrdiff_t idx = pixel_list[p];
      pixels_i[p] = idx % dims[0];
      pixels_j[p] = idx / dims[0];
    }

    free(pixel_list);
    free(visited);
  } else {
    // Full bean: all pixels within hw_px of the sub-track segments
    // Compute bounding box of sub-track ± hw_px
    float bb_i0 = track_i[seg_start], bb_i1 = track_i[seg_start];
    float bb_j0 = track_j[seg_start], bb_j1 = track_j[seg_start];
    for (ptrdiff_t k = seg_start; k <= seg_end + 1; k++) {
      if (track_i[k] < bb_i0) bb_i0 = track_i[k];
      if (track_i[k] > bb_i1) bb_i1 = track_i[k];
      if (track_j[k] < bb_j0) bb_j0 = track_j[k];
      if (track_j[k] > bb_j1) bb_j1 = track_j[k];
    }
    ptrdiff_t i0 = (ptrdiff_t)(bb_i0 - hw_px - 1.0f);
    ptrdiff_t i1 = (ptrdiff_t)(bb_i1 + hw_px + 2.0f);
    ptrdiff_t j0 = (ptrdiff_t)(bb_j0 - hw_px - 1.0f);
    ptrdiff_t j1 = (ptrdiff_t)(bb_j1 + hw_px + 2.0f);
    if (i0 < 0) i0 = 0;
    if (j0 < 0) j0 = 0;
    if (i1 >= dims[0]) i1 = dims[0] - 1;
    if (j1 >= dims[1]) j1 = dims[1] - 1;

    // For each pixel in the bounding box, check distance to sub-track
    for (ptrdiff_t pj = j0; pj <= j1; pj++) {
      for (ptrdiff_t pi = i0; pi <= i1; pi++) {
        float px = (float)pi;
        float py = (float)pj;
        float min_d = FLT_MAX;

        for (ptrdiff_t k = seg_start; k <= seg_end; k++) {
          float proj_i, proj_j, lam;
          float d = point_to_segment_distance(
              px, py, track_i[k], track_j[k],
              track_i[k + 1], track_j[k + 1],
              &proj_i, &proj_j, &lam);
          float ad = fabsf(d);
          if (ad < min_d) min_d = ad;
        }

        if (min_d <= hw_px) {
          pixels_i[n_pixels] = pi;
          pixels_j[n_pixels] = pj;
          n_pixels++;
        }
      }
    }
  }

  free(seed_buf);
  free(best_abs);
  free(signed_dist_map);
  free(nearest_seg);
  free(cum_dist);
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

TOPOTOOLBOX_API
ptrdiff_t swath_invert_distance_map(
    float *restrict centre_line_i,
    float *restrict centre_line_j,
    float *restrict width,
    float *restrict dist_from_boundary_out,
    const float *restrict track_i,
    const float *restrict track_j,
    ptrdiff_t n_track_points,
    ptrdiff_t dims[2],
    float cellsize,
    float half_width)
{
  if (n_track_points < 2)
    return 0;

  ptrdiff_t total = dims[0] * dims[1];
  float hw_px = half_width / cellsize;

  // Step 1: Global distance map
  float *best_abs = (float *)malloc(total * sizeof(float));
  float *signed_dist = (float *)malloc(total * sizeof(float));
  frontier_distance_map(best_abs, signed_dist, NULL,
                        track_i, track_j, n_track_points, dims, hw_px);

  // Step 2: Find boundary pixels and split by side
  // A boundary pixel is inside the mask but has a D4 neighbour outside
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

    if (signed_dist[idx] >= 0.0f) {
      if (n_pos >= pos_cap) {
        pos_cap *= 2;
        pos_seeds = (ptrdiff_t *)realloc(pos_seeds, pos_cap * sizeof(ptrdiff_t));
      }
      pos_seeds[n_pos++] = idx;
    }
    if (signed_dist[idx] <= 0.0f) {
      if (n_neg >= neg_cap) {
        neg_cap *= 2;
        neg_seeds = (ptrdiff_t *)realloc(neg_seeds, neg_cap * sizeof(ptrdiff_t));
      }
      neg_seeds[n_neg++] = idx;
    }
  }

  // Step 3: Dijkstra from positive-side boundary
  float *dist_pos = (float *)malloc(total * sizeof(float));
  boundary_dijkstra(dist_pos, best_abs, pos_seeds, n_pos, dims, hw_px);

  // Step 4: Dijkstra from negative-side boundary
  float *dist_neg = (float *)malloc(total * sizeof(float));
  boundary_dijkstra(dist_neg, best_abs, neg_seeds, n_neg, dims, hw_px);

  free(pos_seeds);
  free(neg_seeds);

  // Step 5: Compute combined distance-from-boundary = min(dist_pos, dist_neg)
  // and extract centre line via ridge detection
  float *dfb = dist_from_boundary_out;
  int free_dfb = 0;
  if (dfb == NULL) {
    dfb = (float *)malloc(total * sizeof(float));
    free_dfb = 1;
  }

  for (ptrdiff_t idx = 0; idx < total; idx++) {
    dfb[idx] = dist_pos[idx] < dist_neg[idx] ? dist_pos[idx] : dist_neg[idx];
  }

  // Centre line: ridge of dfb — local maxima in at least 2 opposing D4
  // directions (i.e., value >= both neighbours along at least one axis)
  ptrdiff_t count = 0;
  for (ptrdiff_t idx = 0; idx < total; idx++) {
    if (best_abs[idx] == FLT_MAX || best_abs[idx] > hw_px)
      continue;
    if (dfb[idx] == FLT_MAX || dfb[idx] == 0.0f)
      continue;

    ptrdiff_t ci = idx % dims[0];
    ptrdiff_t cj = idx / dims[0];
    float val = dfb[idx];

    // Check if ridge along i-axis (left-right)
    int ridge_i = 0;
    if (ci > 0 && ci < dims[0] - 1) {
      float left = dfb[idx - 1];
      float right = dfb[idx + 1];
      if (val >= left && val >= right)
        ridge_i = 1;
    }

    // Check if ridge along j-axis (up-down)
    int ridge_j = 0;
    if (cj > 0 && cj < dims[1] - 1) {
      float up = dfb[idx - dims[0]];
      float down = dfb[idx + dims[0]];
      if (val >= up && val >= down)
        ridge_j = 1;
    }

    if (!ridge_i && !ridge_j)
      continue;

    // This pixel is on the centre line
    centre_line_i[count] = (float)ci;
    centre_line_j[count] = (float)cj;
    width[count] = (dist_pos[idx] + dist_neg[idx]) * cellsize;
    count++;
  }

  // Cleanup
  free(dist_pos);
  free(dist_neg);
  free(best_abs);
  free(signed_dist);
  if (free_dfb) free(dfb);

  return count;
}
