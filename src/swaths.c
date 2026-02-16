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

  // Step 1: Frontier distance map bounded by half_width
  float *best_abs = (float *)malloc(total_pixels * sizeof(float));
  ptrdiff_t *nearest_seg = (ptrdiff_t *)malloc(total_pixels * sizeof(ptrdiff_t));
  frontier_distance_map(best_abs, NULL, nearest_seg,
                        track_i, track_j, n_track_points, dims, hw_px);

  // Step 2: Cumulative along-track distance (meters)
  float *cum_dist = (float *)malloc(n_track_points * sizeof(float));
  cum_dist[0] = 0.0f;
  for (ptrdiff_t k = 1; k < n_track_points; k++) {
    float di = track_i[k] - track_i[k - 1];
    float dj = track_j[k] - track_j[k - 1];
    cum_dist[k] = cum_dist[k - 1] + sqrtf(di * di + dj * dj) * cellsize;
  }
  float total_track_length = cum_dist[n_track_points - 1];

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

  // Step 4: Single pass over all pixels
  for (ptrdiff_t j = 0; j < dims[1]; j++) {
    for (ptrdiff_t i = 0; i < dims[0]; i++) {
      ptrdiff_t idx = j * dims[0] + i;

      if (best_abs[idx] > hw_px) continue;  // unvisited or outside
      if (isnan(dem[idx])) continue;

      ptrdiff_t seg = nearest_seg[idx];
      if (seg < 0 || seg >= n_track_points - 1) continue;

      // Recompute lambda for this pixel's known nearest segment (cheap: 1 segment)
      float px = (float)i;
      float py = (float)j;
      float proj_x, proj_y, lambda;
      point_to_segment_distance(px, py,
          track_i[seg], track_j[seg], track_i[seg + 1], track_j[seg + 1],
          &proj_x, &proj_y, &lambda);

      // Along-track position interpolated from segment endpoints
      float along_pos = cum_dist[seg] +
                         lambda * (cum_dist[seg + 1] - cum_dist[seg]);

      // exclude_extended_bin: clamp pixels projecting past full-track endpoints.
      // Sub-track endpoint clipping is handled naturally by the along-track
      // range check (pixels past sub-track boundaries get correct along_pos
      // from the full-track distance map).
      if (exclude_extended_bin) {
        if (seg == 0 && lambda < 1e-6f)
          along_pos = 0.0f;
        if (seg == n_track_points - 2 && lambda > (1.0f - 1e-6f))
          along_pos = total_track_length;
      }

      // Binary search: find track points within ±binning_distance
      ptrdiff_t k_start = lower_bound(cum_dist, n_track_points,
                                       along_pos - binning_distance);
      ptrdiff_t k_end = upper_bound(cum_dist, n_track_points,
                                     along_pos + binning_distance);

      float value = dem[idx];
      for (ptrdiff_t k = k_start; k < k_end; k++) {
        accumulator_add(&accumulators[k], value);
        if (compute_percentiles) {
          percentile_accumulator_add(&p_accumulators[k], value);
        }
      }
    }
  }

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

  // Compute cumulative along-track distance (meters)
  float *cum_dist = (float *)malloc(n_track_points * sizeof(float));
  cum_dist[0] = 0.0f;
  for (ptrdiff_t k = 1; k < n_track_points; k++) {
    float di = track_i[k] - track_i[k - 1];
    float dj = track_j[k] - track_j[k - 1];
    cum_dist[k] = cum_dist[k - 1] + sqrtf(di * di + dj * dj) * cellsize;
  }

  // Find sub-track range
  float center_dist = cum_dist[point_index];
  ptrdiff_t sub_start = point_index;
  ptrdiff_t sub_end = point_index;

  while (sub_start > 0 &&
         (center_dist - cum_dist[sub_start - 1]) <= binning_distance)
    sub_start--;
  while (sub_end < n_track_points - 1 &&
         (cum_dist[sub_end + 1] - center_dist) <= binning_distance)
    sub_end++;

  if (sub_start == sub_end) {
    if (sub_end < n_track_points - 1) sub_end++;
    else if (sub_start > 0) sub_start--;
  }

  free(cum_dist);

  ptrdiff_t n_sub_segments = sub_end - sub_start;
  if (n_sub_segments < 1)
    return 0;

  // Bounding box of sub-track ± half_width (pixels)
  float hw_pixels = half_width / cellsize;
  float bb_min_i = track_i[sub_start];
  float bb_max_i = track_i[sub_start];
  float bb_min_j = track_j[sub_start];
  float bb_max_j = track_j[sub_start];
  for (ptrdiff_t s = sub_start + 1; s <= sub_end; s++) {
    if (track_i[s] < bb_min_i) bb_min_i = track_i[s];
    if (track_i[s] > bb_max_i) bb_max_i = track_i[s];
    if (track_j[s] < bb_min_j) bb_min_j = track_j[s];
    if (track_j[s] > bb_max_j) bb_max_j = track_j[s];
  }

  ptrdiff_t i_lo = (ptrdiff_t)(bb_min_i - hw_pixels - 1.0f);
  ptrdiff_t i_hi = (ptrdiff_t)(bb_max_i + hw_pixels + 2.0f);
  ptrdiff_t j_lo = (ptrdiff_t)(bb_min_j - hw_pixels - 1.0f);
  ptrdiff_t j_hi = (ptrdiff_t)(bb_max_j + hw_pixels + 2.0f);

  if (i_lo < 0) i_lo = 0;
  if (i_hi > dims[0]) i_hi = dims[0];
  if (j_lo < 0) j_lo = 0;
  if (j_hi > dims[1]) j_hi = dims[1];

  ptrdiff_t n = 0;
  for (ptrdiff_t j = j_lo; j < j_hi; j++) {
    for (ptrdiff_t i = i_lo; i < i_hi; i++) {
      float px = (float)i;
      float py = (float)j;

      float min_dist, proj_x, proj_y, lambda;
      ptrdiff_t seg = find_nearest_segment(
          px, py, track_i, track_j, n_track_points,
          sub_start, sub_end, &min_dist, &proj_x, &proj_y, &lambda);

      if (seg < 0) continue;

      if (exclude_extended_bin) {
        if (seg == sub_start && lambda < 1e-6f && sub_start > 0)
          continue;
        if (seg == sub_end - 1 && lambda > (1.0f - 1e-6f) &&
            sub_end < n_track_points - 1)
          continue;
      }

      float dist_m = fabsf(min_dist) * cellsize;
      if (dist_m > half_width) continue;

      pixels_i[n] = i;
      pixels_j[n] = j;
      n++;
    }
  }

  return n;
}
