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

  1. Bin-based: Aggregates DEM pixels by perpendicular distance to track,
     creating a single averaged cross-sectional profile.

  2. Per-point: Computes individual swath statistics for each track point,
     showing how the profile changes along the track.

  Both methods use Euclidean distance in pixel space scaled by cellsize.
*/

// ============================================================================
// Statistics Accumulator - Welford's online algorithm for numerical stability
// ============================================================================
// Accumulates mean, variance, min, max in a single pass without storing values.
// Uses double precision internally to prevent catastrophic cancellation in
// variance computation even with large elevation values.

typedef struct {
  ptrdiff_t count;
  double sum;       // For mean computation
  double sum_sq;    // For variance: E[X²] - E[X]²
  float min_val;
  float max_val;
} swath_stats_accumulator;

// ============================================================================
// Percentile Accumulator - Stores values for percentile computation
// ============================================================================
// Dynamically stores all values to enable median and percentile calculation.
// Memory grows as needed during accumulation.

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

  // Grow array if needed
  if (acc->count >= acc->capacity) {
    acc->capacity *= 2;
    acc->values = (float *)realloc(acc->values, acc->capacity * sizeof(float));
  }

  acc->values[acc->count++] = value;
}

// Comparison function for qsort
static int compare_floats(const void *a, const void *b) {
  float fa = *(const float *)a;
  float fb = *(const float *)b;
  return (fa > fb) - (fa < fb);
}

// Sort values in preparation for percentile computation
static inline void percentile_accumulator_sort(percentile_accumulator *acc) {
  if (acc->count > 0) {
    qsort(acc->values, acc->count, sizeof(float), compare_floats);
  }
}

// Get percentile from sorted values (p in [0, 100])
// Uses linear interpolation between nearest ranks
static inline float percentile_accumulator_get(const percentile_accumulator *acc,
                                                float p) {
  if (acc->count == 0)
    return NAN;

  // Clamp percentile to valid range
  if (p < 0.0f)
    p = 0.0f;
  if (p > 100.0f)
    p = 100.0f;

  // Convert percentile to array index (0-based)
  float index = p / 100.0f * (acc->count - 1);
  ptrdiff_t lower = (ptrdiff_t)index;
  ptrdiff_t upper = lower + 1;

  // Boundary cases
  if (upper >= acc->count)
    return acc->values[acc->count - 1];

  // Linear interpolation between adjacent values
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
  if (value < acc->min_val)
    acc->min_val = value;
  if (value > acc->max_val)
    acc->max_val = value;
}

static inline float accumulator_mean(const swath_stats_accumulator *acc) {
  return acc->count > 0 ? (float)(acc->sum / acc->count) : NAN;
}

// Computes unbiased sample standard deviation (n-1 denominator)
static inline float accumulator_stddev(const swath_stats_accumulator *acc) {
  if (acc->count <= 1)
    return 0.0f;
  double variance = (acc->sum_sq - acc->sum * acc->sum / acc->count) /
                    (acc->count - 1);
  return variance > 0.0 ? (float)sqrt(variance) : 0.0f;
}

// ============================================================================
// Geometric Distance Computation
// ============================================================================

/*
  Compute signed perpendicular distance from point to line segment.

  Uses vector projection: projects point onto line, clamps to segment endpoints,
  then computes Euclidean distance. Sign determined by cross product (2D rotation):
  - Negative: point is left of segment direction (a→b)
  - Positive: point is right of segment direction

  Returns distance in same units as input coordinates (typically pixel space).
*/
static float point_to_segment_distance(float px, float py, float ax, float ay,
                                       float bx, float by, float *proj_x,
                                       float *proj_y, float *lambda) {
  // Segment direction vector
  float dx = bx - ax;
  float dy = by - ay;

  // Vector from segment start to point
  float dpx = px - ax;
  float dpy = py - ay;

  // Projection parameter: t = (p-a)·(b-a) / |b-a|²
  float seg_length_sq = dx * dx + dy * dy;

  if (seg_length_sq < 1e-10f) {
    // Degenerate segment (endpoints coincide), treat as point distance
    *lambda = 0.0f;
    *proj_x = ax;
    *proj_y = ay;
    float dist = sqrtf(dpx * dpx + dpy * dpy);
    return dist;
  }

  float t = (dpx * dx + dpy * dy) / seg_length_sq;

  // Clamp projection to segment: t ∈ [0,1]
  // t < 0: closest point is 'a'
  // t > 1: closest point is 'b'
  // else:  closest point is on segment interior
  if (t < 0.0f)
    t = 0.0f;
  else if (t > 1.0f)
    t = 1.0f;

  *lambda = t;

  // Compute projection point on segment
  *proj_x = ax + t * dx;
  *proj_y = ay + t * dy;

  // Distance from projection to point
  float to_p_x = px - *proj_x;
  float to_p_y = py - *proj_y;
  float dist = sqrtf(to_p_x * to_p_x + to_p_y * to_p_y);

  // Determine sign via cross product: (p-a) × (b-a)
  // Positive cross product → point is to the right
  float cross = dpx * dy - dpy * dx;

  return cross >= 0.0f ? dist : -dist;
}

/*
  Find nearest track segment to a point within search range.

  Linear search through segment range, comparing absolute distances to handle
  signed distances correctly. Uses projection information for downstream
  interpolation if needed.

  Search range optimization: Typically called with hints from neighbors,
  narrowing search to nearby segments rather than entire track.
*/
static ptrdiff_t find_nearest_segment(float px, float py,
                                      const float *restrict track_i,
                                      const float *restrict track_j,
                                      ptrdiff_t n_track_points,
                                      ptrdiff_t search_start,
                                      ptrdiff_t search_end, float *min_dist,
                                      float *proj_x, float *proj_y) {
  ptrdiff_t best_segment = -1;
  *min_dist = FLT_MAX;

  // Clamp search range to valid segment indices
  if (search_start < 0)
    search_start = 0;
  if (search_end >= n_track_points)
    search_end = n_track_points - 1;

  // Test each segment in range
  for (ptrdiff_t k = search_start; k < search_end; k++) {
    float ax = track_i[k];
    float ay = track_j[k];
    float bx = track_i[k + 1];
    float by = track_j[k + 1];

    float proj_i, proj_j, lambda;
    float dist =
        point_to_segment_distance(px, py, ax, ay, bx, by, &proj_i, &proj_j, &lambda);

    // Compare absolute distances since sign only indicates direction
    if (fabsf(dist) < fabsf(*min_dist)) {
      *min_dist = dist;
      best_segment = k;
      *proj_x = proj_i;
      *proj_y = proj_j;
    }
  }

  return best_segment;
}

// ============================================================================
// Public API Functions
// ============================================================================

TOPOTOOLBOX_API
ptrdiff_t swath_compute_nbins(float half_width, float bin_resolution) {
  if (bin_resolution <= 0.0f || half_width <= 0.0f)
    return 0;
  // Bins span [-half_width, +half_width] with given resolution
  // Formula: 2 * half_width / resolution, rounded to nearest odd number
  // +1 ensures odd count so bin centers align (one bin exactly at distance=0)
  return (ptrdiff_t)(2.0f * half_width / bin_resolution) + 1;
}

// ============================================================================
// Distance Map Computation with Iterative Refinement
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

  // Allocate temporary segment map if not provided by caller
  // This map stores which track segment is nearest to each pixel
  ptrdiff_t total_pixels = dims[0] * dims[1];
  ptrdiff_t *temp_nearest = NULL;
  int need_free = 0;

  if (nearest_segment == NULL) {
    temp_nearest = (ptrdiff_t *)malloc(total_pixels * sizeof(ptrdiff_t));
    need_free = 1;
  } else {
    temp_nearest = nearest_segment;
  }

  // Initialize: -1 indicates no segment assigned yet
  for (ptrdiff_t idx = 0; idx < total_pixels; idx++) {
    temp_nearest[idx] = -1;
  }

  // Seed: assign track pixels to their own segments
  // This provides initial "islands" from which refinement propagates
  for (ptrdiff_t k = 0; k < n_track_points - 1; k++) {
    ptrdiff_t i = (ptrdiff_t)(track_i[k] + 0.5f);
    ptrdiff_t j = (ptrdiff_t)(track_j[k] + 0.5f);

    if (i >= 0 && i < dims[0] && j >= 0 && j < dims[1]) {
      temp_nearest[j * dims[0] + i] = k;
    }
  }

  // Iterative refinement: propagate segment assignments outward
  // Similar to distance transform algorithms, but handles arbitrary track geometry
  // Typically converges in 5-10 iterations for reasonable track spacing
  int max_iterations = 20;
  for (int iter = 0; iter < max_iterations; iter++) {
    int changes = 0;

    // Forward pass: top-left to bottom-right
    for (ptrdiff_t j = 0; j < dims[1]; j++) {
      for (ptrdiff_t i = 0; i < dims[0]; i++) {
        ptrdiff_t idx = j * dims[0] + i;

        ptrdiff_t best_seg = temp_nearest[idx];
        ptrdiff_t search_start = 0;
        ptrdiff_t search_end = n_track_points - 1;

        // Optimization: use neighbor segments as hints to narrow search
        // If a neighbor knows its nearest segment, check nearby segments first
        for (int di = -1; di <= 1; di++) {
          for (int dj = -1; dj <= 1; dj++) {
            if (di == 0 && dj == 0)
              continue;

            ptrdiff_t ni = i + di;
            ptrdiff_t nj = j + dj;

            if (ni >= 0 && ni < dims[0] && nj >= 0 && nj < dims[1]) {
              ptrdiff_t neighbor_seg = temp_nearest[nj * dims[0] + ni];
              if (neighbor_seg >= 0) {
                if (best_seg < 0) {
                  // Neighbor has assignment, search ±2 segments around it
                  search_start = neighbor_seg > 1 ? neighbor_seg - 1 : 0;
                  search_end = neighbor_seg + 2 < n_track_points - 1
                                   ? neighbor_seg + 2
                                   : n_track_points - 1;
                  break;
                }
              }
            }
          }
          if (best_seg >= 0 || search_start > 0)
            break;
        }

        // Fallback: if no hints available, search all segments (rare after first iteration)
        if (best_seg < 0 && search_start == 0) {
          search_start = 0;
          search_end = n_track_points - 1;
        }

        // Work in pixel coordinates for distance computation
        float px = (float)i;
        float py = (float)j;
        float min_dist, proj_x, proj_y;

        ptrdiff_t new_seg = find_nearest_segment(
            px, py, track_i, track_j, n_track_points, search_start, search_end,
            &min_dist, &proj_x, &proj_y);

        // Update if we found a segment and it's different from current assignment
        if (new_seg >= 0 && new_seg != temp_nearest[idx]) {
          temp_nearest[idx] = new_seg;
          // Convert distance from pixels to meters
          distance[idx] = min_dist * cellsize;
          changes++;
        } else if (new_seg >= 0) {
          // Convert distance from pixels to meters
          distance[idx] = min_dist * cellsize;
        }
      }
    }

    // Convergence: no assignments changed, we're done
    if (changes == 0)
      break;
  }

  if (need_free)
    free(temp_nearest);
}

// ============================================================================
// Binned Swath Profile - Averaged Cross-Section
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

  // Determine if we need percentile computation
  int compute_percentiles =
      (bin_medians != NULL || bin_q1 != NULL || bin_q3 != NULL ||
       (percentile_list != NULL && n_percentiles > 0 && bin_percentiles != NULL));

  // Allocate one stats accumulator per distance bin
  swath_stats_accumulator *accumulators =
      (swath_stats_accumulator *)calloc(n_bins, sizeof(swath_stats_accumulator));

  for (ptrdiff_t b = 0; b < n_bins; b++) {
    accumulator_init(&accumulators[b]);
  }

  // Allocate percentile accumulators if needed
  percentile_accumulator *p_accumulators = NULL;
  if (compute_percentiles) {
    p_accumulators = (percentile_accumulator *)calloc(
        n_bins, sizeof(percentile_accumulator));
    for (ptrdiff_t b = 0; b < n_bins; b++) {
      percentile_accumulator_init(&p_accumulators[b], 64);
    }
  }

  // Compute signed perpendicular distance for every pixel
  float *distance = (float *)malloc(dims[0] * dims[1] * sizeof(float));
  swath_distance_map(distance, NULL, track_i, track_j, n_track_points, dims,
                     cellsize);

  // Normalization: compute mean elevation at track center
  // This removes absolute elevation, highlighting relative topography
  float reference_elevation = 0.0f;
  if (normalize) {
    swath_stats_accumulator track_acc;
    accumulator_init(&track_acc);

    // Sample pixels very close to track (within one bin width)
    for (ptrdiff_t j = 0; j < dims[1]; j++) {
      for (ptrdiff_t i = 0; i < dims[0]; i++) {
        ptrdiff_t idx = j * dims[0] + i;
        if (fabsf(distance[idx]) <= bin_resolution) {
          accumulator_add(&track_acc, dem[idx]);
        }
      }
    }
    reference_elevation = accumulator_mean(&track_acc);
  }

  // Bin pixels by perpendicular distance
  // All pixels along entire track length contribute to bins
  for (ptrdiff_t j = 0; j < dims[1]; j++) {
    for (ptrdiff_t i = 0; i < dims[0]; i++) {
      ptrdiff_t idx = j * dims[0] + i;
      float dist = distance[idx];

      // Exclude pixels outside swath width
      if (fabsf(dist) > half_width)
        continue;
      if (isnan(dem[idx]))
        continue;

      // Map distance to bin index
      // Distance range [-half_width, +half_width] → bin index [0, n_bins-1]
      ptrdiff_t bin_idx =
          (ptrdiff_t)((dist + half_width) / bin_resolution + 0.5f);
      if (bin_idx < 0)
        bin_idx = 0;
      if (bin_idx >= n_bins)
        bin_idx = n_bins - 1;

      float value = dem[idx];
      if (normalize) {
        value -= reference_elevation;
      }

      accumulator_add(&accumulators[bin_idx], value);

      // Also store for percentile computation if needed
      if (compute_percentiles) {
        percentile_accumulator_add(&p_accumulators[bin_idx], value);
      }
    }
  }

  // Compute percentiles if requested
  if (compute_percentiles) {
    for (ptrdiff_t b = 0; b < n_bins; b++) {
      percentile_accumulator_sort(&p_accumulators[b]);

      // Compute standard percentiles
      if (bin_medians != NULL) {
        bin_medians[b] = percentile_accumulator_get(&p_accumulators[b], 50.0f);
        if (normalize && !isnan(bin_medians[b])) {
          bin_medians[b] += reference_elevation;
        }
      }

      if (bin_q1 != NULL) {
        bin_q1[b] = percentile_accumulator_get(&p_accumulators[b], 25.0f);
        if (normalize && !isnan(bin_q1[b])) {
          bin_q1[b] += reference_elevation;
        }
      }

      if (bin_q3 != NULL) {
        bin_q3[b] = percentile_accumulator_get(&p_accumulators[b], 75.0f);
        if (normalize && !isnan(bin_q3[b])) {
          bin_q3[b] += reference_elevation;
        }
      }

      // Compute arbitrary percentiles
      if (percentile_list != NULL && n_percentiles > 0 &&
          bin_percentiles != NULL) {
        for (ptrdiff_t p = 0; p < n_percentiles; p++) {
          float pval = percentile_accumulator_get(&p_accumulators[b],
                                                   (float)percentile_list[p]);
          if (normalize && !isnan(pval)) {
            pval += reference_elevation;
          }
          bin_percentiles[b * n_percentiles + p] = pval;
        }
      }
    }
  }

  // Finalize basic statistics for each bin
  for (ptrdiff_t b = 0; b < n_bins; b++) {
    bin_distances[b] = -half_width + b * bin_resolution;

    float mean = accumulator_mean(&accumulators[b]);
    float stddev = accumulator_stddev(&accumulators[b]);

    // Restore absolute elevation if normalized
    if (normalize && !isnan(mean)) {
      mean += reference_elevation;
    }

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
  free(distance);
  free(accumulators);
  if (compute_percentiles) {
    for (ptrdiff_t b = 0; b < n_bins; b++) {
      percentile_accumulator_free(&p_accumulators[b]);
    }
    free(p_accumulators);
  }
}

// ============================================================================
// Per-Point Swath Profile - Along-Track Variation
// ============================================================================
//
// For each track point k, defines a local sub-track of all track points
// within ±binning_distance (measured along the track) of point k. Then uses
// the same perpendicular-distance-to-segments logic as the transverse method
// to gather DEM pixels within half_width of that sub-track.
//
// If binning_distance < cellsize: the sub-track is just the two half-segments
// around point k, producing a single orthogonal line at that point.
//
// If binning_distance >= cellsize: the sub-track encompasses a longer portion
// of the track, gathering pixels as if that sub-track were the entire profile
// for the transverse method.

TOPOTOOLBOX_API
void swath_longitudinal(
    float *restrict point_means, float *restrict point_stddevs,
    float *restrict point_mins, float *restrict point_maxs,
    ptrdiff_t *restrict point_counts, float *restrict point_medians,
    float *restrict point_q1, float *restrict point_q3,
    const int *restrict percentile_list, ptrdiff_t n_percentiles,
    float *restrict point_percentiles, ptrdiff_t *restrict pixel_indices,
    ptrdiff_t *restrict point_offsets, const float *restrict dem,
    const float *restrict track_i, const float *restrict track_j,
    ptrdiff_t n_track_points, ptrdiff_t dims[2], float cellsize,
    float half_width, float binning_distance) {
  if (n_track_points < 2)
    return;

  int compute_percentiles =
      (point_medians != NULL || point_q1 != NULL || point_q3 != NULL ||
       (percentile_list != NULL && n_percentiles > 0 &&
        point_percentiles != NULL));

  int track_pixels = (pixel_indices != NULL && point_offsets != NULL);

  // Precompute cumulative along-track distance (in meters) for each point
  float *cum_dist = (float *)malloc(n_track_points * sizeof(float));
  cum_dist[0] = 0.0f;
  for (ptrdiff_t k = 1; k < n_track_points; k++) {
    float di = track_i[k] - track_i[k - 1];
    float dj = track_j[k] - track_j[k - 1];
    cum_dist[k] = cum_dist[k - 1] + sqrtf(di * di + dj * dj) * cellsize;
  }

  // Convert half_width to pixels for bounding box computation
  float hw_pixels = half_width / cellsize;

  // If tracking pixels: two-pass approach (count then store)
  // First pass counts pixels per point, second pass stores indices
  ptrdiff_t *pixel_counts_temp = NULL;
  if (track_pixels) {
    pixel_counts_temp = (ptrdiff_t *)calloc(n_track_points, sizeof(ptrdiff_t));
    point_offsets[0] = 0;
  }

  // Process each track point
  for (ptrdiff_t k = 0; k < n_track_points; k++) {
    // Find sub-track range: all points within ±binning_distance along track
    float center_dist = cum_dist[k];
    ptrdiff_t sub_start = k;
    ptrdiff_t sub_end = k;

    // Expand sub-track backwards
    while (sub_start > 0 &&
           (center_dist - cum_dist[sub_start - 1]) <= binning_distance) {
      sub_start--;
    }
    // Expand sub-track forwards
    while (sub_end < n_track_points - 1 &&
           (cum_dist[sub_end + 1] - center_dist) <= binning_distance) {
      sub_end++;
    }

    // Ensure at least one segment exists
    // If sub_start == sub_end and it's not the last point, use [k, k+1]
    // If it's the last point, use [k-1, k]
    if (sub_start == sub_end) {
      if (sub_end < n_track_points - 1)
        sub_end++;
      else if (sub_start > 0)
        sub_start--;
    }

    // Number of segments in sub-track
    ptrdiff_t n_sub_segments = sub_end - sub_start;
    if (n_sub_segments < 1) {
      // Degenerate: single point, no segments possible
      point_means[k] = NAN;
      point_stddevs[k] = 0.0f;
      point_mins[k] = NAN;
      point_maxs[k] = NAN;
      point_counts[k] = 0;
      if (point_medians) point_medians[k] = NAN;
      if (point_q1) point_q1[k] = NAN;
      if (point_q3) point_q3[k] = NAN;
      if (percentile_list && n_percentiles > 0 && point_percentiles) {
        for (ptrdiff_t p = 0; p < n_percentiles; p++)
          point_percentiles[k * n_percentiles + p] = NAN;
      }
      if (track_pixels) {
        pixel_counts_temp[k] = 0;
      }
      continue;
    }

    // Compute bounding box of sub-track ± half_width (in pixel coords)
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

    // Expand bounding box by half_width (in pixels) + 1 pixel margin
    ptrdiff_t i_lo = (ptrdiff_t)(bb_min_i - hw_pixels - 1.0f);
    ptrdiff_t i_hi = (ptrdiff_t)(bb_max_i + hw_pixels + 2.0f);
    ptrdiff_t j_lo = (ptrdiff_t)(bb_min_j - hw_pixels - 1.0f);
    ptrdiff_t j_hi = (ptrdiff_t)(bb_max_j + hw_pixels + 2.0f);

    // Clamp to grid bounds
    if (i_lo < 0) i_lo = 0;
    if (i_hi > dims[0]) i_hi = dims[0];
    if (j_lo < 0) j_lo = 0;
    if (j_hi > dims[1]) j_hi = dims[1];

    // Initialize accumulators for this point
    swath_stats_accumulator acc;
    accumulator_init(&acc);

    percentile_accumulator p_acc = {0, 0, NULL};
    if (compute_percentiles) {
      percentile_accumulator_init(&p_acc, 64);
    }

    // Scan pixels in bounding box, compute perpendicular distance to sub-track
    for (ptrdiff_t j = j_lo; j < j_hi; j++) {
      for (ptrdiff_t i = i_lo; i < i_hi; i++) {
        ptrdiff_t idx = j * dims[0] + i;
        if (isnan(dem[idx]))
          continue;

        float px = (float)i;
        float py = (float)j;

        // Find perpendicular distance to sub-track (search only sub-track segments)
        float min_dist, proj_x, proj_y;
        ptrdiff_t seg = find_nearest_segment(
            px, py, track_i, track_j, n_track_points,
            sub_start, sub_end, &min_dist, &proj_x, &proj_y);

        if (seg < 0)
          continue;

        // Convert distance from pixels to meters
        float dist_m = fabsf(min_dist) * cellsize;
        if (dist_m > half_width)
          continue;

        accumulator_add(&acc, dem[idx]);
        if (compute_percentiles) {
          percentile_accumulator_add(&p_acc, dem[idx]);
        }
        if (track_pixels) {
          pixel_counts_temp[k]++;
        }
      }
    }

    // Finalize statistics for this point
    point_means[k] = accumulator_mean(&acc);
    point_stddevs[k] = accumulator_stddev(&acc);
    point_mins[k] = acc.count > 0 ? acc.min_val : NAN;
    point_maxs[k] = acc.count > 0 ? acc.max_val : NAN;
    point_counts[k] = acc.count;

    // Percentiles
    if (compute_percentiles) {
      percentile_accumulator_sort(&p_acc);

      if (point_medians != NULL)
        point_medians[k] = percentile_accumulator_get(&p_acc, 50.0f);
      if (point_q1 != NULL)
        point_q1[k] = percentile_accumulator_get(&p_acc, 25.0f);
      if (point_q3 != NULL)
        point_q3[k] = percentile_accumulator_get(&p_acc, 75.0f);

      if (percentile_list != NULL && n_percentiles > 0 &&
          point_percentiles != NULL) {
        for (ptrdiff_t p = 0; p < n_percentiles; p++) {
          point_percentiles[k * n_percentiles + p] =
              percentile_accumulator_get(&p_acc, (float)percentile_list[p]);
        }
      }

      percentile_accumulator_free(&p_acc);
    }
  }

  // Pixel tracking: build offset array then do second pass to store indices
  if (track_pixels) {
    for (ptrdiff_t k = 0; k < n_track_points; k++) {
      point_offsets[k + 1] = point_offsets[k] + pixel_counts_temp[k];
      pixel_counts_temp[k] = 0; // Reset as insertion counters
    }

    // Second pass: store pixel indices
    for (ptrdiff_t k = 0; k < n_track_points; k++) {
      float center_dist = cum_dist[k];
      ptrdiff_t sub_start = k;
      ptrdiff_t sub_end = k;

      while (sub_start > 0 &&
             (center_dist - cum_dist[sub_start - 1]) <= binning_distance) {
        sub_start--;
      }
      while (sub_end < n_track_points - 1 &&
             (cum_dist[sub_end + 1] - center_dist) <= binning_distance) {
        sub_end++;
      }

      if (sub_start == sub_end) {
        if (sub_end < n_track_points - 1)
          sub_end++;
        else if (sub_start > 0)
          sub_start--;
      }

      ptrdiff_t n_sub_segments = sub_end - sub_start;
      if (n_sub_segments < 1)
        continue;

      // Same bounding box as first pass
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

      for (ptrdiff_t j = j_lo; j < j_hi; j++) {
        for (ptrdiff_t i = i_lo; i < i_hi; i++) {
          ptrdiff_t idx = j * dims[0] + i;
          if (isnan(dem[idx]))
            continue;

          float px = (float)i;
          float py = (float)j;

          float min_dist, proj_x, proj_y;
          ptrdiff_t seg = find_nearest_segment(
              px, py, track_i, track_j, n_track_points,
              sub_start, sub_end, &min_dist, &proj_x, &proj_y);

          if (seg < 0)
            continue;

          float dist_m = fabsf(min_dist) * cellsize;
          if (dist_m > half_width)
            continue;

          ptrdiff_t pos = point_offsets[k] + pixel_counts_temp[k];
          pixel_indices[pos] = idx;
          pixel_counts_temp[k]++;
        }
      }
    }
    free(pixel_counts_temp);
  }

  free(cum_dist);
}
