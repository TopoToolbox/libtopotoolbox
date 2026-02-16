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

        // Convert pixel coordinates to metric space for distance computation
        float px = (float)i * cellsize;
        float py = (float)j * cellsize;
        float min_dist, proj_x, proj_y;

        ptrdiff_t new_seg = find_nearest_segment(
            px, py, track_i, track_j, n_track_points, search_start, search_end,
            &min_dist, &proj_x, &proj_y);

        // Update if we found a segment and it's different from current assignment
        if (new_seg >= 0 && new_seg != temp_nearest[idx]) {
          temp_nearest[idx] = new_seg;
          distance[idx] = min_dist * cellsize;
          changes++;
        } else if (new_seg >= 0) {
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
    const float *restrict dem, const float *restrict track_i,
    const float *restrict track_j, ptrdiff_t n_track_points, ptrdiff_t dims[2],
    float cellsize, float half_width, float bin_resolution, ptrdiff_t n_bins,
    int normalize) {
  if (n_track_points < 2 || n_bins <= 0)
    return;

  // Allocate one accumulator per distance bin
  swath_stats_accumulator *accumulators =
      (swath_stats_accumulator *)calloc(n_bins, sizeof(swath_stats_accumulator));

  for (ptrdiff_t b = 0; b < n_bins; b++) {
    accumulator_init(&accumulators[b]);
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
    }
  }

  // Finalize statistics for each bin
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

  free(distance);
  free(accumulators);
}

// ============================================================================
// Per-Point Swath Profile - Along-Track Variation
// ============================================================================

TOPOTOOLBOX_API
void swath_longitudinal(
    float *restrict point_means, float *restrict point_stddevs,
    float *restrict point_mins, float *restrict point_maxs,
    ptrdiff_t *restrict point_counts, ptrdiff_t *restrict pixel_indices,
    ptrdiff_t *restrict point_offsets, const float *restrict dem,
    const float *restrict track_i, const float *restrict track_j,
    ptrdiff_t n_track_points, ptrdiff_t dims[2], float cellsize,
    float half_width) {
  if (n_track_points < 2)
    return;

  // Allocate one accumulator per track point
  swath_stats_accumulator *accumulators = (swath_stats_accumulator *)calloc(
      n_track_points, sizeof(swath_stats_accumulator));

  for (ptrdiff_t k = 0; k < n_track_points; k++) {
    accumulator_init(&accumulators[k]);
  }

  // Pixel tracking requires two passes: count then store
  ptrdiff_t *pixel_counts_temp = NULL;
  int track_pixels = (pixel_indices != NULL && point_offsets != NULL);

  if (track_pixels) {
    pixel_counts_temp = (ptrdiff_t *)calloc(n_track_points, sizeof(ptrdiff_t));
    if (point_offsets != NULL) {
      point_offsets[0] = 0;
    }
  }

  // First pass: accumulate statistics and count pixels per point
  // Each pixel is assigned to its nearest track segment
  for (ptrdiff_t j = 0; j < dims[1]; j++) {
    for (ptrdiff_t i = 0; i < dims[0]; i++) {
      float px = (float)i * cellsize;
      float py = (float)j * cellsize;

      // Find which track segment this pixel belongs to
      float min_dist, proj_x, proj_y;
      ptrdiff_t nearest_seg =
          find_nearest_segment(px, py, track_i, track_j, n_track_points, 0,
                               n_track_points - 1, &min_dist, &proj_x, &proj_y);

      if (nearest_seg < 0)
        continue;

      float dist = min_dist * cellsize;
      if (fabsf(dist) > half_width)
        continue;

      ptrdiff_t idx = j * dims[0] + i;
      if (isnan(dem[idx]))
        continue;

      // Accumulate to the segment's statistics
      accumulator_add(&accumulators[nearest_seg], dem[idx]);

      if (track_pixels) {
        pixel_counts_temp[nearest_seg]++;
      }
    }
  }

  // Build offset array: convert counts to cumulative offsets
  // point_offsets[k] = start index in pixel_indices for point k
  // point_offsets[k+1] - point_offsets[k] = number of pixels for point k
  if (track_pixels) {
    for (ptrdiff_t k = 0; k < n_track_points; k++) {
      point_offsets[k + 1] = point_offsets[k] + pixel_counts_temp[k];
      pixel_counts_temp[k] = 0; // Reuse as insertion counters
    }
  }

  // Second pass: store linear pixel indices in flat array
  // Must repeat distance computation since we don't cache results
  if (track_pixels) {
    for (ptrdiff_t j = 0; j < dims[1]; j++) {
      for (ptrdiff_t i = 0; i < dims[0]; i++) {
        float px = (float)i * cellsize;
        float py = (float)j * cellsize;

        float min_dist, proj_x, proj_y;
        ptrdiff_t nearest_seg = find_nearest_segment(
            px, py, track_i, track_j, n_track_points, 0, n_track_points - 1,
            &min_dist, &proj_x, &proj_y);

        if (nearest_seg < 0)
          continue;

        float dist = min_dist * cellsize;
        if (fabsf(dist) > half_width)
          continue;

        ptrdiff_t idx = j * dims[0] + i;
        if (isnan(dem[idx]))
          continue;

        // Insert pixel index into flat array at correct position
        ptrdiff_t pos = point_offsets[nearest_seg] + pixel_counts_temp[nearest_seg];
        pixel_indices[pos] = idx;
        pixel_counts_temp[nearest_seg]++;
      }
    }
    free(pixel_counts_temp);
  }

  // Finalize statistics for each track point
  for (ptrdiff_t k = 0; k < n_track_points; k++) {
    point_means[k] = accumulator_mean(&accumulators[k]);
    point_stddevs[k] = accumulator_stddev(&accumulators[k]);
    point_mins[k] = accumulators[k].count > 0 ? accumulators[k].min_val : NAN;
    point_maxs[k] = accumulators[k].count > 0 ? accumulators[k].max_val : NAN;
    point_counts[k] = accumulators[k].count;
  }

  free(accumulators);
}
