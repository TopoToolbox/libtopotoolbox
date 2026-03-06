#define TOPOTOOLBOX_BUILD

#include <float.h>
#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "helpers/deque.h"
#include "helpers/dijkstra.h"
#include "helpers/polyline.h"
#include "helpers/stat_func.h"
#include "topotoolbox.h"

// 8-connected neighbor offsets (row, col) — shared by all D8 operations.
static const int k_di8[8] = {-1, -1, -1, 0, 0, 1, 1, 1};
static const int k_dj8[8] = {-1, 0, 1, -1, 1, -1, 0, 1};

// ============================================================================
// Swath profile — overview
// ============================================================================
//
// This file implements the core C routines for swath profiling. Two profile
// types are supported:
//
//   Transverse — global statistics binned by signed cross-track distance.
//   Longitudinal — windowed statistics computed along the track.
//
// The following functions are exposed through the public API:
//
//   swath_frontier_distance_map   Dijkstra expansion from track seed pixels,
//                                 producing absolute and signed pixel-unit
//                                 distance maps and a nearest-track-point map.
//
//   swath_boundary_dijkstra       Inward D8 Dijkstra from swath-edge pixels,
//                                 producing a distance-from-boundary map.
//
//   voronoi_ridge_to_centreline   Extracts the centreline of a shape from a
//                                 boundary distance map via Voronoi ridges.
//
//   swath_longitudinal            Computes a full longitudinal swath profile.
//                                 Each pixel is assigned to exactly one track
//                                 point (nearest-point rule).
//
//   swath_get_point_pixels        Returns all pixels assigned to a single
//                                 track point (or a window of track points).
//
// Memory allocation for output arrays is the caller's responsibility, except
// for the internal percentile scratch buffer in swath_longitudinal.
// Transverse swath statistics and windowed longitudinal logic are implemented
// in the higher-level layer(s).
//
// Polyline geometry (Bresenham rasterisation, simplification) is in
// helpers/polyline.c.

// ============================================================================
// Boundary Dijkstra — inward D8 wavefront from swath-edge pixels
// ============================================================================
TOPOTOOLBOX_API
void swath_boundary_dijkstra(float *dist_out, const int8_t *swath_mask,
                             const ptrdiff_t *seeds, ptrdiff_t n_seeds,
                             ptrdiff_t dims[2]) {
  GridDijkstra gd;
  grid_dijkstra_init(&gd, dist_out, dims);
  for (ptrdiff_t s = 0; s < n_seeds; s++)
    grid_dijkstra_seed(&gd, seeds[s], 0.0f, NULL, NULL);
  grid_dijkstra_run(&gd, swath_mask, dijkstra_cost_d8, NULL, NULL);
  grid_dijkstra_free(&gd);
}

// ============================================================================
// Frontier distance map — perpendicular distance from track via Dijkstra
// ============================================================================
//
// Context and callbacks for the swath-specific frontier cost function.
// The cost uses point_to_segment_distance against the segment adjacent to the
// parent's nearest track point (clamped projection prevents phantom rays).
// on_update writes the nearest point index and signed distance.

typedef struct {
  const float *track_i, *track_j;
  ptrdiff_t n_track_points;
  ptrdiff_t nrows;
  ptrdiff_t *nearest_point;
  float *signed_dist;
  float max_dist_px;
  // Scratch: set by frontier_cost, read by frontier_on_update.
  ptrdiff_t _last_point;
  float _last_signed;
} FrontierCtx;

static float frontier_cost(ptrdiff_t from, ptrdiff_t to, float from_dist,
                           int dir, void *ctx_ptr) {
  (void)from_dist;
  (void)dir;
  FrontierCtx *ctx = (FrontierCtx *)ctx_ptr;
  ptrdiff_t k = ctx->nearest_point[from];
  if (k < 0) return FLT_MAX;
  float to_i = (float)(to % ctx->nrows);
  float to_j = (float)(to / ctx->nrows);
  // Use the segment adjacent to track point k.
  // k < n-1: segment k→k+1; last point: segment k-1→k.
  // point_to_segment_distance has clamped projection (t∈[0,1]), which
  // prevents phantom rays at every bend and at both endpoints.
  ptrdiff_t ka = (k < ctx->n_track_points - 1) ? k : k - 1;
  ptrdiff_t kb = ka + 1;
  float proj_x, proj_y, lam;
  float signed_d = point_to_segment_distance(
      to_i, to_j, ctx->track_i[ka], ctx->track_j[ka], ctx->track_i[kb],
      ctx->track_j[kb], &proj_x, &proj_y, &lam);
  float abs_d = fabsf(signed_d);
  if (abs_d > ctx->max_dist_px) return FLT_MAX;
  ctx->_last_point = k;
  ctx->_last_signed = signed_d;
  return abs_d;
}

static void frontier_on_update(ptrdiff_t idx, float new_dist, void *ctx_ptr) {
  (void)new_dist;
  FrontierCtx *ctx = (FrontierCtx *)ctx_ptr;
  if (ctx->nearest_point) ctx->nearest_point[idx] = ctx->_last_point;
  if (ctx->signed_dist) ctx->signed_dist[idx] = ctx->_last_signed;
}

TOPOTOOLBOX_API
void swath_frontier_distance_map(float *restrict best_abs,
                                 float *restrict signed_dist,
                                 ptrdiff_t *restrict nearest_point,
                                 const float *restrict track_i,
                                 const float *restrict track_j,
                                 ptrdiff_t n_track_points, ptrdiff_t dims[2],
                                 float max_dist_px, const int8_t *mask) {
  ptrdiff_t total = dims[0] * dims[1];
  if (total <= 0) return;

  if (nearest_point)
    for (ptrdiff_t i = 0; i < total; i++) nearest_point[i] = -1;
  if (signed_dist)
    for (ptrdiff_t i = 0; i < total; i++) signed_dist[i] = 0.0f;

  FrontierCtx ctx = {.track_i = track_i,
                     .track_j = track_j,
                     .n_track_points = n_track_points,
                     .nrows = dims[0],
                     .nearest_point = nearest_point,
                     .signed_dist = signed_dist,
                     .max_dist_px = max_dist_px,
                     ._last_point = -1,
                     ._last_signed = 0.0f};

  GridDijkstra gd;
  grid_dijkstra_init(&gd, best_abs, dims);

  // Seed: one seed per track point at its rounded pixel position.
  // Track is assumed D8-contiguous (guaranteed by prepare_track in Python).
  for (ptrdiff_t k = 0; k < n_track_points; k++) {
    ptrdiff_t pi = (ptrdiff_t)(track_i[k] + 0.5f);
    ptrdiff_t pj = (ptrdiff_t)(track_j[k] + 0.5f);
    if (pi < 0 || pi >= dims[0] || pj < 0 || pj >= dims[1]) continue;
    ptrdiff_t idx = pj * dims[0] + pi;
    if (mask && !mask[idx]) continue;
    ctx._last_point = k;
    ctx._last_signed = 0.0f;
    grid_dijkstra_seed(&gd, idx, 0.0f, frontier_on_update, &ctx);
  }

  grid_dijkstra_run(&gd, mask, frontier_cost, frontier_on_update, &ctx);
  grid_dijkstra_free(&gd);
}

// ============================================================================
// Public API
// ============================================================================

// ============================================================================
// Voronoi ridge — centre-line pixel extraction from two boundary DT fields
// ============================================================================
//
// A pixel is on the ridge if dist_pos <= dist_neg (it is equidistant or
// closer to the positive boundary) AND at least one 8-neighbour has
// dist_pos > dist_neg (is on the negative side).  End-cap pixels whose
// nearest-segment projection falls outside the track endpoints are excluded.
//
// centre_line_i/j/width must be caller-allocated (size dims[0]*dims[1]).
// Returns the number of ridge pixels written.
TOPOTOOLBOX_API
ptrdiff_t voronoi_ridge_to_centreline(
    float *centre_line_i, float *centre_line_j, float *centre_width,
    const float *dist_pos, const float *dist_neg, const float *best_abs,
    float hw_px, const ptrdiff_t *nearest_point, const float *track_i,
    const float *track_j, ptrdiff_t n_track_points, ptrdiff_t dims[2],
    float cellsize) {
  ptrdiff_t n_centre = 0;
  ptrdiff_t total = dims[0] * dims[1];
  for (ptrdiff_t idx = 0; idx < total; idx++) {
    if (best_abs[idx] == FLT_MAX || best_abs[idx] > hw_px) continue;
    float dp = dist_pos[idx], dn = dist_neg[idx];
    if (dp == FLT_MAX || dn == FLT_MAX) continue;

    int my_sign = (dp < dn) ? 1 : (dp > dn) ? -1 : 0;
    if (my_sign < 0) continue;

    ptrdiff_t ci = idx % dims[0];
    ptrdiff_t cj = idx / dims[0];

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
        }
      }
    }
    if (!on_cl) continue;

    // End-cap clipping: exclude pixels outside the track endpoints.
    ptrdiff_t k = nearest_point[idx];
    if (k == 0 || k == n_track_points - 1) {
      float di = (float)ci - track_i[k];
      float dj = (float)cj - track_j[k];
      float tdi, tdj;
      if (k == 0) {
        tdi = track_i[1] - track_i[0];
        tdj = track_j[1] - track_j[0];
      } else {
        tdi = track_i[k] - track_i[k - 1];
        tdj = track_j[k] - track_j[k - 1];
      }
      float along = di * tdi + dj * tdj;
      if (k == 0 && along < 0.0f) continue;
      if (k == n_track_points - 1 && along > 0.0f) continue;
    }

    centre_line_i[n_centre] = (float)ci;
    centre_line_j[n_centre] = (float)cj;
    centre_width[n_centre] = (dp + dn) * cellsize;
    n_centre++;
  }
  return n_centre;
}

// ============================================================================
// Longitudinal swath — per-track-point statistics
// ============================================================================
//
// Case 2 (binning_distance > 0):
//   Uses caller-supplied nearest_point map (from swath_frontier_distance_map)
//   to assign pixels to track points, then applies a sliding window of width
//   ±binning_distance metres along the cumulative track distance, using
//   Lemire monotonic deques for O(1) amortised min/max per output point.

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
    const ptrdiff_t *restrict nearest_point, const float *restrict cum_dist,
    ptrdiff_t skip, float *result_track_i, float *result_track_j) {
  if (n_track_points < 2) return 0;
  ptrdiff_t effective_skip = (skip < 1) ? 1 : skip;

  int compute_percentiles =
      (point_medians != NULL || point_q1 != NULL || point_q3 != NULL ||
       (percentile_list != NULL && n_percentiles > 0 &&
        point_percentiles != NULL));

  // --- Step 1: cumulative along-track distance ---
  // Caller-supplied (size n_track_points). cum_dist[k] = arc-length in meters
  // from track point 0 to track point k.

  // Percentile accumulators: one raw-value buffer per output track point.
  // These are populated lazily during the sliding-window pass (Step 4).
  percentile_accumulator *p_accumulators = NULL;
  if (compute_percentiles) {
    p_accumulators = (percentile_accumulator *)calloc(
        n_track_points, sizeof(percentile_accumulator));
    for (ptrdiff_t k = 0; k < n_track_points; k++)
      percentile_accumulator_init(&p_accumulators[k], 64);
  }

  {
    // --- Step 2: assign each DEM pixel to its nearest track point (block) ---
    // nearest_point[idx] was precomputed by swath_frontier_distance_map.
    // Here we accumulate per-block stats (sum, sum², count, min, max) for
    // every pixel within half_width of the track.  The sliding window in
    // Step 3 will aggregate these blocks into a ±binning_distance window.
    ptrdiff_t total = dims[0] * dims[1];

    stats_accumulator *blk_acc =
        (stats_accumulator *)calloc(n_track_points, sizeof(stats_accumulator));
    for (ptrdiff_t k = 0; k < n_track_points; k++)
      accumulator_init(&blk_acc[k]);

    // For percentiles we also need to remember which pixel indices belong to
    // each block, so we can replay them into the per-output-point percentile
    // accumulator during the sliding-window pass.  We store them in a CSR
    // (compressed sparse row) layout: pt_pixels[pt_offset[k]..pt_offset[k+1]]
    // are the flat pixel indices assigned to block k.
    ptrdiff_t *pt_offset = NULL, *pt_pixels = NULL, *pt_fill_arr = NULL;
    if (compute_percentiles) {
      // First pass: count pixels per block to size the CSR arrays.
      ptrdiff_t *cnt = (ptrdiff_t *)calloc(n_track_points, sizeof(ptrdiff_t));
      ptrdiff_t n_assigned = 0;
      for (ptrdiff_t idx = 0; idx < total; idx++) {
        ptrdiff_t a = nearest_point[idx];
        if (a < 0) continue;
        float d = distance_from_track[idx];
        if (isnan(d) || fabsf(d) > half_width) continue;
        if (isnan(dem[idx])) continue;
        cnt[a]++;
        n_assigned++;
      }
      // Build CSR row-start offsets from counts.
      pt_offset = (ptrdiff_t *)malloc((n_track_points + 1) * sizeof(ptrdiff_t));
      pt_offset[0] = 0;
      for (ptrdiff_t k = 0; k < n_track_points; k++)
        pt_offset[k + 1] = pt_offset[k] + cnt[k];
      free(cnt);
      pt_pixels = (ptrdiff_t *)malloc(n_assigned * sizeof(ptrdiff_t));
      pt_fill_arr = (ptrdiff_t *)calloc(n_track_points, sizeof(ptrdiff_t));
      // Second pass: fill block stats and CSR pixel list together.
      for (ptrdiff_t idx = 0; idx < total; idx++) {
        ptrdiff_t a = nearest_point[idx];
        if (a < 0) continue;
        float d = distance_from_track[idx];
        if (isnan(d) || fabsf(d) > half_width) continue;
        if (isnan(dem[idx])) continue;
        accumulator_add(&blk_acc[a], dem[idx]);
        pt_pixels[pt_offset[a] + pt_fill_arr[a]] = idx;
        pt_fill_arr[a]++;
      }
    } else {
      // No percentiles: single pass, block stats only.
      for (ptrdiff_t idx = 0; idx < total; idx++) {
        ptrdiff_t a = nearest_point[idx];
        if (a < 0) continue;
        float d = distance_from_track[idx];
        if (isnan(d) || fabsf(d) > half_width) continue;
        if (isnan(dem[idx])) continue;
        accumulator_add(&blk_acc[a], dem[idx]);
      }
    }

    // --- Step 3: sliding window along the track ---
    // The window [lo..hi] is a contiguous range of block indices such that
    // both cum_dist[lo] and cum_dist[hi] are within binning_distance of
    // cum_dist[pt].  As pt advances, lo and hi move monotonically right,
    // giving O(n) total work.
    //
    // sum/sum² are maintained incrementally (add on right, subtract on left).
    // min/max use SlidingExtremaDQ (Lemire) for O(1) amortised per step.
    ptrdiff_t *min_idx =
        (ptrdiff_t *)malloc(n_track_points * sizeof(ptrdiff_t));
    ptrdiff_t *max_idx =
        (ptrdiff_t *)malloc(n_track_points * sizeof(ptrdiff_t));
    float *min_val = (float *)malloc(n_track_points * sizeof(float));
    float *max_val = (float *)malloc(n_track_points * sizeof(float));
    SlidingExtremaDQ min_dq, max_dq;
    sedq_init(&min_dq, min_idx, min_val);
    sedq_init(&max_dq, max_idx, max_val);
    double win_sum = 0.0, win_sum2 = 0.0;
    ptrdiff_t win_count = 0;
    ptrdiff_t lo = 0, hi = -1;

    for (ptrdiff_t pt = 0; pt < n_track_points; pt++) {
      // Expand right: absorb blocks whose distance from pt is ≤
      // binning_distance.
      while (hi < n_track_points - 1 &&
             cum_dist[hi + 1] - cum_dist[pt] <= binning_distance) {
        hi++;
        win_sum += blk_acc[hi].sum;
        win_sum2 += blk_acc[hi].sum_sq;
        win_count += blk_acc[hi].count;
        if (blk_acc[hi].count > 0) {
          sedq_push_min(&min_dq, hi, blk_acc[hi].min_val);
          sedq_push_max(&max_dq, hi, blk_acc[hi].max_val);
        }
      }
      // Shrink left: evict blocks that have fallen more than binning_distance
      // behind pt.
      while (lo <= hi && cum_dist[pt] - cum_dist[lo] > binning_distance) {
        win_sum -= blk_acc[lo].sum;
        win_sum2 -= blk_acc[lo].sum_sq;
        win_count -= blk_acc[lo].count;
        sedq_evict(&min_dq, lo);
        sedq_evict(&max_dq, lo);
        lo++;
      }
      // Write output for kept track points (honouring skip).
      if (effective_skip <= 1 || pt % effective_skip == 0) {
        ptrdiff_t out = pt / effective_skip;
        if (win_count > 0) {
          double mean = win_sum / (double)win_count;
          // Var = E[x²] - E[x]² (online formula, stable for float precision
          // here).
          double var = win_sum2 / (double)win_count - mean * mean;
          point_counts[out] = win_count;
          point_means[out] = (float)mean;
          point_stddevs[out] = (float)sqrt(var > 0.0 ? var : 0.0);
          point_mins[out] =
              !sedq_empty(&min_dq) ? sedq_front_val(&min_dq) : NAN;
          point_maxs[out] =
              !sedq_empty(&max_dq) ? sedq_front_val(&max_dq) : NAN;
        } else {
          point_counts[out] = 0;
          point_means[out] = NAN;
          point_stddevs[out] = NAN;
          point_mins[out] = NAN;
          point_maxs[out] = NAN;
        }
        if (result_track_i) result_track_i[out] = track_i[pt];
        if (result_track_j) result_track_j[out] = track_j[pt];
      }
      // For percentiles: replay all raw pixel values in the current window
      // [lo..hi] into this output point's percentile accumulator.
      // This is more expensive than mean/std (O(pixels in window) per point)
      // but unavoidable since percentiles don't decompose additively.
      if (compute_percentiles &&
          (effective_skip <= 1 || pt % effective_skip == 0)) {
        for (ptrdiff_t k = lo; k <= hi; k++)
          for (ptrdiff_t q = pt_offset[k]; q < pt_offset[k + 1]; q++)
            percentile_accumulator_add(&p_accumulators[pt], dem[pt_pixels[q]]);
      }
    }

    free(blk_acc);
    free(min_idx);
    free(max_idx);
    free(min_val);
    free(max_val);
    if (compute_percentiles) {
      free(pt_offset);
      free(pt_pixels);
      free(pt_fill_arr);
    }
  }

  // --- Step 4: percentile finalization ---
  // Sort each output point's raw-value buffer and query the requested ranks.
  if (compute_percentiles) {
    for (ptrdiff_t k = 0; k < n_track_points; k++) {
      if (effective_skip > 1 && k % effective_skip != 0) continue;
      ptrdiff_t out = k / effective_skip;
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
        for (ptrdiff_t p = 0; p < n_percentiles; p++)
          point_percentiles[out * n_percentiles + p] =
              percentile_accumulator_get(&p_accumulators[k],
                                         (float)percentile_list[p]);
      }
    }
    for (ptrdiff_t k = 0; k < n_track_points; k++)
      percentile_accumulator_free(&p_accumulators[k]);
    free(p_accumulators);
  }

  // Count written output points (= ceil(n_track_points / effective_skip)).
  ptrdiff_t n_result = 0;
  for (ptrdiff_t k = 0; k < n_track_points; k++) {
    if (effective_skip > 1 && k % effective_skip != 0) continue;
    n_result = k / effective_skip + 1;
  }

  return n_result;
}

// ============================================================================
// Per-point pixel retrieval — mirrors swath_longitudinal (vanilla)
// ============================================================================
//
// Case 1 (binning_distance <= 0): Bresenham cross-section filtered by
//   distance_from_track.
// Case 2 (binning_distance  > 0): uses caller-supplied nearest_point map;
//   gathers all pixels whose nearest track point falls within [pt_lo, pt_hi]
//   (the ±binning_distance window around point_index along cumulative
//   distance).

TOPOTOOLBOX_API
ptrdiff_t swath_get_point_pixels(
    ptrdiff_t *restrict pixels_i, ptrdiff_t *restrict pixels_j,
    const float *restrict track_i, const float *restrict track_j,
    ptrdiff_t n_track_points, ptrdiff_t point_index,
    const float *restrict distance_from_track, ptrdiff_t dims[2],
    float cellsize, float half_width, float binning_distance,
    const ptrdiff_t *restrict nearest_point, const float *restrict cum_dist) {
  if (n_track_points < 2) return 0;
  if (point_index < 0 || point_index >= n_track_points) return 0;

  ptrdiff_t n_pixels = 0;

  ptrdiff_t pt_lo = point_index, pt_hi = point_index;
  while (pt_lo > 0 &&
         cum_dist[point_index] - cum_dist[pt_lo - 1] <= binning_distance)
    pt_lo--;
  while (pt_hi < n_track_points - 1 &&
         cum_dist[pt_hi + 1] - cum_dist[point_index] <= binning_distance)
    pt_hi++;

  for (ptrdiff_t pj = 0; pj < dims[1]; pj++) {
    for (ptrdiff_t pi = 0; pi < dims[0]; pi++) {
      ptrdiff_t idx = pj * dims[0] + pi;
      ptrdiff_t a = nearest_point[idx];
      if (a < pt_lo || a > pt_hi) continue;
      float d = distance_from_track[idx];
      if (isnan(d) || fabsf(d) > half_width) continue;
      pixels_i[n_pixels] = pi;
      pixels_j[n_pixels] = pj;
      n_pixels++;
    }
  }

  return n_pixels;
}
