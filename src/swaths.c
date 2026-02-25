#define TOPOTOOLBOX_BUILD

#include <float.h>
#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "dijkstra.h"
#include "polyline.h"
#include "stat_helper.h"
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

  Distance maps — swath_frontier_distance_map (raw pixel-unit Dijkstra outputs).

  Polyline geometry (Bresenham, tangent, simplification) lives in polyline.c.
*/


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
// The cost function computes perpendicular distance to the nearest track
// segment in a ±2 window around the parent's nearest segment.
// on_update writes the nearest segment index and signed distance.

typedef struct {
  const float *track_i, *track_j;
  ptrdiff_t n_track_points;
  ptrdiff_t nrows;
  ptrdiff_t *nearest_seg;
  float *signed_dist;
  float max_dist_px;
  // Scratch: set by frontier_cost, read by frontier_on_update.
  ptrdiff_t _last_seg;
  float _last_signed;
} FrontierCtx;

static float frontier_cost(ptrdiff_t from, ptrdiff_t to, float from_dist,
                           int dir, void *ctx_ptr) {
  (void)from_dist;
  (void)dir;
  FrontierCtx *ctx = (FrontierCtx *)ctx_ptr;
  float px = (float)(to % ctx->nrows);
  float py = (float)(to / ctx->nrows);

  ptrdiff_t seg = ctx->nearest_seg ? ctx->nearest_seg[from] : 0;
  ptrdiff_t seg_lo = seg > 1 ? seg - 2 : 0;
  ptrdiff_t seg_hi =
      seg + 2 < ctx->n_track_points - 1 ? seg + 2 : ctx->n_track_points - 2;

  float best_ad = FLT_MAX, best_sd = 0.0f;
  ptrdiff_t best_seg = -1;
  for (ptrdiff_t sk = seg_lo; sk <= seg_hi; sk++) {
    float proj_x, proj_y, lam;
    float d = point_to_segment_distance(
        px, py, ctx->track_i[sk], ctx->track_j[sk], ctx->track_i[sk + 1],
        ctx->track_j[sk + 1], &proj_x, &proj_y, &lam);
    if (fabsf(d) < best_ad) {
      best_ad = fabsf(d);
      best_sd = d;
      best_seg = sk;
    }
  }
  if (best_ad > ctx->max_dist_px) return FLT_MAX;
  ctx->_last_seg = best_seg;
  ctx->_last_signed = best_sd;
  return best_ad;
}

static void frontier_on_update(ptrdiff_t idx, float new_dist, void *ctx_ptr) {
  (void)new_dist;
  FrontierCtx *ctx = (FrontierCtx *)ctx_ptr;
  if (ctx->nearest_seg) ctx->nearest_seg[idx] = ctx->_last_seg;
  if (ctx->signed_dist) ctx->signed_dist[idx] = ctx->_last_signed;
}

TOPOTOOLBOX_API
void swath_frontier_distance_map(
    float *restrict best_abs, float *restrict signed_dist,
    ptrdiff_t *restrict nearest_seg, const float *restrict track_i,
    const float *restrict track_j, ptrdiff_t n_track_points, ptrdiff_t dims[2],
    float max_dist_px, const float *dem, const int8_t *mask) {
  ptrdiff_t total = dims[0] * dims[1];

  int free_best = 0;
  if (best_abs == NULL) {
    best_abs = (float *)malloc(total * sizeof(float));
    free_best = 1;
  }
  if (nearest_seg)
    for (ptrdiff_t i = 0; i < total; i++) nearest_seg[i] = -1;
  if (signed_dist)
    for (ptrdiff_t i = 0; i < total; i++) signed_dist[i] = 0.0f;

  // Pre-compute pixel mask from dem/mask inputs.
  int8_t *pixel_mask = NULL;
  if (dem || mask) {
    pixel_mask = (int8_t *)malloc(total * sizeof(int8_t));
    for (ptrdiff_t i = 0; i < total; i++) {
      pixel_mask[i] = 1;
      if (dem && isnan(dem[i])) pixel_mask[i] = 0;
      if (mask && !mask[i]) pixel_mask[i] = 0;
    }
  }

  FrontierCtx ctx = {.track_i = track_i,
                     .track_j = track_j,
                     .n_track_points = n_track_points,
                     .nrows = dims[0],
                     .nearest_seg = nearest_seg,
                     .signed_dist = signed_dist,
                     .max_dist_px = max_dist_px,
                     ._last_seg = -1,
                     ._last_signed = 0.0f};

  GridDijkstra gd;
  grid_dijkstra_init(&gd, best_abs, dims);

  // Seed: rasterize track segments at ~0.5px spacing.
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
      if (pi < 0 || pi >= dims[0] || pj < 0 || pj >= dims[1]) continue;
      ptrdiff_t idx = pj * dims[0] + pi;
      if (pixel_mask && !pixel_mask[idx]) continue;

      float px = (float)pi, py = (float)pj;
      float best_ad = FLT_MAX, best_sd = 0.0f;
      ptrdiff_t best_seg = -1;
      for (ptrdiff_t sk = seg_lo; sk <= seg_hi; sk++) {
        float proj_x, proj_y, lam;
        float d = point_to_segment_distance(px, py, track_i[sk], track_j[sk],
                                            track_i[sk + 1], track_j[sk + 1],
                                            &proj_x, &proj_y, &lam);
        if (fabsf(d) < best_ad) {
          best_ad = fabsf(d);
          best_sd = d;
          best_seg = sk;
        }
      }
      if (best_ad <= max_dist_px) {
        ctx._last_seg = best_seg;
        ctx._last_signed = best_sd;
        grid_dijkstra_seed(&gd, idx, best_ad, frontier_on_update, &ctx);
      }
    }
  }

  grid_dijkstra_run(&gd, pixel_mask, frontier_cost, frontier_on_update, &ctx);
  grid_dijkstra_free(&gd);
  free(pixel_mask);
  if (free_best) free(best_abs);
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
    float hw_px, const ptrdiff_t *nearest_seg, const float *track_i,
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

    // End-cap clipping: recompute unclamped projection t
    ptrdiff_t seg = nearest_seg[idx];
    float si = track_i[seg + 1] - track_i[seg];
    float sj = track_j[seg + 1] - track_j[seg];
    float len2 = si * si + sj * sj;
    if (len2 > 0.0f) {
      float t = ((ci - track_i[seg]) * si + (cj - track_j[seg]) * sj) / len2;
      if (seg == 0 && t < 0.0f) continue;
      if (seg == n_track_points - 2 && t > 1.0f) continue;
    }

    centre_line_i[n_centre] = (float)ci;
    centre_line_j[n_centre] = (float)cj;
    centre_width[n_centre] = (dp + dn) * cellsize;
    n_centre++;
  }
  return n_centre;
}


// ============================================================================
// Transverse swath — DEM pixels binned by signed perpendicular distance
// ============================================================================
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
  stats_accumulator *accumulators =
      (stats_accumulator *)calloc(n_bins, sizeof(stats_accumulator));
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
    stats_accumulator track_acc;
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
  stats_accumulator *accumulators = NULL;
  if (binning_distance <= 0.0f) {
    accumulators =
        (stats_accumulator *)calloc(n_track_points, sizeof(stats_accumulator));
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
