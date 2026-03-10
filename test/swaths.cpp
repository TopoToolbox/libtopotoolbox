/**
   @file swaths.cpp
   @brief Property tests for the swath profiling functions:
     - swath_frontier_distance_map
     - swath_boundary_dijkstra
     - voronoi_ridge_to_centreline
     - swath_longitudinal  (binning_distance=0 and binning_distance=10)
     - swath_get_point_pixels

   All tests use the same geometry: a straight horizontal track at row 64
   on a 128x128 grid, running from column 10 to column 110.  The half-width
   is 20 pixels.  Because the track is axis-aligned, perpendicular distances
   and swath statistics are analytically exact, making it possible to check
   function outputs against known values rather than just sanity bounds.
 */

#undef NDEBUG
#include <cassert>
#include <cfloat>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <vector>

namespace tt {
extern "C" {
#include "topotoolbox.h"
}
}  // namespace tt

// Grid and track parameters shared across tests.
static const ptrdiff_t NROWS = 128;
static const ptrdiff_t NCOLS = 128;
static const ptrdiff_t TRACK_ROW = 64;
static const ptrdiff_t COL_START = 10;
static const ptrdiff_t COL_END = 110;
static const ptrdiff_t N_TRACK = COL_END - COL_START + 1;  // 101
static const float HW_PX = 20.0f;

// Build the straight horizontal track arrays.
static void make_track(std::vector<float> &ti, std::vector<float> &tj) {
  ti.resize(N_TRACK);
  tj.resize(N_TRACK);
  for (ptrdiff_t k = 0; k < N_TRACK; k++) {
    ti[k] = (float)TRACK_ROW;
    tj[k] = (float)(COL_START + k);
  }
}

// Collect the outer boundary of the full swath mask (pixels inside the mask
// that are adjacent to outside-mask or out-of-bounds via 4-connectivity),
// then filter by the sign of signed_dist:
//   sign_filter > 0  ->  keep pixels with signed_dist >= 0  (positive half)
//   sign_filter < 0  ->  keep pixels with signed_dist <= 0  (negative half)
// This mirrors the Python implementation that seeds the two Dijkstra waves
// for the Voronoi ridge from opposite outer edges of the swath.
static std::vector<ptrdiff_t> swath_boundary_by_sign(
    const std::vector<int8_t> &mask, const std::vector<float> &signed_dist,
    int sign_filter, ptrdiff_t dims[2]) {
  static const int DI4[4] = {-1, 1, 0, 0};
  static const int DJ4[4] = {0, 0, -1, 1};
  std::vector<ptrdiff_t> seeds;
  for (ptrdiff_t pj = 0; pj < dims[1]; pj++) {
    for (ptrdiff_t pi = 0; pi < dims[0]; pi++) {
      ptrdiff_t idx = pj * dims[0] + pi;
      if (!mask[idx]) continue;
      float sd = signed_dist[idx];
      if (sign_filter > 0 && sd < 0.0f) continue;
      if (sign_filter < 0 && sd > 0.0f) continue;
      for (int d = 0; d < 4; d++) {
        ptrdiff_t ni = pi + DI4[d], nj = pj + DJ4[d];
        if (ni < 0 || ni >= dims[0] || nj < 0 || nj >= dims[1] ||
            !mask[nj * dims[0] + ni]) {
          seeds.push_back(idx);
          break;
        }
      }
    }
  }
  return seeds;
}

// Test swath_frontier_distance_map and swath_boundary_dijkstra.
//
// For a horizontal track at TRACK_ROW, the perpendicular distance from any
// pixel (pi, pj) to the nearest track segment is exactly |pi - TRACK_ROW|
// for interior columns (no endpoint clamping).  We verify best_abs against
// this value within a tolerance of 1.5 px (to allow for D8 stepping and
// segment clamping at column boundaries).
//
// The swath mask is built from best_abs < HW_PX.  Seeds for
// swath_boundary_dijkstra are all outer-boundary pixels of the full swath
// (both positive and negative halves).  For interior pixels, the D8 distance
// to the nearest boundary pixel should be approximately
//   HW_PX - |pi - TRACK_ROW| - 1
// i.e. it decreases from ~HW_PX at the track to 0 at the swath edge.
// We check a loose lower bound (expected - 1) to tolerate D8 asymmetry.
void test_distance_maps() {
  ptrdiff_t dims[2] = {NROWS, NCOLS};
  ptrdiff_t total = NROWS * NCOLS;

  std::vector<float> ti, tj;
  make_track(ti, tj);

  std::vector<float> best_abs(total);
  std::vector<float> signed_dist(total);
  std::vector<ptrdiff_t> nearest_pt(total);

  tt::swath_frontier_distance_map(best_abs.data(), signed_dist.data(),
                                  nearest_pt.data(), ti.data(), tj.data(),
                                  N_TRACK, dims, HW_PX, nullptr);

  // For interior columns (30-90) the expected frontier distance is |pi-64|.
  // Pixels not reached by the Dijkstra (best_abs >= FLT_MAX) are outside the
  // swath and skipped.
  for (ptrdiff_t pj = 30; pj <= 90; pj++) {
    for (ptrdiff_t pi = 0; pi < NROWS; pi++) {
      ptrdiff_t idx = pj * NROWS + pi;
      float abs_d = best_abs[idx];
      if (abs_d >= FLT_MAX) continue;
      float expected = (float)std::abs(pi - TRACK_ROW);
      assert(std::abs(abs_d - expected) <= 1.5f);
    }
  }
  std::cout << "swath_frontier_distance_map: OK\n";

  // Build the swath mask (pixels whose frontier distance is < HW_PX).
  std::vector<int8_t> swath_mask(total);
  for (ptrdiff_t i = 0; i < total; i++)
    swath_mask[i] = (best_abs[i] < HW_PX) ? 1 : 0;

  // Merge positive and negative outer-boundary seeds into a single seed set.
  std::vector<ptrdiff_t> pos_seeds =
      swath_boundary_by_sign(swath_mask, signed_dist, +1, dims);
  std::vector<ptrdiff_t> neg_seeds =
      swath_boundary_by_sign(swath_mask, signed_dist, -1, dims);
  std::vector<ptrdiff_t> all_seeds;
  all_seeds.insert(all_seeds.end(), pos_seeds.begin(), pos_seeds.end());
  all_seeds.insert(all_seeds.end(), neg_seeds.begin(), neg_seeds.end());
  assert(!all_seeds.empty());

  std::vector<float> bdy_dist(total);
  tt::swath_boundary_dijkstra(bdy_dist.data(), swath_mask.data(),
                              all_seeds.data(), (ptrdiff_t)all_seeds.size(),
                              dims);

  // For interior columns (30-90) check that bdy_dist[pi,pj] is at least
  // (HW_PX - |pi-TRACK_ROW| - 1) - 1.0 (double tolerance for D8 routing).
  // Rows within 2 px of the swath boundary are excluded to avoid edge effects.
  for (ptrdiff_t pj = 30; pj <= 90; pj++) {
    for (ptrdiff_t pi = TRACK_ROW - (ptrdiff_t)HW_PX + 2;
         pi <= TRACK_ROW + (ptrdiff_t)HW_PX - 2; pi++) {
      ptrdiff_t idx = pj * NROWS + pi;
      if (!swath_mask[idx]) continue;
      float expected = HW_PX - (float)std::abs(pi - TRACK_ROW) - 1.0f;
      assert(bdy_dist[idx] >= expected - 1.0f);
    }
  }
  std::cout << "swath_boundary_dijkstra: OK\n";
}

// Test voronoi_ridge_to_centreline.
//
// The Voronoi ridge is found by running two boundary Dijkstra waves from
// opposite outer edges of the swath (positive and negative half separately)
// and marking pixels where the two wave-fronts meet.  For a symmetric
// horizontal track, the ridge should coincide with the track row.
//
// We verify:
//   - At least one ridge pixel is returned.
//   - Every ridge pixel's row coordinate is within ±2 px of TRACK_ROW.
//     (The ±2 tolerance allows for D8 asymmetry near the segment endpoints.)
void test_voronoi_centreline() {
  ptrdiff_t dims[2] = {NROWS, NCOLS};
  ptrdiff_t total = NROWS * NCOLS;

  std::vector<float> ti, tj;
  make_track(ti, tj);

  std::vector<float> best_abs(total);
  std::vector<float> signed_dist(total);
  std::vector<ptrdiff_t> nearest_pt(total);

  tt::swath_frontier_distance_map(best_abs.data(), signed_dist.data(),
                                  nearest_pt.data(), ti.data(), tj.data(),
                                  N_TRACK, dims, HW_PX, nullptr);

  std::vector<int8_t> swath_mask(total);
  for (ptrdiff_t i = 0; i < total; i++)
    swath_mask[i] = (best_abs[i] < HW_PX) ? 1 : 0;

  // Seeds are the outer boundary of the full swath, split by sign of
  // signed_dist.  Using only the outer boundary (not the boundary between the
  // two halves) ensures the two Dijkstra waves start from opposite outer edges
  // and meet in the middle, producing the correct centreline.
  std::vector<ptrdiff_t> pos_seeds =
      swath_boundary_by_sign(swath_mask, signed_dist, +1, dims);
  std::vector<ptrdiff_t> neg_seeds =
      swath_boundary_by_sign(swath_mask, signed_dist, -1, dims);
  assert(!pos_seeds.empty() && !neg_seeds.empty());

  std::vector<float> dist_pos(total), dist_neg(total);
  tt::swath_boundary_dijkstra(dist_pos.data(), swath_mask.data(),
                              pos_seeds.data(), (ptrdiff_t)pos_seeds.size(),
                              dims);
  tt::swath_boundary_dijkstra(dist_neg.data(), swath_mask.data(),
                              neg_seeds.data(), (ptrdiff_t)neg_seeds.size(),
                              dims);

  std::vector<float> cl_i(total), cl_j(total), cl_w(total);
  ptrdiff_t n_ridge = tt::voronoi_ridge_to_centreline(
      cl_i.data(), cl_j.data(), cl_w.data(), dist_pos.data(), dist_neg.data(),
      best_abs.data(), HW_PX, nearest_pt.data(), ti.data(), tj.data(), N_TRACK,
      dims, /*cellsize=*/1.0f);

  std::cout << "voronoi_ridge_to_centreline: " << n_ridge << " pixels\n";
  assert(n_ridge > 0);

  // For a symmetric swath, all ridge pixels should sit on or immediately
  // adjacent to the original track row.
  for (ptrdiff_t k = 0; k < n_ridge; k++)
    assert(std::abs(cl_i[k] - (float)TRACK_ROW) <= 2.0f);

  std::cout << "voronoi_ridge_to_centreline: OK\n";
}

// Test swath_longitudinal (two window sizes) and swath_get_point_pixels.
//
// DEM: z[pi, pj] = MAX_ELEV - |pi - TRACK_ROW|
// A tent (pyramid) shape that is constant in the column direction.  Every
// column of the swath has the same elevation profile, so swath statistics
// are independent of which column(s) end up in a window.
//
// --- binning_distance = 0 ---
// Each output point aggregates only pixels whose nearest_point equals that
// track index (no window).  For any interior block:
//   count  > 0                       (at least the track-row pixel itself)
//   max   <= MAX_ELEV                (cannot exceed the DEM peak)
//   min   >= MAX_ELEV - HW_PX - 1   (loose bound; see BD=10 for exact value)
//   min  <= mean <= max              (ordering sanity)
//
// --- binning_distance = 10 ---
// The sliding window spans ±10 track points.  For interior points (margin 21)
// the window always contains:
//   - The track-row pixel of every block  => max = MAX_ELEV exactly
//   - Pixels at perpendicular distance HW_PX (rows 44 and 84), because both
//     frontier_cost and swath_longitudinal use strict >, so distance = HW_PX
//     is included.  Elevation at those rows = MAX_ELEV - HW_PX exactly.
//     A pixel at distance HW_PX has zero lateral offset (sqrt(20²-20²)=0), so
//     it is always assigned to the column-aligned track point and always falls
//     within BD=10.  => min = MAX_ELEV - HW_PX exactly
//   - counts >= those of binning_distance=0 (window is strictly wider)
//
// --- swath_get_point_pixels (binning_distance = 0) ---
// Returns pixels assigned to a single track point by the Dijkstra.  For a
// horizontal track, tie-breaking at column boundaries makes the exact count
// unpredictable, so we only assert a non-empty result and valid row coords.
void test_longitudinal_and_get_pixels() {
  ptrdiff_t dims[2] = {NROWS, NCOLS};
  ptrdiff_t total = NROWS * NCOLS;

  const float MAX_ELEV = 100.0f;

  // Build pyramid DEM.
  std::vector<float> dem(total);
  for (ptrdiff_t pj = 0; pj < NCOLS; pj++)
    for (ptrdiff_t pi = 0; pi < NROWS; pi++)
      dem[pj * NROWS + pi] = MAX_ELEV - (float)std::abs(pi - TRACK_ROW);

  std::vector<float> ti, tj;
  make_track(ti, tj);

  std::vector<float> best_abs(total);
  std::vector<float> signed_dist(total);
  std::vector<ptrdiff_t> nearest_pt(total);

  tt::swath_frontier_distance_map(best_abs.data(), signed_dist.data(),
                                  nearest_pt.data(), ti.data(), tj.data(),
                                  N_TRACK, dims, HW_PX, nullptr);

  // Cumulative distance along the horizontal track with cellsize=1: cum_dist[k]
  // = k.
  std::vector<float> cum_dist(N_TRACK);
  for (ptrdiff_t k = 0; k < N_TRACK; k++) cum_dist[k] = (float)k;

  // --- binning_distance = 0 ---
  std::vector<float> means(N_TRACK), stddevs(N_TRACK);
  std::vector<float> mins(N_TRACK), maxs(N_TRACK);
  std::vector<ptrdiff_t> counts(N_TRACK);

  ptrdiff_t n_out = tt::swath_longitudinal(
      means.data(), stddevs.data(), mins.data(), maxs.data(), counts.data(),
      nullptr, nullptr, nullptr, nullptr, 0, nullptr, dem.data(), ti.data(),
      tj.data(), N_TRACK, signed_dist.data(), dims, /*cellsize=*/1.0f,
      /*half_width=*/HW_PX, /*binning_distance=*/0.0f, nearest_pt.data(),
      cum_dist.data(), /*skip=*/1, nullptr, nullptr);

  assert(n_out == N_TRACK);

  // Interior track points (5 away from each column endpoint).
  for (ptrdiff_t k = 5; k < N_TRACK - 5; k++) {
    assert(counts[k] > 0);
    assert(maxs[k] <= MAX_ELEV + 0.5f);
    assert(mins[k] >= MAX_ELEV - HW_PX - 1.0f);
    assert(means[k] >= mins[k] && means[k] <= maxs[k]);
  }
  std::cout << "swath_longitudinal binning_distance=0: OK\n";

  // --- binning_distance = 10 ---
  // Window spans [k-10, k+10].  For the pyramid DEM, max and min are
  // analytically determined regardless of which exact pixels end up in each
  // block (see function-level comment above).
  // Safe interior margin: BD(10) + max lateral nearest_pt displacement for
  // distance-HW_PX pixels (0) + column-end buffer (5) + a few extra = 21.
  const float BD = 10.0f;
  std::vector<float> means5(N_TRACK), stddevs5(N_TRACK);
  std::vector<float> mins5(N_TRACK), maxs5(N_TRACK);
  std::vector<ptrdiff_t> counts5(N_TRACK);

  ptrdiff_t n_out5 = tt::swath_longitudinal(
      means5.data(), stddevs5.data(), mins5.data(), maxs5.data(),
      counts5.data(), nullptr, nullptr, nullptr, nullptr, 0, nullptr,
      dem.data(), ti.data(), tj.data(), N_TRACK, signed_dist.data(), dims,
      /*cellsize=*/1.0f, /*half_width=*/HW_PX, /*binning_distance=*/BD,
      nearest_pt.data(), cum_dist.data(), /*skip=*/1, nullptr, nullptr);

  assert(n_out5 == N_TRACK);

  for (ptrdiff_t k = 21; k < N_TRACK - 21; k++) {
    // Wider window must aggregate at least as many pixels.
    assert(counts5[k] >= counts[k]);
    // Track-row pixel always in window: max = MAX_ELEV.
    assert(std::abs(maxs5[k] - MAX_ELEV) <= 0.5f);
    // Pixels at distance exactly HW_PX always in window: min = MAX_ELEV -
    // HW_PX.
    assert(std::abs(mins5[k] - (MAX_ELEV - HW_PX)) <= 0.5f);
    assert(means5[k] >= mins5[k] && means5[k] <= maxs5[k]);
  }
  std::cout << "swath_longitudinal binning_distance=" << BD << ": OK\n";

  // --- swath_get_point_pixels (binning_distance = 0) ---
  // Query the middle track point.  Because binning_distance=0 restricts
  // results to pixels whose nearest_point equals exactly mid, and because
  // D8 Dijkstra tie-breaking splits contested column boundaries arbitrarily,
  // the exact count is not analytically fixed.  We only check:
  //   - at least one pixel is returned
  //   - every returned pixel's row is within HW_PX of the track row
  ptrdiff_t mid = N_TRACK / 2;
  std::vector<ptrdiff_t> px_i(total), px_j(total);
  ptrdiff_t n_px = tt::swath_get_point_pixels(
      px_i.data(), px_j.data(), ti.data(), tj.data(), N_TRACK, mid,
      signed_dist.data(), dims, /*cellsize=*/1.0f, /*half_width=*/HW_PX,
      /*binning_distance=*/0.0f, nearest_pt.data(), cum_dist.data());

  std::cout << "swath_get_point_pixels: " << n_px << " pixels\n";
  assert(n_px > 0);
  for (ptrdiff_t k = 0; k < n_px; k++)
    assert(std::abs(px_i[k] - TRACK_ROW) <= (ptrdiff_t)HW_PX);
  std::cout << "swath_get_point_pixels: OK\n";
}

int main() {
  test_distance_maps();
  test_voronoi_centreline();
  test_longitudinal_and_get_pixels();
  std::cout << "All swath tests passed.\n";
  return 0;
}
