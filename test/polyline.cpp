/**
   @file polyline.cpp
   @brief Property tests for rasterize_path, thin_rasterised_line_to_D8,
   and simplify_line.

   All three tests reuse the same 10-point zigzag reference path on a 128x128
   grid.  The path has an aggressive 1:2 slope alternating left/right, which
   exercises diagonal Bresenham steps and produces many removable elbows.
 */

#undef NDEBUG
#include <cassert>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <vector>

namespace tt {
extern "C" {
#include "topotoolbox.h"
}
}  // namespace tt

// Reference control points for a 10-segment zigzag path on a 128x128 grid.
// Segments alternate between going right+down and right+up with roughly 1:2
// row-to-column slope, so Bresenham produces a dense diagonal raster.
// The last segment (110,8)->(118,118) is nearly horizontal and exercises the
// long-run case.  Expected D8 pixel count per segment:
//   (10,5)->(20,25) : max(10,20)=20     (35,8)->(45,30): max(10,22)=22
//   (20,25)->(35,8) : max(15,17)=17     (45,30)->(60,5): max(15,25)=25
//   (60,5)->(70,28) : max(10,23)=23     (70,28)->(85,5): max(15,23)=23
//   (85,5)->(95,30) : max(10,25)=25     (95,30)->(110,8):max(15,22)=22
//   (110,8)->(118,118): max(8,110)=110
// Total ≈ 287, checked loosely as > 200 and < 500.
static const ptrdiff_t N_REFS = 10;
static const ptrdiff_t ref_i[N_REFS] = {10, 20, 35, 45,  60,
                                        70, 85, 95, 110, 118};
static const ptrdiff_t ref_j[N_REFS] = {5, 25, 8, 30, 5, 28, 5, 30, 8, 118};

// Test rasterize_path.
//
// Checks:
//   - Output pixel count is in the expected range (>200, <500) based on the
//     per-segment D8 bound estimates above.
//   - Every output pixel lies within the 128x128 grid bounds.
//   - First output pixel matches the first reference point.
//   - Last output pixel matches the last reference point.
void test_rasterize_path() {
  std::vector<ptrdiff_t> out_i(4096), out_j(4096);

  ptrdiff_t n =
      tt::rasterize_path(out_i.data(), out_j.data(), ref_i, ref_j, N_REFS,
                         /*close_loop=*/0, /*use_d4=*/0);

  std::cout << "rasterize_path D8: " << n << " pixels\n";
  // Expected total D8 pixel count based on per-segment max(|di|,|dj|) sums.
  assert(n > 200 && n < 500);

  // All pixels must be inside the 128x128 grid.
  for (ptrdiff_t k = 0; k < n; k++) {
    assert(out_i[k] >= 0 && out_i[k] < 128);
    assert(out_j[k] >= 0 && out_j[k] < 128);
  }

  // Endpoints must be exactly the first and last reference points.
  assert(out_i[0] == ref_i[0] && out_j[0] == ref_j[0]);
  assert(out_i[n - 1] == ref_i[N_REFS - 1] &&
         out_j[n - 1] == ref_j[N_REFS - 1]);
}

// Test thin_rasterised_line_to_D8.
//
// A manually constructed staircase is used instead of the zigzag path because
// the thinning geometry needs to be exact.  The staircase alternates a
// cardinal-down step with a cardinal-right step:
//   (0,0) -> (1,0) -> (1,1) -> (2,1) -> (2,2) -> ... -> (K,K)
// Each "elbow" pixel at position (k+1, k) has D8-adjacent neighbours on both
// sides (the pixel before going diagonally, the pixel after going right), so
// thinning can remove it.  With K=20 steps the path has 41 pixels and about
// 20 elbows, giving a theoretical removal of ~49%.
//
// Checks:
//   - The thinned count is positive (path is not entirely removed).
//   - At least 20% of the original pixels were removed (n_thin < 0.8 * n).
void test_thin_rasterised_line_to_D8() {
  ptrdiff_t dims[2] = {128, 128};

  const ptrdiff_t K = 20;
  const ptrdiff_t n = 2 * K + 1;
  std::vector<float> fi(n), fj(n), fw(n, 0.0f);
  for (ptrdiff_t k = 0; k < K; k++) {
    fi[2 * k] = (float)k;
    fj[2 * k] = (float)k;
    fi[2 * k + 1] = (float)(k + 1);
    fj[2 * k + 1] = (float)k;
  }
  fi[2 * K] = (float)K;
  fj[2 * K] = (float)K;

  ptrdiff_t n_thin =
      tt::thin_rasterised_line_to_D8(fi.data(), fj.data(), fw.data(), n, dims);

  std::cout << "thin_rasterised_line_to_D8: " << n << " -> " << n_thin
            << " pixels\n";
  assert(n_thin > 0);
  // At least 20% of the 41-pixel staircase should be removed as redundant
  // elbows (expected ~49% removal; the 20% threshold is conservative).
  assert(n_thin < n * 80 / 100);
}

// Test simplify_line (all three methods) on the rasterized zigzag path.
//
// Method 0 (FIXED_N / IEF with target count):
//   tolerance = target number of output points.
//   Checks: output count <= tolerance, and a smaller target gives fewer points.
//
// Method 1 (KNEEDLE / IEF with automatic knee detection):
//   tolerance is ignored; the method finds the knee of the IEF curve itself.
//   Checks: output has at least 2 points and fewer than the input.
//
// Method 2 (VW_AREA / Visvalingam-Whyatt area threshold):
//   tolerance = minimum triangle area to keep a point.
//   Checks: higher tolerance produces fewer or equal output points (monotone).
void test_simplify_line() {
  std::vector<ptrdiff_t> out_i(4096), out_j(4096);

  ptrdiff_t n =
      tt::rasterize_path(out_i.data(), out_j.data(), ref_i, ref_j, N_REFS,
                         /*close_loop=*/0, /*use_d4=*/0);

  std::vector<float> fi(n), fj(n), si(n), sj(n);
  for (ptrdiff_t k = 0; k < n; k++) {
    fi[k] = (float)out_i[k];
    fj[k] = (float)out_j[k];
  }

  // Method 0: tolerance is the target output count.
  // n10 requests ≤10 points, n20 requests ≤20 points.
  ptrdiff_t n10 = tt::simplify_line(si.data(), sj.data(), fi.data(), fj.data(),
                                    n, 10.0f, 0);
  ptrdiff_t n20 = tt::simplify_line(si.data(), sj.data(), fi.data(), fj.data(),
                                    n, 20.0f, 0);
  std::cout << "simplify_line method 0: n10=" << n10 << " n20=" << n20 << "\n";
  assert(n10 <= 10);
  assert(n20 <= 20);
  // Larger target must produce at least as many points as smaller target.
  assert(n10 < n20);

  // Method 1: automatic knee detection; result is between 2 and n-1.
  ptrdiff_t n_knee =
      tt::simplify_line(si.data(), sj.data(), fi.data(), fj.data(), n, 0.0f, 1);
  std::cout << "simplify_line method 1 (kneedle): " << n_knee << "\n";
  assert(n_knee >= 2 && n_knee < n);

  // Method 2: Visvalingam-Whyatt.  Higher area threshold removes more points.
  ptrdiff_t n_vw_lo =
      tt::simplify_line(si.data(), sj.data(), fi.data(), fj.data(), n, 1.0f, 2);
  ptrdiff_t n_vw_hi = tt::simplify_line(si.data(), sj.data(), fi.data(),
                                        fj.data(), n, 50.0f, 2);
  std::cout << "simplify_line method 2 (VW): lo=" << n_vw_lo
            << " hi=" << n_vw_hi << "\n";
  assert(n_vw_lo >= 2);
  assert(n_vw_hi >= 2);
  // Higher tolerance should produce fewer or equal output points.
  assert(n_vw_hi <= n_vw_lo);
}

int main() {
  test_rasterize_path();
  test_thin_rasterised_line_to_D8();
  test_simplify_line();
  std::cout << "All polyline tests passed.\n";
  return 0;
}
