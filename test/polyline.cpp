#undef NDEBUG
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <vector>

namespace tt {
extern "C" {
#include "topotoolbox.h"
}
}  // namespace tt

// 128x128 grid.
// Aggressive saw-tooth pattern: segments with ~1:2 slope alternating left/right.
// Bresenham on a 1:2 slope produces alternating diagonal+cardinal steps, so
// roughly every other pixel is a removable elbow.
//
// Points (row, col) and per-segment D8 bound = max(|di|,|dj|):
//   (10, 5)->(20,25): max(10,20)=20
//   (20,25)->(35, 8): max(15,17)=17
//   (35, 8)->(45,30): max(10,22)=22
//   (45,30)->(60, 5): max(15,25)=25
//   (60, 5)->(70,28): max(10,23)=23
//   (70,28)->(85, 5): max(15,23)=23
//   (85, 5)->(95,30): max(10,25)=25
//   (95,30)->(110,8): max(15,22)=22
//   (110,8)->(118,118): max(8,110)=110
// Sum ~ 287. Loose check: > 200 and < 500.

static const ptrdiff_t N_REFS = 10;
static const ptrdiff_t ref_i[N_REFS] = {10, 20, 35, 45, 60, 70, 85, 95, 110, 118};
static const ptrdiff_t ref_j[N_REFS] = {5, 25, 8, 30, 5, 28, 5, 30, 8, 118};

void test_rasterize_path() {
  ptrdiff_t dims[2] = {128, 128};
  ptrdiff_t max_out = 4096;
  std::vector<ptrdiff_t> out_i(max_out), out_j(max_out);

  ptrdiff_t n =
      tt::rasterize_path(out_i.data(), out_j.data(), ref_i, ref_j, N_REFS,
                         /*close_loop=*/0, /*use_d4=*/0);

  std::cout << "rasterize_path D8: " << n << " pixels\n";
  assert(n > 200 && n < 500);

  // All output pixels must lie within the grid.
  for (ptrdiff_t k = 0; k < n; k++) {
    assert(out_i[k] >= 0 && out_i[k] < dims[0]);
    assert(out_j[k] >= 0 && out_j[k] < dims[1]);
  }

  // First and last pixel must match the first and last reference point.
  assert(out_i[0] == ref_i[0] && out_j[0] == ref_j[0]);
  assert(out_i[n - 1] == ref_i[N_REFS - 1] &&
         out_j[n - 1] == ref_j[N_REFS - 1]);
}

void test_thin_rasterised_line_to_D8() {
  ptrdiff_t dims[2] = {128, 128};

  // Manually built staircase: alternating cardinal-down + cardinal-right steps.
  // Pattern: (0,0),(1,0),(1,1),(2,1),(2,2),(3,2),...,(K,K)
  // Each corner pixel at odd indices has D8-adjacent neighbours on both sides,
  // so thinning removes it. With 20 steps -> 41 pixels, ~20 removed.
  const int K = 20;
  const int n = 2 * K + 1;
  std::vector<float> fi(n), fj(n), fw(n, 0.0f);
  for (int k = 0; k < K; k++) {
    fi[2 * k]     = (float)k;       fj[2 * k]     = (float)k;
    fi[2 * k + 1] = (float)(k + 1); fj[2 * k + 1] = (float)k;
  }
  fi[2 * K] = (float)K;
  fj[2 * K] = (float)K;

  ptrdiff_t n_thin =
      tt::thin_rasterised_line_to_D8(fi.data(), fj.data(), fw.data(), n, dims);

  std::cout << "thin_rasterised_line_to_D8: " << n << " -> " << n_thin
            << " pixels\n";
  assert(n_thin > 0);
  // Every other pixel is a removable elbow: expect at least 20% gone.
  assert(n_thin < n * 80 / 100);
}

void test_simplify_line() {
  ptrdiff_t dims[2] = {128, 128};
  ptrdiff_t max_out = 4096;
  std::vector<ptrdiff_t> out_i(max_out), out_j(max_out);

  ptrdiff_t n =
      tt::rasterize_path(out_i.data(), out_j.data(), ref_i, ref_j, N_REFS,
                         /*close_loop=*/0, /*use_d4=*/0);

  std::vector<float> fi(n), fj(n);
  for (ptrdiff_t k = 0; k < n; k++) {
    fi[k] = (float)out_i[k];
    fj[k] = (float)out_j[k];
  }

  std::vector<float> si(n), sj(n);

  // Method 0 (FIXED_N): tolerance = target point count.
  ptrdiff_t n10 = tt::simplify_line(si.data(), sj.data(), fi.data(), fj.data(),
                                    n, 10.0f, 0);
  ptrdiff_t n20 = tt::simplify_line(si.data(), sj.data(), fi.data(), fj.data(),
                                    n, 20.0f, 0);
  std::cout << "simplify_line method 0: n10=" << n10 << " n20=" << n20 << "\n";
  assert(n10 <= 10);
  assert(n20 <= 20);
  assert(n10 < n20);

  // Method 1 (KNEEDLE): tolerance ignored, automatic.
  ptrdiff_t n_knee = tt::simplify_line(si.data(), sj.data(), fi.data(),
                                       fj.data(), n, 0.0f, 1);
  std::cout << "simplify_line method 1 (kneedle): " << n_knee << "\n";
  assert(n_knee >= 2 && n_knee < n);

  // Method 2 (VW_AREA): tolerance = triangle area threshold.
  // Higher tolerance -> fewer points kept.
  ptrdiff_t n_vw_lo = tt::simplify_line(si.data(), sj.data(), fi.data(),
                                        fj.data(), n, 1.0f, 2);
  ptrdiff_t n_vw_hi = tt::simplify_line(si.data(), sj.data(), fi.data(),
                                        fj.data(), n, 50.0f, 2);
  std::cout << "simplify_line method 2 (VW): lo=" << n_vw_lo
            << " hi=" << n_vw_hi << "\n";
  assert(n_vw_lo >= 2);
  assert(n_vw_hi >= 2);
  assert(n_vw_hi <= n_vw_lo);
}

int main() {
  test_rasterize_path();
  test_thin_rasterised_line_to_D8();
  test_simplify_line();
  std::cout << "All polyline tests passed.\n";
  return 0;
}
