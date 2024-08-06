#undef NDEBUG
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iostream>

extern "C" {
#include "topotoolbox.h"
#include "utils.h"
}

float upwind_gradient(float *u, ptrdiff_t row, ptrdiff_t col, float cellsize,
                      ptrdiff_t dims[2]) {
  float north_gradient =
      (u[col * dims[0] + row] - u[col * dims[0] + row - 1]) / cellsize;
  float south_gradient =
      (u[col * dims[0] + row + 1] - u[col * dims[0] + row]) / cellsize;

  float ns_gradient =
      std::fmaxf(0.0f, std::fmaxf(north_gradient, -south_gradient));

  float west_gradient =
      (u[col * dims[0] + row] - u[(col - 1) * dims[0] + row]) / cellsize;
  float east_gradient =
      (u[(col + 1) * dims[0] + row] - u[col * dims[0] + row]) / cellsize;

  float ew_gradient =
      std::fmaxf(0.0f, std::fmaxf(west_gradient, -east_gradient));
  float g = sqrtf(ns_gradient * ns_gradient + ew_gradient * ew_gradient);

  return g;
}

int32_t test_excess_constraint(float *excess, float *dem, ptrdiff_t dims[2]) {
  for (ptrdiff_t j = 0; j < dims[1]; j++) {
    for (ptrdiff_t i = 0; i < dims[0]; i++) {
      assert(excess[i + j * dims[0]] <= dem[i + j * dims[0]]);
    }
  }
  return 0;
}

int32_t test_upwind_gradient(float *excess, float *threshold, float cellsize,
                             ptrdiff_t dims[2]) {
  for (ptrdiff_t j = 0; j < dims[1]; j++) {
    for (ptrdiff_t i = 0; i < dims[0]; i++) {
      if ((i > 0) && (i < dims[0] - 1) && (j > 0) && (j < dims[1] - 1)) {
        assert(upwind_gradient(excess, i, j, cellsize, dims) -
                   threshold[j * dims[0] + i] <
               1e-4);
      }
    }
  }
  return 0;
}

int32_t test_excess_isfinite(float *excess, ptrdiff_t dims[2]) {
  for (ptrdiff_t j = 0; j < dims[1]; j++) {
    for (ptrdiff_t i = 0; i < dims[0]; i++) {
      assert(excess[i + j * dims[0]] != INFINITY);
    }
  }
  return 0;
}

int32_t random_dem_test(ptrdiff_t dims[2], ptrdiff_t nlayers, uint32_t seed) {
  float *dem = new float[dims[0] * dims[1]];
  float *fmm_excess = new float[dims[0] * dims[1]];
  float *fsm_excess = new float[dims[0] * dims[1]];
  float *fmm_excess3d = new float[dims[0] * dims[1]];
  float *lithstack = new float[nlayers * dims[0] * dims[1]];
  float *threshold_slopes3d = new float[nlayers];
  ptrdiff_t *heap = new ptrdiff_t[dims[0] * dims[1]];
  ptrdiff_t *heap3d = new ptrdiff_t[dims[0] * dims[1]];
  ptrdiff_t *back = new ptrdiff_t[dims[0] * dims[1]];
  ptrdiff_t *back3d = new ptrdiff_t[dims[0] * dims[1]];
  float *threshold = new float[dims[0] * dims[1]];

  for (uint32_t col = 0; col < dims[1]; col++) {
    for (uint32_t row = 0; row < dims[0]; row++) {
      ptrdiff_t idx = col * dims[0] + row;
      dem[idx] = 100.0f * pcg4d(row, col, seed, 0);
      fsm_excess[idx] = dem[idx];
      fmm_excess[idx] = dem[idx];
      threshold[idx] = pcg4d(row, col, seed, 1);

      // Initialize lithstack with nlayers equally-spaced layers up to
      // 10 m greater than the maximum elevation.
      for (ptrdiff_t layer = 0; layer < nlayers; layer++) {
        lithstack[(col * dims[0] + row) * nlayers + layer] =
            110.0f * (layer + 1) / nlayers;
      }
    }
  }
  float cellsize = 30.0;

  for (ptrdiff_t layer = 0; layer < nlayers; layer++) {
    // Alternate hard and soft layers in the three-dimensional excess topography
    threshold_slopes3d[layer] = (layer % 2 == 0) ? 1.0f : 0.5f;
  }

  excesstopography_fsm2d(fsm_excess, dem, threshold, cellsize, dims);

  test_excess_constraint(fsm_excess, dem, dims);
  test_upwind_gradient(fsm_excess, threshold, cellsize, dims);

  excesstopography_fmm2d(fmm_excess, heap, back, dem, threshold, cellsize,
                         dims);

  test_excess_constraint(fmm_excess, dem, dims);
  test_upwind_gradient(fmm_excess, threshold, cellsize, dims);

  excesstopography_fmm3d(fmm_excess3d, heap3d, back3d, dem, lithstack,
                         threshold_slopes3d, cellsize, dims, nlayers);

  test_excess_constraint(fmm_excess3d, dem, dims);

  delete[] dem;
  delete[] fmm_excess;
  delete[] fsm_excess;
  delete[] fmm_excess3d;
  delete[] lithstack;
  delete[] threshold_slopes3d;
  delete[] heap;
  delete[] heap3d;
  delete[] back;
  delete[] back3d;
  delete[] threshold;

  return 0;
}

// This test represents an excesstopography use case that produces
// rounding errors that affect the eikonal solver. It uses random
// threshold slopes and fixes the boundary pixels to zero, but the
// interior of the DEM is made infinite and thus unconstrained. When
// slopes are very low and elevations relatively high, this can lead
// to taking the square root of negative numbers, so that the proposed
// elevation is NaN. This proposal is never accepted, and the pixels
// are left at their infinite values.
int32_t eikonal_numerics_test(ptrdiff_t dims[2], uint32_t test) {
  float *dem = new float[dims[0] * dims[1]];
  float *excess = new float[dims[0] * dims[1]];
  ptrdiff_t *heap = new ptrdiff_t[dims[0] * dims[1]];
  ptrdiff_t *back = new ptrdiff_t[dims[0] * dims[1]];
  float *threshold = new float[dims[0] * dims[1]];

  for (uint32_t j = 0; j < dims[1]; j++) {
    for (uint32_t i = 0; i < dims[0]; i++) {
      threshold[i + j * dims[0]] = pcg4d(i, j, test, 0);
      if (i > 0 && i < dims[0] - 1 && j > 0 && j < dims[1] - 1) {
        dem[i + j * dims[0]] = INFINITY;
      } else {
        dem[i + j * dims[0]] = 0.0f;
      }
    }
  }

  excesstopography_fmm2d(excess, heap, back, dem, threshold, 30.0, dims);

  test_upwind_gradient(excess, threshold, 30.0, dims);
  test_excess_constraint(excess, dem, dims);
  test_excess_isfinite(excess, dims);

  delete[] dem;
  delete[] excess;
  delete[] heap;
  delete[] back;
  delete[] threshold;

  return 0;
}

int main(int argc, char *argv[]) {
  ptrdiff_t dims[2] = {200, 100};

  for (uint32_t test = 0; test < 100; test++) {
    int32_t result = random_dem_test(dims, 8, test);
    if (result < 0) {
      return result;
    }

    result = eikonal_numerics_test(dims, test);
    if (result < 0) {
      return result;
    }
  }

  return 0;
}
