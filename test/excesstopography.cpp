#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iostream>

extern "C" {
#include "topotoolbox.h"
#include "utils.h"
}

float upwind_gradient(float *u, ptrdiff_t row, ptrdiff_t col, float cellsize,
                      ptrdiff_t nrows, ptrdiff_t ncols) {
  float north_gradient =
      (u[col * nrows + row] - u[col * nrows + row - 1]) / cellsize;
  float south_gradient =
      (u[col * nrows + row + 1] - u[col * nrows + row]) / cellsize;

  float ns_gradient =
      std::fmaxf(0.0f, std::fmaxf(north_gradient, -south_gradient));

  float west_gradient =
      (u[col * nrows + row] - u[(col - 1) * nrows + row]) / cellsize;
  float east_gradient =
      (u[(col + 1) * nrows + row] - u[col * nrows + row]) / cellsize;

  float ew_gradient =
      std::fmaxf(0.0f, std::fmaxf(west_gradient, -east_gradient));
  float g = sqrtf(ns_gradient * ns_gradient + ew_gradient * ew_gradient);

  return g;
}

int32_t random_dem_test(ptrdiff_t nrows, ptrdiff_t ncols, ptrdiff_t nlayers,
                        uint32_t seed) {
  float *dem = new float[nrows * ncols];
  float *fmm_excess = new float[nrows * ncols];
  float *fsm_excess = new float[nrows * ncols];
  float *fmm_excess3d = new float[nrows * ncols];
  float *lithstack = new float[nlayers * nrows * ncols];
  float *threshold_slopes3d = new float[nlayers];
  ptrdiff_t *heap = new ptrdiff_t[nrows * ncols];
  ptrdiff_t *heap3d = new ptrdiff_t[nrows * ncols];
  ptrdiff_t *back = new ptrdiff_t[nrows * ncols];
  ptrdiff_t *back3d = new ptrdiff_t[nrows * ncols];
  float *threshold = new float[nrows * ncols];

  for (uint32_t col = 0; col < ncols; col++) {
    for (uint32_t row = 0; row < nrows; row++) {
      ptrdiff_t idx = col * nrows + row;
      dem[idx] = 100.0f * pcg4d(row, col, seed, 0);
      fsm_excess[idx] = dem[idx];
      fmm_excess[idx] = dem[idx];
      threshold[idx] = pcg4d(row, col, seed, 1);

      // Initialize lithstack with nlayers equally-spaced layers up to
      // 10 m greater than the maximum elevation.
      for (ptrdiff_t layer = 0; layer < nlayers; layer++) {
        lithstack[(col * nrows + row) * nlayers + layer] =
            110.0f * (layer + 1) / nlayers;
      }
    }
  }
  float cellsize = 30.0;

  for (ptrdiff_t layer = 0; layer < nlayers; layer++) {
    // Alternate hard and soft layers in the three-dimensional excess topography
    threshold_slopes3d[layer] = (layer % 2 == 0) ? 1.0f : 0.5f;
  }

  excesstopography_fsm2d(fsm_excess, dem, threshold, cellsize, nrows, ncols);
  excesstopography_fmm2d(fmm_excess, heap, back, dem, threshold, cellsize,
                         nrows, ncols);
  excesstopography_fmm3d(fmm_excess3d, heap3d, back3d, dem, lithstack,
                         threshold_slopes3d, cellsize, nrows, ncols, nlayers);
  for (ptrdiff_t col = 0; col < ncols; col++) {
    for (ptrdiff_t row = 0; row < nrows; row++) {
      ptrdiff_t idx = col * nrows + row;

      if (fmm_excess[idx] > dem[idx]) {
        std::cout << "Seed (" << seed << "): Pixel (" << row << ", " << col
                  << ") excess topography is greater than DEM" << std::endl;
        return -1;
      };

      if (fsm_excess[idx] > dem[idx]) {
        std::cout << "Seed (" << seed << "): Pixel (" << row << ", " << col
                  << ") FSM excess topography is greater than DEM" << std::endl;
        return -1;
      }

      if (fmm_excess[idx] > dem[idx]) {
        std::cout << "Seed (" << seed << "): Pixel (" << row << ", " << col
                  << ") 3D excess topograhy is greater than DEM" << std::endl;
      }

      // Test the upwind gradient
      if (col > 0 && col < ncols - 1 && row > 0 && row < nrows - 1) {
        float g = upwind_gradient(fmm_excess, row, col, cellsize, nrows, ncols);

        if ((g - threshold[col * nrows + row]) > 1e-4) {
          std::cout << "Seed (" << seed << "): Pixel (" << row << ", " << col
                    << ") discrete gradient " << g
                    << " is greater than threshold" << std::endl;
          return -1;
        };
      }

      // Test the FSM upwind gradient
      if (col > 0 && col < ncols - 1 && row > 0 && row < nrows - 1) {
        float g = upwind_gradient(fsm_excess, row, col, cellsize, nrows, ncols);

        if ((g - threshold[col * nrows + row]) > 1e-4) {
          std::cout << "Seed (" << seed << "): Pixel (" << row << ", " << col
                    << ") FSM discrete gradient " << g
                    << " is greater than threshold" << std::endl;
          return -1;
        };
      }
    }
  }

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
int32_t eikonal_numerics_test(ptrdiff_t nrows, ptrdiff_t ncols, uint32_t test) {
  float *dem = new float[nrows * ncols];
  float *excess = new float[nrows * ncols];
  ptrdiff_t *heap = new ptrdiff_t[nrows * ncols];
  ptrdiff_t *back = new ptrdiff_t[nrows * ncols];
  float *threshold = new float[nrows * ncols];

  for (uint32_t j = 0; j < ncols; j++) {
    for (uint32_t i = 0; i < nrows; i++) {
      threshold[i + j * nrows] = pcg4d(i, j, test, 0);
      if (i > 0 && i < nrows - 1 && j > 0 && j < ncols - 1) {
        dem[i + j * nrows] = INFINITY;
      } else {
        dem[i + j * nrows] = 0.0f;
      }
    }
  }

  excesstopography_fmm2d(excess, heap, back, dem, threshold, 30.0, nrows,
                         ncols);

  for (ptrdiff_t j = 0; j < ncols; j++) {
    for (ptrdiff_t i = 0; i < nrows; i++) {
      if (excess[i + j * ncols] == INFINITY) {
        std::cout << "Pixel (" << i << ", " << j << ") is infinite"
                  << std::endl;
        return -1;
      }
    }
  }

  delete[] dem;
  delete[] excess;
  delete[] heap;
  delete[] back;
  delete[] threshold;

  return 0;
}

int main(int argc, char *argv[]) {
  ptrdiff_t nrows = 200;
  ptrdiff_t ncols = 100;

  for (uint32_t test = 0; test < 100; test++) {
    int32_t result = random_dem_test(nrows, ncols, 8, test);
    if (result < 0) {
      return result;
    }

    result = eikonal_numerics_test(nrows, ncols, test);
    if (result < 0) {
      return result;
    }
  }

  return 0;
}
