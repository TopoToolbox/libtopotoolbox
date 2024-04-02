#include <cstddef>
#include <cstdint>
#include <iostream>

extern "C" {
#include "topotoolbox.h"
}

/*
  PCG4D hash function


  Jarzynski, Mark and Olano, Marc. (2020). Hash functions for GPU
  rendering. Journal of Computer Graphics Techniques. Vol 9,
  No. 3. 21-38.
 */
float pcg4d(uint32_t a, uint32_t b, uint32_t c, uint32_t d) {
  uint32_t x = a * 1664525u + 1013904223u;
  uint32_t y = b * 1664525u + 1013904223u;
  uint32_t z = c * 1664525u + 1013904223u;
  uint32_t w = d * 1664525u + 1013904223u;

  x += y * w;
  y += z * x;
  z += x * y;
  w += y * z;

  x ^= x >> 16;
  y ^= y >> 16;
  z ^= z >> 16;
  w ^= w >> 16;

  x += y * w;
  y += z * x;
  z += x * y;
  w += y * z;

  return (float)(w >> 8) * 0x1.0p-24;
}

int32_t random_dem_test(ptrdiff_t nrows, ptrdiff_t ncols, uint64_t seed) {
  // Initialize a random DEM
  float *dem = new float[nrows * ncols];
  float *output = new float[nrows * ncols];

  for (ptrdiff_t col = 0; col < ncols; col++) {
    for (ptrdiff_t row = 0; row < nrows; row++) {
      dem[col * nrows + row] = 100.0 * pcg4d(row, col, seed, 1);
    }
  }

  fillsinks(output, dem, nrows, ncols);

  ptrdiff_t col_offset[8] = {-1, -1, -1, 0, 1, 1, 1, 0};
  ptrdiff_t row_offset[8] = {1, 0, -1, -1, -1, 0, 1, 1};

  for (ptrdiff_t col = 1; col < ncols - 1; col++) {
    for (ptrdiff_t row = 1; row < nrows - 1; row++) {
      // Each pixel of the filled raster should be >= the DEM
      float z = output[col * nrows + row];

      if (z < dem[col * nrows + row]) {
        std::cout << "Pixel (" << row << ", " << col << ") is below the DEM"
                  << std::endl;
        std::cout << "Value: " << z << std::endl;
        std::cout << "DEM: " << dem[col * nrows + row] << std::endl;
        return -1;
      }

      // No pixel of the filled raster should be surrounded by
      // neighbors that are higher than it
      int32_t count = 0;
      for (ptrdiff_t neighbor = 0; neighbor < 8; neighbor++) {
        ptrdiff_t neighbor_row = row + row_offset[neighbor];
        ptrdiff_t neighbor_col = col + col_offset[neighbor];

        if (z < output[neighbor_col * nrows + neighbor_row]) {
          count++;
        }
      }

      if (count == 8) {
        std::cout << "Pixel (" << row << ", " << col << ") is a sink"
                  << std::endl;
        return -1;
      }
    }
  }

  return 0;
}

int main(int argc, char *argv[]) {
  ptrdiff_t nrows = 100;
  ptrdiff_t ncols = 200;

  for (uint64_t test = 0; test < 100; test++) {
    int32_t result = random_dem_test(nrows, ncols, test);
    if (result < 0) {
      return result;
    }
  }
}
