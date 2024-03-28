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
double pcg4d(uint64_t a, uint64_t b, uint64_t c, uint64_t d) {
  uint64_t x = a * 0x5851f42d4c957f2d + 0x14057b7ef767814f;
  uint64_t y = b * 0x5851f42d4c957f2d + 0x14057b7ef767814f;
  uint64_t z = c * 0x5851f42d4c957f2d + 0x14057b7ef767814f;
  uint64_t w = d * 0x5851f42d4c957f2d + 0x14057b7ef767814f;

  x += y * w;
  y += z * x;
  z += x * y;
  w += y * z;

  x ^= x >> 32;
  y ^= y >> 32;
  z ^= z >> 32;
  w ^= w >> 32;

  x += y * w;
  y += z * x;
  z += x * y;
  w += y * z;

  return (double)(w >> 11) * 0x1.0p-53;
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
      double z = output[col * nrows + row];

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
