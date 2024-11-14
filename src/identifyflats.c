#define TOPOTOOLBOX_BUILD

#include <math.h>
#include <stddef.h>
#include <stdint.h>

#include "topotoolbox.h"

/*
  Identify flat regions and sills in a digital elevation model.

  The arrays pointed to by `output` and `dem` should represent
  two-dimensional arrays of size (dims[0], dims[1]). The output is a
  bitfield where

  Bit 0 is set if the pixel is a flat
  Bit 1 is set if the pixel is a sill
  Bit 2 is set if the pixel is a presill (i.e. a flat neighboring a sill)

  To test if pixel p is a flat, use

  output[p] & 1

  To test if pixel p is a sill, use

  output[p] & 2

  To test if pixel p is a presill, use

  output[p] & 4
 */
TOPOTOOLBOX_API
ptrdiff_t identifyflats(int32_t *output, float *dem, ptrdiff_t dims[2]) {
  ptrdiff_t j_offset[8] = {-1, -1, -1, 0, 0, 1, 1, 1};
  ptrdiff_t i_offset[8] = {-1, 0, 1, -1, 1, -1, 0, 1};

  ptrdiff_t count_flats = 0;

  // A flat is a pixel whose elevation is equal to the minimum
  // elevation of all of its neighbors.
  for (ptrdiff_t j = 0; j < dims[1]; j++) {
    for (ptrdiff_t i = 0; i < dims[0]; i++) {
      // Zero the output for all non-flat/sill pixels
      output[j * dims[0] + i] = 0;

      // Skip border pixels
      if (j == 0 || j == dims[1] - 1 || i == 0 || i == dims[0] - 1) {
        continue;
      }

      float dem_height = dem[j * dims[0] + i];
      float min_height = dem_height;

      // Compute the minimum height in the neighborhood around the
      // current pixel
      for (int32_t neighbor = 0; neighbor < 8; neighbor++) {
        ptrdiff_t neighbor_i = i + i_offset[neighbor];
        ptrdiff_t neighbor_j = j + j_offset[neighbor];

        // neighbor_i and neighbor_j are valid indices because we
        // skipped border pixels above
        float neighbor_height = dem[neighbor_j * dims[0] + neighbor_i];

        neighbor_height = isnan(neighbor_height) ? -INFINITY : neighbor_height;

        min_height = fminf(min_height, neighbor_height);
      }

      // If dem_height is a NaN, this will automatically fail, and the
      // pixel will not be counted as a flat
      if (dem_height == min_height) {
        // Pixel is a flat
        output[j * dims[0] + i] |= 1;
        count_flats++;
      }
    }
  }

  // A sill is a pixel that
  // 1. is not a flat
  // 2. borders at least one flat
  // 3. has the same elevation as a flat that it touches
  for (ptrdiff_t j = 0; j < dims[1]; j++) {
    for (ptrdiff_t i = 0; i < dims[0]; i++) {
      if (output[j * dims[0] + i] & 1) {
        // Pixel is a flat, skip it
        continue;
      }

      float dem_height = dem[j * dims[0] + i];

      for (int32_t neighbor = 0; neighbor < 8; neighbor++) {
        ptrdiff_t neighbor_i = i + i_offset[neighbor];
        ptrdiff_t neighbor_j = j + j_offset[neighbor];

        if (neighbor_i < 0 || neighbor_i >= dims[0] || neighbor_j < 0 ||
            neighbor_j >= dims[1]) {
          // Skip neighbors outside the image
          continue;
        }

        if (!(output[neighbor_j * dims[0] + neighbor_i] & 1)) {
          // Neighbor is not a flat, skip it
          continue;
        }

        float neighbor_height = dem[neighbor_j * dims[0] + neighbor_i];
        if (neighbor_height == dem_height) {
          // flat neighbor has a height equal to that of the current pixel.
          // Current pixel is a sill
          output[j * dims[0] + i] |= 2;

          // Neighboring pixel is a presill
          output[neighbor_j * dims[0] + neighbor_i] |= 4;
          continue;
        }
      }
    }
  }
  return count_flats;
}
