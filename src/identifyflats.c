#define TOPOTOOLBOX_BUILD

#include <stddef.h>
#include <stdint.h>

#include "topotoolbox.h"

/*
  Identify flat regions and sills in a digital elevation model.

  The arrays pointed to by `output` and `dem` should represent
  two-dimensional arrays of size (nrows, ncols). The output is a
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
ptrdiff_t identifyflats(int32_t *output, float *dem, ptrdiff_t nrows,
                        ptrdiff_t ncols) {
  ptrdiff_t col_offset[8] = {-1, -1, -1, 0, 0, 1, 1, 1};
  ptrdiff_t row_offset[8] = {-1, 0, 1, -1, 1, -1, 0, 1};

  ptrdiff_t count_flats = 0;

  // A flat is a pixel whose elevation is equal to the minimum
  // elevation of all of its neighbors.
  for (ptrdiff_t col = 0; col < ncols; col++) {
    for (ptrdiff_t row = 0; row < nrows; row++) {
      // Zero the output for all non-flat/sill pixels
      output[col * nrows + row] = 0;

      // Skip border pixels
      if (col == 0 || col == ncols - 1 || row == 0 || row == nrows - 1) {
        continue;
      }

      float dem_height = dem[col * nrows + row];
      float min_height = dem_height;

      // Compute the minimum height in the neighborhood around the
      // current pixel
      for (int32_t neighbor = 0; neighbor < 8; neighbor++) {
        ptrdiff_t neighbor_row = row + row_offset[neighbor];
        ptrdiff_t neighbor_col = col + col_offset[neighbor];

        // neighbor_row and neighbor_col are valid indices because we
        // skipped border pixels above
        float neighbor_height = dem[neighbor_col * nrows + neighbor_row];
        min_height =
            min_height < neighbor_height ? min_height : neighbor_height;
      }

      if (dem_height == min_height) {
        // Pixel is a flat
        output[col * nrows + row] |= 1;
        count_flats++;
      }
    }
  }

  // A sill is a pixel that
  // 1. is not a flat
  // 2. borders at least one flat
  // 3. has the same elevation as a flat that it touches
  for (ptrdiff_t col = 0; col < ncols; col++) {
    for (ptrdiff_t row = 0; row < nrows; row++) {
      if (output[col * nrows + row] & 1) {
        // Pixel is a flat, skip it
        continue;
      }

      float dem_height = dem[col * nrows + row];

      for (int32_t neighbor = 0; neighbor < 8; neighbor++) {
        ptrdiff_t neighbor_row = row + row_offset[neighbor];
        ptrdiff_t neighbor_col = col + col_offset[neighbor];

        if (neighbor_row < 0 || neighbor_row >= nrows || neighbor_col < 0 ||
            neighbor_col >= ncols) {
          // Skip neighbors outside the image
          continue;
        }

        if (!(output[neighbor_col * nrows + neighbor_row] & 1)) {
          // Neighbor is not a flat, skip it
          continue;
        }

        float neighbor_height = dem[neighbor_col * nrows + neighbor_row];
        if (neighbor_height == dem_height) {
          // flat neighbor has a height equal to that of the current pixel.
          // Current pixel is a sill
          output[col * nrows + row] |= 2;

          // Neighboring pixel is a presill
          output[neighbor_col * nrows + neighbor_row] |= 4;
          continue;
        }
      }
    }
  }
  return count_flats;
}
