#include "reconstruct.h"

#include <stddef.h>
#include <stdint.h>

/*
  Perform a partial reconstruction by scanning in the forward
  direction.

  Returns the number of pixels that were modified in the current scan.

  The forward scan replaces every pixel in `marker` with the maximum
  of a neighborhood consisting of the 4 neighbors already visited by
  the raster scan, denoted by 'x' in the diagram:

  x x .
   \|
  x-o .
   /
  x . .


  The new marker pixel value is constrained to lie below the
  corresponding pixel in `mask`.
 */
ptrdiff_t forward_scan(float *marker, float *mask, ptrdiff_t dims[2]) {
  // Offsets for the four neighbors
  ptrdiff_t j_offset[4] = {-1, -1, -1, 0};
  ptrdiff_t i_offset[4] = {1, 0, -1, -1};

  ptrdiff_t count = 0;  // Number of modified pixels

  for (ptrdiff_t j = 0; j < dims[1]; j++) {
    for (ptrdiff_t i = 0; i < dims[0]; i++) {
      ptrdiff_t p = j * dims[0] + i;

      // Compute the maximum of the marker at the current pixel and all
      // of its previously visisted neighbors
      float max_height = marker[p];
      for (ptrdiff_t neighbor = 0; neighbor < 4; neighbor++) {
        ptrdiff_t neighbor_i = i + i_offset[neighbor];
        ptrdiff_t neighbor_j = j + j_offset[neighbor];

        ptrdiff_t q = neighbor_j * dims[0] + neighbor_i;

        // Skip pixels outside the boundary
        if (neighbor_i < 0 || neighbor_i >= dims[0] || neighbor_j < 0 ||
            neighbor_j >= dims[1]) {
          continue;
        }

        max_height = max_height > marker[q] ? max_height : marker[q];
      }

      // Set the marker at the current pixel to the minimum of the
      // maximum height of the neighborhood and the mask at the current
      // pixel.

      float z = max_height < mask[p] ? max_height : mask[p];

      if (z != marker[p]) {
        // Increment count only if we change the current pixel
        count++;
        marker[p] = z;
      }
    }
  }
  return count;
}

/*
  Perform a partial reconstruction by scanning in the backward
  direction.

  Returns the number of pixels that were modified in the current scan.

  The backward scan replaces every pixel in `marker` with the maximum
  of a neighborhood consisting of the 4 neighbors already visited by
  the raster scan, denoted by 'x' in the diagram:

  . . x
     /
  . o-x
    |\
  . x x


  The new marker pixel value is constrained to lie below the
  corresponding pixel in `mask`.
 */
ptrdiff_t backward_scan(float *marker, float *mask, ptrdiff_t dims[2]) {
  // Offsets for the four neighbors
  ptrdiff_t j_offset[4] = {1, 1, 1, 0};
  ptrdiff_t i_offset[4] = {-1, 0, 1, 1};

  ptrdiff_t count = 0;  // Number of modified pixels

  // Note that the loop decreases. p must have a signed type for this
  // to work correctly.
  for (ptrdiff_t j = dims[1] - 1; j >= 0; j--) {
    for (ptrdiff_t i = dims[0] - 1; i >= 0; i--) {
      ptrdiff_t p = j * dims[0] + i;

      // Compute the maximum of the marker at the current pixel and all
      // of its previously visited neighbors
      float max_height = marker[p];
      for (ptrdiff_t neighbor = 0; neighbor < 4; neighbor++) {
        ptrdiff_t neighbor_i = i + i_offset[neighbor];
        ptrdiff_t neighbor_j = j + j_offset[neighbor];
        ptrdiff_t q = neighbor_j * dims[0] + neighbor_i;

        // Skip pixels outside the boundary
        if (neighbor_i < 0 || neighbor_i >= dims[0] || neighbor_j < 0 ||
            neighbor_j >= dims[1]) {
          continue;
        }
        max_height = max_height > marker[q] ? max_height : marker[q];
      }

      // Set the marker at the current pixel to the minimum of the
      // maximum height of the neighborhood and the mask at the current
      // pixel.

      float z = max_height < mask[p] ? max_height : mask[p];
      if (z != marker[p]) {
        // Increment count only if we change the current pixel
        count++;
        marker[p] = z;
      }
    }
  }
  return count;
}

/*
  Grayscale reconstruction

  Performs a grayscale reconstruction of the `mask` image by the
  `marker` image using the sequential reconstruction algorithm of Vincent
  (1993).

  Both `marker` and `mask` should point to two-dimensional arrays of
  size (dims[0], dims[1]) with the first dimension (dims[0]) changing
  fastest. The `marker` array is updated with the result in-place.

  The algorithm alternately scans the marker image in the forward and
  reverse directions performing a partial reconstruction in each
  direction. It repeats these scans until no change is detected or
  until a maximum iteration threshold (currently 1000) is reached.

  Vincent, Luc. (1993). Morphological grayscale reconstruction in
  image analysis: applications and efficient algorithms. IEEE
  Transactions on Image Processing, Vol. 2, No. 2.
  https://doi.org/10.1109/83.217222
 */
void reconstruct(float *marker, float *mask, ptrdiff_t dims[2]) {
  ptrdiff_t n = dims[0] * dims[1];

  const int32_t max_iterations = 1000;
  for (int32_t iteration = 0; iteration < max_iterations && n > 0;
       iteration++) {
    n = forward_scan(marker, mask, dims);
    n += backward_scan(marker, mask, dims);
  }
}
