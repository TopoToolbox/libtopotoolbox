#include "reconstruct.h"

#include <assert.h>
#include <math.h>
#include <stddef.h>
#include <stdint.h>

// FIFO circular queue implementation
//
// The queue can hold a maximum of length - 1 elements.
typedef struct {
  ptrdiff_t *buffer;
  ptrdiff_t length;
  ptrdiff_t head;
  ptrdiff_t tail;
} PixelQueue;

// Add v to the queue.
// Returns 1 if successful, 0 if the queue is full.
int32_t enqueue(PixelQueue *q, ptrdiff_t v) {
  int32_t res = 0;
  ptrdiff_t next = (q->tail + 1) % q->length;
  if (next != q->head) {
    // Queue is not full
    q->buffer[q->tail] = v;
    q->tail = next;
    res = 1;
  }
  return res;
}

// Remove and return the first element of the queue.
//
// Returns a pointer to that element or NULL if the queue is empty.
ptrdiff_t *dequeue(PixelQueue *q) {
  ptrdiff_t *p = 0;
  if (q->head == q->tail) {
    // Queue is empty
    return p;
  }
  p = &q->buffer[q->head];
  q->head = (q->head + 1) % q->length;
  return p;
}

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
ptrdiff_t backward_scan(float *marker, PixelQueue *queue, float *mask,
                        ptrdiff_t dims[2]) {
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

      if (queue) {
        // Scan the neighborhood again to check if the pixel should be
        // added to the queue
        for (ptrdiff_t neighbor = 0; neighbor < 4; neighbor++) {
          ptrdiff_t neighbor_i = i + i_offset[neighbor];
          ptrdiff_t neighbor_j = j + j_offset[neighbor];
          ptrdiff_t q = neighbor_j * dims[0] + neighbor_i;

          // Skip pixels outside the boundary
          if (neighbor_i < 0 || neighbor_i >= dims[0] || neighbor_j < 0 ||
              neighbor_j >= dims[1]) {
            continue;
          }

          if (marker[q] < marker[p] && marker[q] < mask[q]) {
            if (enqueue(queue, p) == 0) {
              // In hybrid mode, we don't need to count the changes:
              // Instead, we use the return value to signal if we need
              // to repeat the scan because the queue filled up.
              return -1;
            };
          }
        }
      }
    }
  }
  return count;
}

/*
  Propagates changes via a breadth-first search of image.
 */
int32_t propagate(float *marker, PixelQueue *queue, float *mask,
                  ptrdiff_t dims[2]) {
  int32_t repeat_flag = 0;

  ptrdiff_t j_offset[8] = {-1, -1, -1, 0, 0, 1, 1, 1};
  ptrdiff_t i_offset[8] = {-1, 0, 1, -1, 1, -1, 0, 1};

  // p will be NULL if the queue is empty
  ptrdiff_t *p = dequeue(queue);
  while (p) {
    ptrdiff_t i = (*p) % dims[0];
    ptrdiff_t j = (*p) / dims[0];
    float pz = marker[*p];

    for (ptrdiff_t neighbor = 0; neighbor < 8; neighbor++) {
      ptrdiff_t neighbor_i = i + i_offset[neighbor];
      ptrdiff_t neighbor_j = j + j_offset[neighbor];
      ptrdiff_t q = neighbor_j * dims[0] + neighbor_i;

      if (neighbor_i < 0 || neighbor_i >= dims[0] || neighbor_j < 0 ||
          neighbor_j >= dims[1]) {
        continue;
      }

      if ((marker[q] < pz) && (mask[q] != marker[q])) {
        // Update the neighbor only if it meets the above criteria
        marker[q] = fminf(pz, mask[q]);

        if (enqueue(queue, q) == 0) {
          // If enqueuing a pixel fails because the queue is full, set
          // the repeat flag and skip this pixel.
          repeat_flag = 1;
        }
      }
    }

    p = dequeue(queue);
  }

  return repeat_flag;
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
    n += backward_scan(marker, NULL, mask, dims);
  }
}

/*
  Grayscale reconstruction using the hybrid algorithm of Vincent(1993).

  Requires a ptrdiff_t array of the same size as the marker and mask
  arrays for use as the backing store of a PixelQueue. The PixelQueue
  can only store dims[0]*dims[1] - 1 pixels. In the unlikely event
  that it fills up, either during the `backward_scan` step or the
  `propagate` step, the `repeat` flag is set, and the process is
  repeated to ensure that the algorithm converges.
 */
void reconstruct_hybrid(float *marker, ptrdiff_t *queue, float *mask,
                        ptrdiff_t dims[2]) {
  PixelQueue q = {0};
  q.buffer = queue;
  q.length = dims[0] * dims[1];

  int32_t repeat = 0;
  do {
    forward_scan(marker, mask, dims);
    // Backward scan returns -1 if the queue fills up, in which case
    // we need to repeat the scan.
    repeat |= (backward_scan(marker, &q, mask, dims) == -1);
    repeat |= propagate(marker, &q, mask, dims);
  } while (repeat);
}
