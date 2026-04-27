#define TOPOTOOLBOX_BUILD

#include "topotoolbox.h"

/////////////////////////////////////////
// Bit maps for storing weighted edge sets
//
// Bit maps are a combination of a uint8_t array to identify outgoing
// edges of each pixel and a float array of corresponding edge
// weights. The k-th of byte [i, j] in the array is set if the k-th
// outgoing edge is present. The edges are numbered like
//
// 3 | 2 | 1
// 4 | x | 0
// 5 | 6 | 7
//
// This is the ordering used by the Dinf algorithm. Note that in
// column-major, i increases down and j increases right. In row-major
// i increases right and j increases down. The bitmaps don't care
// about row or column major in the usual way of libtopotoolbox
// functions: reverse the order of the dims array if you are using
// row-major arrays.
//
// The edge weights are stored in an appropriately sized float
// array. The offset of each pixel's edge weights within this array is
// given by the running sum of the number of edges as you scan through
// the bitmap. The edges for each pixel are stored in the order of the
// edge numbering (counterclockwise from right).

// Count the number of bits set in a byte.
static int bitcount(uint8_t v) {
  int c;
  for (c = 0; v; c++) v &= v - 1;
  return c;
}

// Count all the edges in the bitmap. This is used to preallocate
// arrays of edges.
TOPOTOOLBOX_API
ptrdiff_t edgeset_count(uint8_t *bitmap, ptrdiff_t dims[2]) {
  ptrdiff_t n = 0;
  for (ptrdiff_t j = 0; j < dims[1]; j++) {
    for (ptrdiff_t i = 0; i < dims[0]; i++) {
      n += bitcount(bitmap[j * dims[0] + i]);
    }
  }
  return n;
}

// Prefix sum (scan) of the number of edges. The scan array contains
// offsets into the list of weights for each pixel. Returns the number
// of edges, so you can use scan_edges to prepare the scan array and
// count edges for allocating weights.
TOPOTOOLBOX_API
ptrdiff_t edgeset_scan(ptrdiff_t *scan, uint8_t *bitmap, ptrdiff_t dims[2]) {
  ptrdiff_t n = 0;

  for (ptrdiff_t j = 0; j < dims[1]; j++) {
    for (ptrdiff_t i = 0; i < dims[0]; i++) {
      scan[j * dims[0] + i] = n;
      n += bitcount(bitmap[j * dims[0] + i]);
    }
  }
  return n;
}

// Actually storing the weights in the edge array is the
// responsibility of the caller because they may need to compute the
// weights or other quantities.

// Count the total number of edges in the set merge of b1 and b2 This
// is used to preallocate the array of weights for the merged
// set. Duplicated edges are counted only once: this is not
// a multiset.
TOPOTOOLBOX_API
ptrdiff_t edgeset_count_merged(uint8_t *b1, uint8_t *b2, ptrdiff_t dims[2]) {
  ptrdiff_t n = 0;
  for (ptrdiff_t j = 0; j < dims[1]; j++) {
    for (ptrdiff_t i = 0; i < dims[0]; i++) {
      n += bitcount(b1[j * dims[0] + i] | b2[j * dims[0] + i]);
    }
  }
  return n;
}

// Merge the two edge sets (b1, w1) and (b2, w2). The edge weights are
// stored in the array w0 at the offsets in s0. If an edge exists in
// both b1 and b2, the weight from b1 is taken. The offset arrays of
// b1 and b2 are not necessary because the edge indices are computed
// along the way. The b1 array is updated with the new edges.
TOPOTOOLBOX_API
ptrdiff_t edgeset_merge(float *w0, ptrdiff_t *s0, uint8_t *b1, float *w1,
                        uint8_t *b2, float *w2, ptrdiff_t dims[2]) {
  ptrdiff_t n = 0;
  ptrdiff_t n1 = 0;
  ptrdiff_t n2 = 0;
  for (ptrdiff_t j = 0; j < dims[1]; j++) {
    for (ptrdiff_t i = 0; i < dims[0]; i++) {
      s0[j * dims[0] + i] = n;

      // Now loop over the neighbors
      for (ptrdiff_t k = 0; k < 8; k++) {
        if (b1[j * dims[0] + i] & (1 << k)) {
          w0[n++] = w1[n1++];

          // If the bit is also set in b2, skip that weight.
          if (b2[j * dims[0] + i] & (1 << k)) n2++;
        } else if (b2[j * dims[0] + i] & (1 << k)) {
          w0[n++] = w2[n2++];

          // Set the bit in b1
          b1[j * dims[0] + i] |= (1 << k);
        }
      }
    }
  }
  return n;
}
