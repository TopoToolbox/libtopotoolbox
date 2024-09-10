#define TOPOTOOLBOX_BUILD

#include <assert.h>
#include <stddef.h>
#include <stdint.h>

#include "topotoolbox.h"

TOPOTOOLBOX_API
void flow_accumulation_edgelist(float *acc, ptrdiff_t *source,
                                ptrdiff_t *target, float *fraction,
                                float *weights, ptrdiff_t edge_count,
                                ptrdiff_t dims[2]) {
  // Initialize with the weights
  for (ptrdiff_t j = 0; j < dims[1]; j++) {
    for (ptrdiff_t i = 0; i < dims[0]; i++) {
      acc[j * dims[0] + i] =
          (weights == NULL) ? 1.0f : weights[j * dims[0] + i];
    }
  }

  for (ptrdiff_t edge = 0; edge < edge_count; edge++) {
    ptrdiff_t src = source[edge];
    ptrdiff_t tgt = target[edge];
    float w = fraction[edge];

    // Ensure that src and tgt are valid pixels before accumulating.
    if (src >= 0 && src < dims[0] * dims[1] && tgt >= 0 &&
        tgt < dims[0] * dims[1]) {
      acc[tgt] += w * acc[src];
    }
  }
}

TOPOTOOLBOX_API
void flow_accumulation(float *acc, ptrdiff_t *source, uint8_t *direction,
                       float *weights, ptrdiff_t dims[2]) {
  ptrdiff_t i_offset[8] = {0, 1, 1, 1, 0, -1, -1, -1};
  ptrdiff_t j_offset[8] = {1, 1, 0, -1, -1, -1, 0, 1};

  // Initialize
  for (ptrdiff_t j = 0; j < dims[1]; j++) {
    for (ptrdiff_t i = 0; i < dims[0]; i++) {
      float w = (weights == NULL) ? 1.0f : weights[j * dims[0] + i];
      acc[j * dims[0] + i] = w;
    }
  }

  for (ptrdiff_t j = 0; j < dims[1]; j++) {
    for (ptrdiff_t i = 0; i < dims[0]; i++) {
      ptrdiff_t src = source[j * dims[0] + i];
      ptrdiff_t src_i = src % dims[0];
      ptrdiff_t src_j = src / dims[0];

      assert(src_i >= 0);
      assert(src_i < dims[0]);
      assert(src_j >= 0);
      assert(src_j < dims[1]);

      float w = acc[src_j * dims[0] + src_i];

      uint8_t flowdir = direction[src];
      if (flowdir == 0) {
        continue;
      }

      uint8_t v = flowdir;
      uint8_t r = 0;
      while (v >>= 1) {
        r++;
      }
      assert(r < 8);
      ptrdiff_t neighbor_i = src_i + i_offset[r];
      ptrdiff_t neighbor_j = src_j + j_offset[r];

      assert(neighbor_i >= 0);
      assert(neighbor_i < dims[0]);
      assert(neighbor_j >= 0);
      assert(neighbor_j < dims[1]);

      acc[neighbor_j * dims[0] + neighbor_i] += w;
    }
  }
}
