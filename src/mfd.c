#define TOPOTOOLBOX_BUILD

#include <assert.h>
#include <math.h>
#include <stddef.h>
#include <stdint.h>

#include "topotoolbox.h"

#define SQRT2 1.41421356237309504880f

TOPOTOOLBOX_API
void flow_routing_mfd_directions(uint8_t *direction, float *totalgradient,
                                 float *dem, ptrdiff_t dims[2], int order) {
  ptrdiff_t e[2][8] = {{0, -1, -1, -1, 0, 1, 1, 1},
                       {1, 1, 0, -1, -1, -1, 0, 1}};

  float d[8] = {1.0, SQRT2, 1.0, SQRT2, 1.0, SQRT2, 1.0, SQRT2};

  for (ptrdiff_t j = 0; j < dims[1]; j++) {
    for (ptrdiff_t i = 0; i < dims[0]; i++) {
      float z = dem[j * dims[0] + i];
      for (int n = 0; n < 8; n++) {
        ptrdiff_t i1 = i + e[order & 1][n];
        ptrdiff_t j1 = j + e[(order ^ 1) & 1][n];

        if (i1 < 0 || i1 >= dims[0] || j1 < 0 || j1 >= dims[1]) continue;

        float z2 = dem[j1 * dims[0] + i1];
        direction[j * dims[0] + i] |= (z > z2) << n;
        totalgradient[j * dims[0] + i] += fmaxf(0.0, z - z2) / d[n];
      }
    }
  }
}

TOPOTOOLBOX_API
void flow_routing_mfd_weights(float *weight, float *totalgradient,
                              uint8_t *direction, float *dem, ptrdiff_t dims[2],
                              int order) {
  ptrdiff_t e[2][8] = {{0, -1, -1, -1, 0, 1, 1, 1},
                       {1, 1, 0, -1, -1, -1, 0, 1}};

  float d[8] = {1.0, SQRT2, 1.0, SQRT2, 1.0, SQRT2, 1.0, SQRT2};

  ptrdiff_t edge = 0;
  for (ptrdiff_t j = 0; j < dims[1]; j++) {
    for (ptrdiff_t i = 0; i < dims[0]; i++) {
      ptrdiff_t idx = j * dims[0] + i;
      float z = dem[idx];
      float g = totalgradient[idx];

      for (int n = 0; n < 8; n++) {
        if (direction[idx] & (1 << n)) {
          ptrdiff_t i1 = i + e[order & 1][n];
          ptrdiff_t j1 = j + e[(order ^ 1) & 1][n];
          float z2 = dem[j1 * dims[0] + i1];

          weight[edge++] = ((z - z2) / d[n]) / g;
        }
      }
    }
  }
}
