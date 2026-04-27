#define TOPOTOOLBOX_BUILD

#include "topotoolbox.h"

#define SQRT2 1.41421356237309504880f

void flow_routing_d8_directions(uint8_t *direction, float *dem,
                                ptrdiff_t dims[2], int order) {
  // Basic D8 flow routing: flow to the maximum downstream neighbor until you
  // hit a flat or a sink.

  // These are the offsets for the appropriate neighbors.
  ptrdiff_t e[2][8] = {{0, -1, -1, -1, 0, 1, 1, 1},
                       {1, 1, 0, -1, -1, -1, 0, 1}};

  float chamfer[8] = {1.0f, SQRT2, 1.0f, SQRT2, 1.0f, SQRT2, 1.0f, SQRT2};

  for (ptrdiff_t j = 0; j < dims[1]; j++) {
    for (ptrdiff_t i = 0; i < dims[0]; i++) {
      float z = dem[j * dims[0] + i];
      float g = 0.0;
      direction[j * dims[0] + i] = 0;

      for (ptrdiff_t n = 0; n < 8; n++) {
        ptrdiff_t in = i + e[order & 1][n];
        ptrdiff_t jn = j + e[(order ^ 1) & 1][n];

        if (in < 0 || in >= dims[0] || jn < 0 || jn >= dims[1]) {
          continue;
        }

        float gn = (z - dem[jn * dims[0] + in]) / chamfer[n];

        if (gn > g) {
          direction[j * dims[0] + i] = 1 << n;
          g = gn;
        }
      }
    }
  }
}

void flow_routing_d8_weights(float *weight, ptrdiff_t count) {
  for (ptrdiff_t e = 0; e < count; e++) {
    weight[e] = 1.0;
  }
}
