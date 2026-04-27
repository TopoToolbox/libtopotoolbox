#define TOPOTOOLBOX_BUILD

#include <assert.h>
#include <math.h>
#include <stddef.h>
#include <stdint.h>

#include "topotoolbox.h"

#define SQRT2 1.41421356237309504880f
#define PI 3.14159265358979323846f

static float replicate_boundaries(float *dem, ptrdiff_t i, ptrdiff_t j,
                                  ptrdiff_t ei, ptrdiff_t ej,
                                  ptrdiff_t dims[2]) {
  float z = dem[j * dims[0] + i];
  if (i + ei < 0 || i + ei >= dims[0]) {
    if (j + ej < 0 || j + ej >= dims[1]) {
      // Corners: z = dem[i, j]
    } else {
      // Top and bottom boundaries
      z = dem[(j + ej) * dims[0] + i];
    }
  } else if (j + ej < 0 || j + ej >= dims[1]) {
    // Left and right boundaries. The ei value is valid because we
    // skipped the previous branch.
    z = dem[j * dims[0] + i + ei];
  } else {
    // Interior
    z = dem[(j + ej) * dims[0] + i + ei];
  }
  return z;
}

typedef struct {
  int n;
  float a;
} facet;

static facet dinf_flow_prop(float v1, float v2, int order) {
  facet out = {-1, 0};
  if (v1 == 0.0 && v2 == 0.0) {
    return out;
  }
  // Compute the angle
  //
  // Despite appearances, this is correct. One must rotate the
  // coordinate axes to line up the (v1, v2) vector with the (ei, ej)
  // vectors. This depends on the memory order of the array. Then we
  // work in units of PI/4 so that we only have to do this division
  // once. We shift the angles from [-PI,PI] to [0,2PI] with the mod(x
  // + 8, 8) (8 is 2PI in units of PI/4).
  if (order) {
    out.a = fmodf(atan2f(-v2, v1) / (PI / 4) + 8, 8);
  } else {
    out.a = fmodf(atan2f(-v1, v2) / (PI / 4) + 8, 8);
  }

  // The result is a number between 0 and 8. The integer part is the
  // facet number. The fractional part is 1 - a, the correct flow
  // proportion proportion.
  out.n = (int)floorf(out.a);
  out.a = (out.n + 1) - out.a;

  return out;
}

TOPOTOOLBOX_API
void flow_routing_dinf_directions(uint8_t *direction, float *prop, float *dem,
                                  ptrdiff_t dims[2], int order) {
  // Compute and store the direction bitmap and the flow proportions
  // rather than the 2D flow vectors.
  ptrdiff_t e[2][8] = {{0, -1, -1, -1, 0, 1, 1, 1},
                       {1, 1, 0, -1, -1, -1, 0, 1}};

  for (ptrdiff_t j = 0; j < dims[1]; j++) {
    for (ptrdiff_t i = 0; i < dims[0]; i++) {
      ptrdiff_t idx0 = j * dims[0] + i;
      float z0 = dem[idx0];

      direction[idx0] = 0;
      prop[idx0] = 0.0;

      if (isnan(z0)) {
        prop[idx0] = NAN;
        continue;
      }

      float maxslope = 0.0;
      for (int n = 0; n < 8; n++) {
        ptrdiff_t i1 = e[order & 1][(n + (order & 1)) & 7];
        ptrdiff_t j1 = e[(order ^ 1) & 1][(n + (order & 1)) & 7];
        float z1 = replicate_boundaries(dem, i, j, i1, j1, dims);

        ptrdiff_t i2 = e[order & 1][(n + 1 - (order & 1)) & 7];
        ptrdiff_t j2 = e[(order ^ 1) & 1][(n + 1 - (order & 1)) & 7];
        float z2 = replicate_boundaries(dem, i, j, i2, j2, dims);

        float s1 = j1 * (z2 - z0) - j2 * (z1 - z0);
        float s2 = i1 * (z0 - z2) - i2 * (z0 - z1);

        float mag = sqrtf(s1 * s1 + s2 * s2);
        float d1 = sqrtf((float)(i1 * i1 + j1 * j1));
        float d2 = sqrtf((float)(i2 * i2 + j2 * j2));

        if (i1 * s2 - j1 * s1 < 0) {
          mag = (s1 * i1 + s2 * j1) / d1;
          s1 = mag * i1 / d1;
          s2 = mag * j1 / d1;
        } else if (j2 * s1 - i2 * s2 < 0) {
          mag = (s1 * i2 + s2 * j2) / d2;
          s1 = mag * i2 / d2;
          s2 = mag * j2 / d2;
        }

        if (mag > maxslope) {
          facet flowdir = dinf_flow_prop(s1, s2, order);

          prop[idx0] = flowdir.a;

          direction[idx0] = 0;
          if (flowdir.a > 0) {
            direction[idx0] |= 1 << (flowdir.n & 7);
          }

          if (flowdir.a < 1) {
            direction[idx0] |= 1 << ((flowdir.n + 1) & 7);
          }
          maxslope = mag;
        }
      }
    }
  }
}

TOPOTOOLBOX_API
void flow_routing_dinf_weights(float *weight, uint8_t *direction, float *prop,
                               ptrdiff_t dims[2]) {
  ptrdiff_t e = 0;
  for (ptrdiff_t j = 0; j < dims[1]; j++) {
    for (ptrdiff_t i = 0; i < dims[0]; i++) {
      ptrdiff_t idx = j * dims[0] + i;
      float a = prop[idx];

      if (direction[idx] == 129) {
        // Swap the order of the edges
        if (a < 1) {
          weight[e++] = 1 - a;
        }

        if (a > 0) {
          weight[e++] = a;
        }
      } else if (direction[idx] > 0) {
        if (a > 0) {
          weight[e++] = a;
        }

        if (a < 1) {
          weight[e++] = 1 - a;
        }
      }
    }
  }
}
