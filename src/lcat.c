#define TOPOTOOLBOX_BUILD

#include "topotoolbox.h"

#define SQRT2 1.41421356237309504880f

void resolve_flats_lcat(uint8_t *direction, uint8_t *resolved, float *aux,
                        float *dem, ptrdiff_t dims[2], int order) {
  // Least cost auxiliary topography carving. The hard part is
  // computing the auxiliary topography in aux, but you do that
  // first. Otherwise this is D8 flow routing over unresolved pixels
  // using the auxiliary topography.
  ptrdiff_t e[2][8] = {{0, -1, -1, -1, 0, 1, 1, 1},
                       {1, 1, 0, -1, -1, -1, 0, 1}};

  float chamfer[8] = {1.0f, SQRT2, 1.0f, SQRT2, 1.0f, SQRT2, 1.0f, SQRT2};

  for (ptrdiff_t j = 0; j < dims[1]; j++) {
    for (ptrdiff_t i = 0; i < dims[0]; i++) {
      if (resolved[j * dims[0] + i]) continue;

      float z = dem[j * dims[0] + i];
      float t = aux[j * dims[0] + i];
      float g = 0.0;

      for (ptrdiff_t n = 0; n < 8; n++) {
        ptrdiff_t in = i + e[order & 1][n];
        ptrdiff_t jn = j + e[(order ^ 1) & 1][n];

        if (in < 0 || in >= dims[0] || jn < 0 || jn >= dims[1]) {
          continue;
        }

        // Skip pixels that are greater than the current elevation
        if (dem[jn * dims[0] + in] > z) {
          continue;
        }

        float gn = (t - aux[jn * dims[0] + in]) / chamfer[n];

        if (gn > g) {
          direction[j * dims[0] + i] = 1 << n;
          g = gn;
        }
      }
    }
  }
}

void resolve_flats_lcat_weights(float *weight, ptrdiff_t count) {
  for (ptrdiff_t e = 0; e < count; e++) {
    weight[e] = 1.0;
  }
}
