#define TOPOTOOLBOX_BUILD

#include <math.h>
#include <stddef.h>
#include <stdint.h>

#include "topotoolbox.h"

/**
 * output: zeros like dem
 * dz_avg: zeros like dem
 * acv: zeros like dem
 * padded_dem: dem but padded with 2 rows/lines
 * dem: normal dem
 * dims: the dims
 */

TOPOTOOLBOX_API
void acv(float *output, float *dz_avg, float *anisotropic_cov,
         float *padded_dem, float *dem, ptrdiff_t dims[2]) {
  // average filter
  int filter[5][5] = {{1, 0, 1, 0, 1},
                      {0, 0, 0, 0, 0},
                      {1, 0, 0, 0, -1},
                      {0, 0, 0, 0, 0},
                      {-1, 0, -1, 0, -1}};
  ptrdiff_t i;
  for (i = 0; i < dims[0]; i++) {
    for (ptrdiff_t j = 0; j < dims[1]; j++) {
      ptrdiff_t location = i * dims[1] + j;
      // TODO: conv2 needs fix
      dz_avg[location] = conv2(padded_dem, filter, 'valid') / 4;
    }
  }

  int filter[4][5][5] = {{// filter 0
                          {0, 0, 0, 0, 0},
                          {0, 0, 0, 0, 0},
                          {1, 0, 0, 0, -1},
                          {0, 0, 0, 0, 0},
                          {0, 0, 0, 0, 0}},
                         {// filter 1
                          {1, 0, 0, 0, 0},
                          {0, 0, 0, 0, 0},
                          {0, 0, 0, 0, 0},
                          {0, 0, 0, 0, 0},
                          {0, 0, 0, 0, -1}},
                         {// filter 2
                          {0, 0, -1, 0, 0},
                          {0, 0, 0, 0, 0},
                          {0, 0, 0, 0, 0},
                          {0, 0, 0, 0, 0},
                          {0, 0, -1, 0, 0}},
                         {// filter 3
                          {0, 0, 0, 0, -1},
                          {0, 0, 0, 0, 0},
                          {0, 0, 0, 0, 0},
                          {0, 0, 0, 0, 0},
                          {1, 0, 0, 0, 0}}};

  for (int n = 0; i < 4; i++) {
    ptrdiff_t i;
    for (i = 0; i < dims[0]; i++) {
      for (ptrdiff_t j = 0; j < dims[1]; j++) {
        ptrdiff_t location = i * dims[1] + j;
        // TODO: conv2 needs fix
        anisotropic_cov[location] =
            anisotropic_cov[location] +
            pow(conv2(padded_dem, filter[n], 'valid') - dz_avg[location], 2.0f);
      }
    }
  }

  int filter[4][3][3] = {{// filter 0
                          {0, 0, 0},
                          {1, 0, -1},
                          {0, 0, 0}},
                         {// filter 1
                          {1, 0, 0},
                          {0, 0, 0},
                          {0, 0, -1}},
                         {// filter 2
                          {0, 1, 0},
                          {0, 0, 0},
                          {0, -1, 0}},
                         {// filter 3
                          {0, 0, 1},
                          {0, 0, 0},
                          {-1, 0, 0}}};

  for (int n = 0; i < 4; i++) {
    ptrdiff_t i;
    for (i = 0; i < dims[0]; i++) {
      for (ptrdiff_t j = 0; j < dims[1]; j++) {
        ptrdiff_t location = i * dims[1] + j;
        // TODO: conv2 need fix
        anisotropic_cov[location] =
            anisotropic_cov[location] +
            pow(conv2(dem, filter[n]) - dz_avg[location], 2.0f);
      }
    }
  }
  ptrdiff_t i;
  for (i = 0; i < dims[0]; i++) {
    for (ptrdiff_t j = 0; j < dims[1]; j++) {
      ptrdiff_t location = i * dims[1] + j;
      // TODO: max function need fix
      dz_avg[location] = maxf(fabsf(dz_avg[location]), 0.001f);
    }
  }
  ptrdiff_t i;
  for (i = 0; i < dims[0]; i++) {
    for (ptrdiff_t j = 0; j < dims[1]; j++) {
      ptrdiff_t location = i * dims[1] + j;
      output[location] =
          logf(1 + sqrtf(anisotropic_cov[location] / 8) / dz_avg[location]);
    }
  }
}