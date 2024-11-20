#define TOPOTOOLBOX_BUILD

#include <math.h>
#include <stddef.h>
#include <stdint.h>
// TODO: remove stdio.h after testing
#include <stdio.h>

#if TOPOTOOLBOX_OPENMP_VERSION > 0
#include <omp.h>
#endif

#include "topotoolbox.h"

/*
The anisotropic coefficient of variation (ACV) describes the general
geometry of the local land surface and can be used to destinguish elongated
from oval land forms.

output:   Empty (or filled with zeros), same shape as dem
dem:      DEM matrix, passed from python in order='C'
use_mp:   If 1 then multiprocessing will be used to compute result
dims[2]:  dimensions of DEM

References:
    Olaya, V. 2009: Basic land-surface parameters. In: Geomorphometry.
    Concepts, Software, Applications, Hengl, T. & Reuter, H. I. (Eds.),
    Elsevier, 33, 141-169.

Author:
    Theophil Bringezu (theophil.bringezu[at]uni-potsdam.de)

Original MATLAB version by:
    Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)

*/

TOPOTOOLBOX_API
void acv(float *output, float *dem, int use_mp, ptrdiff_t dims[2]) {
  /*
  filter_1 :
  {{ 1, 0,  1, 0,  1},
   { 0, 0,  0, 0,  0},
   { 1, 0,  0, 0, -1},
   { 0, 0,  0, 0,  0},
   {-1, 0, -1, 0, -1}}
  */
  float filter_1[25] = {1, 0,  1, 0, -1, 0, 0, 0, 0, 0,  1, 0, 0,
                        0, -1, 0, 0, 0,  0, 0, 1, 0, -1, 0, -1};
  /* filter_2
  {0, 0,  0, 0,  0},
  {0, 0,  0, 0,  0},
  {1, 0,  0, 0, -1},
  {0, 0,  0, 0,  0},
  {0, 0,  0, 0,  0}},

 {{1, 0,  0, 0,  0},
  {0, 0,  0, 0,  0},
  {0, 0,  0, 0,  0},
  {0, 0,  0, 0,  0},
  {0, 0,  0, 0, -1}},

 {{0, 0, -1, 0,  0},
  {0, 0,  0, 0,  0},
  {0, 0,  0, 0,  0},
  {0, 0,  0, 0,  0},
  {0, 0, -1, 0,  0}},

 {{0, 0,  0, 0, -1},
  {0, 0,  0, 0,  0},
  {0, 0,  0, 0,  0},
  {0, 0,  0, 0,  0},
  {1, 0,  0, 0,  0}}}
  */
  float filter_2[4][25] = {{
                               0, 0, 1, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,
                               0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0,
                           },
                           {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1},
                           {
                               0, 0,  0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0,
                               0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0,
                           },
                           {0, 0, 0, 0, 1, 0, 0, 0,  0, 0, 0, 0, 0,
                            0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0}};
  /*
  filter_3
  {{{ 0,  0,  0},
    { 1,  0, -1},
    { 0,  0,  0}},

   {{ 1,  0,  0},
    { 0,  0,  0},
    { 0,  0, -1}},

   {{ 0,  1,  0},
    { 0,  0,  0},
    { 0, -1,  0}},

   {{ 0,  0,  1},
    { 0,  0,  0},
    {-1,  0,  0}}};
  */
  float filter_3[4][9] = {{0, 1, 0, 0, 0, 0, 0, -1, 0},
                          {1, 0, 0, 0, 0, 0, 0, 0, -1},
                          {0, 0, 0, 1, 0, -1, 0, 0, 0},
                          {0, 0, -1, 0, 0, 0, 1, 0, 0}};

  // ACV:
  ptrdiff_t col;
#pragma omp parallel for if (use_mp)
  for (col = 0; col < dims[1]; col++) {
    for (ptrdiff_t row = 0; row < dims[0]; row++) {
      float sum = 0.0f;
      float dz_avg = 0.0f;
      float anisotropic_cov = 0.0f;

      // filter_1
      for (ptrdiff_t k_col = -2; k_col <= 2; k_col++) {
        for (ptrdiff_t k_row = -2; k_row <= 2; k_row++) {
          // Add 2 to k_values to counteract the -2 start value
          ptrdiff_t k_index = (k_col + 2) * 5 + (k_row + 2);
          if (filter_1[k_index] == 0.0f) continue;

          // If row or col would be out of bounds, set it to appropriate border
          ptrdiff_t true_row = row + k_row;
          if (true_row < 0) true_row = 0;
          if (true_row >= dims[0]) true_row = dims[0] - 1;
          ptrdiff_t true_col = col + k_col;
          if (true_col < 0) true_col = 0;
          if (true_col >= dims[1]) true_col = dims[1] - 1;

          // Multiply kernel value by respective dem value
          ptrdiff_t true_index = true_row * dims[0] + true_col;
          sum += filter_1[k_index] * dem[true_index];
        }
      }
      // dz_AVG  = conv2(dem,k,'valid')/4;
      dz_avg = sum / 4.0f;

      // filter_2
      for (ptrdiff_t n = 0; n < 5; n++) {
        sum = 0.0f;
        for (ptrdiff_t k_col = -2; k_col <= 2; k_col++) {
          for (ptrdiff_t k_row = -2; k_row <= 2; k_row++) {
            ptrdiff_t k_index = (k_col + 2) * 5 + (k_row + 2);
            if (filter_2[n][k_index] == 0.0f) continue;

            ptrdiff_t true_row = row + k_row;
            if (true_row < 0) true_row = 0;
            if (true_row >= dims[0]) true_row = dims[0] - 1;
            ptrdiff_t true_col = col + k_col;
            if (true_col < 0) true_col = 0;
            if (true_col >= dims[1]) true_col = dims[1] - 1;

            ptrdiff_t true_index = true_row * dims[0] + true_col;
            sum += filter_2[n][k_index] * dem[true_index];
          }
        }
        anisotropic_cov += powf(sum - dz_avg, 2.0f);
      }

      // filter_3
      for (ptrdiff_t n = 0; n < 5; n++) {
        sum = 0.0f;
        for (ptrdiff_t k_col = -1; k_col <= 1; k_col++) {
          for (ptrdiff_t k_row = -1; k_row <= 1; k_row++) {
            ptrdiff_t k_index = (k_col + 2) * 5 + (k_row + 2);
            if (filter_3[n][k_index] == 0.0f) continue;

            ptrdiff_t true_row = row + k_row;
            if (true_row < 0) true_row = 0;
            if (true_row >= dims[0]) true_row = dims[0] - 1;
            ptrdiff_t true_col = col + k_col;
            if (true_col < 0) true_col = 0;
            if (true_col >= dims[1]) true_col = dims[1] - 1;

            ptrdiff_t true_index = true_row * dims[0] + true_col;
            sum += filter_3[n][k_index] * dem[true_index];
          }
        }
        anisotropic_cov += powf(sum - dz_avg, 2.0f);
      }
      // dz_AVG = max(abs(dz_AVG),0.001);
      dz_avg = fmaxf(0.001f, fabsf(dz_avg));

      // C = log(1 + sqrt(ACV./8)./dz_AVG);
      ptrdiff_t index = col * dims[0] + row;
      output[index] = logf(1.0f + sqrtf(anisotropic_cov / 8.0f) / dz_avg);
    }
  }
}
