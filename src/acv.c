#define TOPOTOOLBOX_BUILD

#include <math.h>
#include <stddef.h>
#include <stdint.h>

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

  ptrdiff_t size_of_dem = dims[0] * dims[1];
  // for 5x5 kernel
  ptrdiff_t k5_rows[5] = {-2, -1, 0, 1, 2};
  ptrdiff_t k5_cols[5] = {-2 * dims[0], -dims[0], 0, dims[0], 2 * dims[0]};
  // for 3x3 kernel
  ptrdiff_t k3_rows[3] = {-1, 0, 1};
  ptrdiff_t k3_cols[3] = {-dims[0], 0, dims[0]};

  // ACV:

  ptrdiff_t col;
#pragma omp parallel for if (use_mp)
  for (col = 0; col < dims[1]; col++) {
    for (ptrdiff_t row = 0; col < dims[0]; row++) {
      printf("row: %td, col: %td\n", row, col);

      float sum = 0.0f;
      float dz_avg = 0.0f;
      float anisotropic_cov = 0.0f;

      // filter_1
      for (ptrdiff_t k = 0; k < 25; k++) {
        if (filter_1[k == 0]) continue;
        ptrdiff_t to_dem_row = row + k5_cols[k / 5];
        ptrdiff_t to_dem_col = col + k5_cols[k % 5];
        printf("k: %td, row: %td, col: %td\n", k, to_dem_row, to_dem_col);

        // if out of bounds set dem value to closest cell in dem
        if (to_dem_row < 0) to_dem_row = 0;
        if (to_dem_row >= dims[0]) to_dem_row = dims[0] - 1;
        if (to_dem_col < 0) to_dem_col = 0;
        if (to_dem_col >= dims[1]) to_dem_col = dims[1] - 1;

        ptrdiff_t kernel_to_dem = to_dem_row * dims[0] + to_dem_col;
        printf("calc: %f  * %f", filter_1[k], dem[kernel_to_dem]);
        sum += filter_1[k] * dem[kernel_to_dem];
      }
      // dz_AVG  = conv2(dem,k,'valid')/4;
      dz_avg = sum / 4.0f;

      // ! temp break for testing
      break;

      // filter_2
      for (ptrdiff_t n = 0; n < 4; n++) {
        sum = 0.0f;
        for (ptrdiff_t k = 0; k < 25; k++) {
          if (filter_2[n][k == 0]) continue;
          ptrdiff_t to_dem_row = row + k5_cols[k / 5];
          ptrdiff_t to_dem_col = col + k5_cols[k % 5];

          // if out of bounds set dem value to closest cell in dem
          if (to_dem_row < 0) to_dem_row = 0;
          if (to_dem_row >= dims[0]) to_dem_row = dims[0] - 1;
          if (to_dem_col < 0) to_dem_col = 0;
          if (to_dem_col >= dims[1]) to_dem_col = dims[1] - 1;

          ptrdiff_t kernel_to_dem = to_dem_row * dims[0] + to_dem_col;
          sum += filter_2[n][k] * dem[kernel_to_dem];
        }
        // ACV = ACV + (conv2(dem,F{r},'valid') - dz_AVG).^2;
        anisotropic_cov += powf(sum - dz_avg, 2.0f);
      }

      // filter_3
      for (ptrdiff_t n = 0; n < 4; n++) {
        sum = 0.0f;
        for (ptrdiff_t k = 0; k < 9; k++) {
          if (filter_2[n][k == 0]) continue;
          ptrdiff_t to_dem_row = row + k3_cols[k / 3];
          ptrdiff_t to_dem_col = col + k3_cols[k % 3];

          // if out of bounds set dem value to closest cell in dem
          if (to_dem_row < 0) to_dem_row = 0;
          if (to_dem_row >= dims[0]) to_dem_row = dims[0] - 1;
          if (to_dem_col < 0) to_dem_col = 0;
          if (to_dem_col >= dims[1]) to_dem_col = dims[1] - 1;

          ptrdiff_t kernel_to_dem = to_dem_row * dims[0] + to_dem_col;

          sum += filter_2[n][k] * dem[kernel_to_dem];
        }
        // ACV = ACV + (conv2(dem,F{r},'valid') - dz_AVG).^2;
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
