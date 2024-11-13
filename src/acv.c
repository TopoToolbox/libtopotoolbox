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

  ptrdiff_t i;
#pragma omp parallel for if (use_mp)
  for (i = 0; i < size_of_dem; i++) {
    ptrdiff_t row = i / dims[0];  // automatic floor because of datatype
    ptrdiff_t col = i % dims[0];

    float sum = 0.0f;
    float dz_avg = 0.0f;
    float anisotropic_cov = 0.0f;

    // filter_1
    for (ptrdiff_t k = 0; k < 25; k++) {
      ptrdiff_t kernel_to_dem = i + k5_rows[k / 25] + k5_cols[k % 25];
      // TODO: skips if out of bounds (assumes bounds are zero)
      if (kernel_to_dem < 0 || kernel_to_dem > size_of_dem) continue;
      if (filter_1[k == 0]) continue;
      sum += filter_1[k] * dem[kernel_to_dem];
    }
    // dz_AVG  = conv2(dem,k,'valid')/4;
    dz_avg = sum / 4.0f;

    // filter_2
    for (ptrdiff_t n = 0; n < 4; n++) {
      sum = 0.0f;
      for (ptrdiff_t k = 0; k < 25; k++) {
        ptrdiff_t kernel_to_dem = i + k5_rows[k / 25] + k5_cols[k % 25];
        // TODO: skips if out of bounds (assumes bounds are zero)
        if (kernel_to_dem < 0 || kernel_to_dem >= size_of_dem) continue;
        if (filter_2[n][k == 0]) continue;
        sum += filter_2[n][k] * dem[kernel_to_dem];
      }
      // ACV = ACV + (conv2(dem,F{r},'valid') - dz_AVG).^2;
      anisotropic_cov += powf(sum - dz_avg, 2.0f);
    }

    // filter_3
    for (ptrdiff_t n = 0; n < 4; n++) {
      sum = 0.0f;
      for (ptrdiff_t k = 0; k < 9; k++) {
        ptrdiff_t kernel_to_dem = i + k3_rows[k / 9] + k3_cols[k % 9];
        // TODO: skips if out of bounds (assumes bounds are zero)
        if (kernel_to_dem < 0 || kernel_to_dem >= size_of_dem) continue;
        if (filter_2[n][k == 0]) continue;
        sum += filter_2[n][k] * dem[kernel_to_dem];
      }
      // ACV = ACV + (conv2(dem,F{r},'valid') - dz_AVG).^2;
      anisotropic_cov += powf(sum - dz_avg, 2.0f);
    }
    // dz_AVG = max(abs(dz_AVG),0.001);
    dz_avg = fmaxf(0.001f, fabsf(dz_avg));

    // C = log(1 + sqrt(ACV./8)./dz_AVG);
    output[i] = logf(1.0f + sqrtf(anisotropic_cov / 8.0f) / dz_avg);
  }
}

/*
  ptrdiff_t i;
#pragma omp parallel for if (use_mp)
  for (i = 0; i < dims[0]; i++) {
    for (ptrdiff_t j = 0; j < dims[1]; j++) {
      ptrdiff_t location = i * dims[1] + j;
      float sum = 0.0;
      float dz_avg = 0.0;
      float anisotropic_cov = 0.0;

      // Filter 1 : Apply 5x5 filter
      for (ptrdiff_t m = 0; m < 5; m++) {
        for (ptrdiff_t n = 0; n < 5; n++) {
          if (m + i - 2 < 0 || n + j - 2 < 0 || m + i - 2 >= dims[0] ||
              n + j - 2 >= dims[1]) {
            continue;
          }
          if (filter_1[m][n] == 0) {
            continue;
          }
          sum += filter_1[m][n] * dem[(i + m - 2) * dims[1] + (j + n - 2)];
        }
      }
      // dz_AVG  = conv2(dem,k,'valid')/4;
      dz_avg = sum / 4;

      // Filter 2 : Apply all four 5x5 filters, 'f_num' to index filters
      for (ptrdiff_t f_num = 0; f_num < 4; f_num++) {
        sum = 0.0;
        for (ptrdiff_t m = 0; m < 5; m++) {
          for (ptrdiff_t n = 0; n < 5; n++) {
            if (m + i - 2 < 0 || n + j - 2 < 0 || m + i - 2 >= dims[0] ||
                n + j - 2 >= dims[1]) {
              continue;
            }
            if (filter_2[f_num][m][n] == 0) {
              continue;
            }
            sum += filter_2[f_num][m][n] *
                   dem[(i + m - 2) * dims[1] + (j + n - 2)];
          }
        }
        // ACV = ACV + (conv2(dem,F{r},'valid') - dz_AVG).^2;
        anisotropic_cov += powf(sum - dz_avg, 2.0f);
      }
      // Filter 3 : Apply all four 3x3 filters, 'f_num' to index filters
      for (ptrdiff_t f_num = 0; f_num < 4; f_num++) {
        sum = 0.0;
        for (ptrdiff_t m = 0; m < 3; m++) {
          for (ptrdiff_t n = 0; n < 3; n++) {
            if (m + i - 1 < 0 || n + j - 1 < 0 || m + i - 1 >= dims[0] ||
                n + j - 1 >= dims[1]) {
              continue;
            }
            if (filter_3[f_num][m][n] == 0) {
              continue;
            }
            sum += filter_3[f_num][m][n] *
                   dem[(i + m - 1) * dims[1] + (j + n - 1)];
          }
        }
        // ACV = ACV + (conv2(dem,F{r},'valid') - dz_AVG).^2;
        anisotropic_cov += powf(sum - dz_avg, 2.0f);
      }
      // dz_AVG = max(abs(dz_AVG),0.001);
      dz_avg = fmaxf(0.001f, fabsf(dz_avg));

      // C = log(1 + sqrt(ACV./8)./dz_AVG);
      output[location] = logf(1.0f + sqrtf(anisotropic_cov / 8.0f) / dz_avg);
    }
  }
}*/
