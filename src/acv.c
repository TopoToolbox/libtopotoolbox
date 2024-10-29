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

output:   zeros like DEM
dz_avg:   zeros like DEM
acv:      zeros like DEm
dem:      DEM matrix
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
  /*float filter_1[5][5] = {{1, 0, 1, 0, 1},
                          {0, 0, 0, 0, 0},
                          {1, 0, 0, 0, -1},
                          {0, 0, 0, 0, 0},
                          {-1, 0, -1, 0, -1}};*/

  float filter_1[5][5] = {{1, 2, 3, 4, 5},
                          {6, 7, 8, 9, 10},
                          {11, 12, 13, 14, 15},
                          {16, 17, 18, 19, 20},
                          {21, 22, 23, 24, 25}};

  float filter_2[4][5][5] = {{// f_num 0
                              {0, 0, 0, 0, 0},
                              {0, 0, 0, 0, 0},
                              {1, 0, 0, 0, -1},
                              {0, 0, 0, 0, 0},
                              {0, 0, 0, 0, 0}},
                             {// f_num 1
                              {1, 0, 0, 0, 0},
                              {0, 0, 0, 0, 0},
                              {0, 0, 0, 0, 0},
                              {0, 0, 0, 0, 0},
                              {0, 0, 0, 0, -1}},
                             {// f_num 2
                              {0, 0, -1, 0, 0},
                              {0, 0, 0, 0, 0},
                              {0, 0, 0, 0, 0},
                              {0, 0, 0, 0, 0},
                              {0, 0, -1, 0, 0}},
                             {// f_num 3
                              {0, 0, 0, 0, -1},
                              {0, 0, 0, 0, 0},
                              {0, 0, 0, 0, 0},
                              {0, 0, 0, 0, 0},
                              {1, 0, 0, 0, 0}}};

  float filter_3[4][3][3] = {{// f_num 0
                              {0, 0, 0},
                              {1, 0, -1},
                              {0, 0, 0}},
                             {// f_num 1
                              {1, 0, 0},
                              {0, 0, 0},
                              {0, 0, -1}},
                             {// f_num 2
                              {0, 1, 0},
                              {0, 0, 0},
                              {0, -1, 0}},
                             {// f_num 3
                              {0, 0, 1},
                              {0, 0, 0},
                              {-1, 0, 0}}};

  ptrdiff_t i;
#pragma omp parallel for if (use_mp)
  for (i = 0; i < dims[0]; i++) {
    for (ptrdiff_t j = 0; j < dims[1]; j++) {
      ptrdiff_t location = i * dims[1] + j;
      float sum = 0.0;
      float dz_avg = 0.0;
      float anisotropic_cov = 0.0;
      printf("Cell value: i=%td, j=%td -> %f \n", i, j, dem[location]);

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

          printf("Filter value: m=%td, n=%td -> %f\n", m, n, filter_1[m][n]);
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
            printf("[%td][%td] + %f -> sum: %f\n", i + m, j + n,
                   filter_2[f_num][m][n] *
                       dem[(i + m - 2) * dims[1] + (j + n - 2)],
                   sum);
          }
        }
        // ACV = ACV + (conv2(dem,F{r},'valid') - dz_AVG).^2;
        anisotropic_cov += pow(sum - dz_avg, 2.0f);
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
        anisotropic_cov += pow(sum - dz_avg, 2.0f);
      }

      // TODO: maybe use fmaxf() function
      // dz_AVG = max(abs(dz_AVG),0.001);
      if (fabsf(dz_avg) > 0.001f) {
        dz_avg = fabsf(dz_avg);
      } else {
        dz_avg = 0.001f;
      }

      // C = log(1 + sqrt(ACV./8)./dz_AVG);
      output[location] = logf(1.0f + sqrtf(anisotropic_cov / 8.0f) / dz_avg);
    }
  }
}
