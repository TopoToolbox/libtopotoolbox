#define TOPOTOOLBOX_BUILD

#include <math.h>
#include <stddef.h>
#include <stdint.h>

#if TOPOTOOLBOX_OPENMP_VERSION > 0
#include <omp.h>
#endif

#include "topotoolbox.h"

/**
 * output: zeros like dem
 * dz_avg: zeros like dem
 * acv: zeros like dem
 * dem: dem matrix with shape of dims[0] x dims[1]
 * dims: the dims
 */

TOPOTOOLBOX_API
void acv(float *output, float *dz_avg, float *anisotropic_cov, float *dem,
         ptrdiff_t dims[2]) {
  float filter_1[5][5] = {{1, 0, 1, 0, 1},
                          {0, 0, 0, 0, 0},
                          {1, 0, 0, 0, -1},
                          {0, 0, 0, 0, 0},
                          {-1, 0, -1, 0, -1}};

  float filter_2[4][5][5] = {{// filter 0
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

  float filter_3[4][3][3] = {{// filter 0
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

  // Loop through all pixels, the value of each pixel will be calculated based
  // on the results of previous filter applications. So the result for each
  // cell is not dependent on the result of any other cell.
  ptrdiff_t i;
  for (i = 0; i < dims[0]; i++) {
    for (ptrdiff_t j = 0; j < dims[1]; j++) {
      ptrdiff_t location = i * dims[1] + j;
      float sum = 0.0;
      // Filter 1 : Apply 5x5 filter, therefor starting at -2 and ending at 2
      for (ptrdiff_t m = -2; m <= 2; m++) {
        for (ptrdiff_t n = -2; n <= 2; n++) {
          // TODO: ensure the right borders are checked here
          if (m + i < 0 || n + j < 0 || m + i >= dims[0] || n + j >= dims[1]) {
            continue;
          }
          sum += filter_1[m][n] * dem[(i + m) * dims[1] + (j + n)];
        }
      }
      // dz_AVG  = conv2(dem,k,'valid')/4;
      dz_avg[location] = sum / 4;

      // Filter 2 : Apply all four 5x5 filters, 'f_num' to index filters
      for (ptrdiff_t f_num = 0; f_num < 4; f_num++) {
        // starting at '-2' and ending at '2' since it's a 5x5 filter
        for (ptrdiff_t m = -2; m <= 2; m++) {
          for (ptrdiff_t n = -2; n <= 2; n++) {
            // TODO: ensure the right borders are checked here
            if (m + i < 0 || n + j < 0 || m + i >= dims[0] ||
                n + j >= dims[1]) {
              continue;
            }
            sum += filter_2[f_num][m][n] * dem[(i + m) * dims[1] + (j + n)];
          }
        }
        // ACV = ACV + (conv2(dem,F{r},'valid') - dz_AVG).^2;
        anisotropic_cov[location] += pow(sum - dz_avg[location], 2.0f);
      }
      // Filter 3 : Apply all four 3x3 filters, 'f_num' to index filters
      for (ptrdiff_t f_num = 0; f_num < 4; f_num++) {
        sum = 0.0;
        // starting at '-1' and ending at '1' since it's a 3x3 filter
        for (ptrdiff_t m = -1; m <= 1; m++) {
          for (ptrdiff_t n = -1; n <= 1; n++) {
            // TODO: ensure the right borders are checked here
            if (m + i < 0 || n + j < 0 || m + i >= dims[0] ||
                n + j >= dims[1]) {
              continue;  // Skip out-of-bound pixels (same as += 0)
            }
            sum += filter_3[f_num][m + 1][n + 1] *
                   dem[(i + m) * dims[1] + (j + n)];
          }
        }
        // Sum up the results of each filter layer
        // ACV = ACV + (conv2(dem,F{r},'valid') - dz_AVG).^2;
        anisotropic_cov[location] += pow(sum - dz_avg[location], 2.0f);
      }

      // dz_AVG = max(abs(dz_AVG),0.001);
      if (fabsf(dz_avg[location]) > 0.001f) {
        dz_avg[location] = fabsf(dz_avg[location]);
      } else {
        dz_avg[location] = 0.001f;
      }

      // C = log(1 + sqrt(ACV./8)./dz_AVG);
      output[location] = logf(1.0f + sqrtf(anisotropic_cov[location] / 8.0f) /
                                         dz_avg[location]);
    }
  }
}
