/**
   @file filter.cpp
   @brief This module contains property-based tests for value filters
   and morphological filters. Descriptions of individual tests
   are given at their definition.
 */

// TODO: no tests with real 3d SEs performed!

#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <iostream>

namespace tt {
extern "C" {
#include "topotoolbox.h"
}
}  // namespace tt

namespace utils {
extern "C" {
#include "utils.h"
}
}  // namespace utils

/**
 * @brief Test that the min_filter function doesn't result
   in values larger than input cells.

   @note Structuring element is of size 3 by 3 and filled with ones
 *
 * @return int
 */
int test_no_larger_values(float *dem, float *output, uint8_t se[9],
                          ptrdiff_t dims[2], ptrdiff_t se_dims[3]) {
  tt::min_filter(output, dem, se, dims, se_dims);

  for (ptrdiff_t row = 0; row < dims[1]; row++) {
    for (ptrdiff_t col = 0; col < dims[0]; col++) {
      ptrdiff_t index = col + row * dims[0];

      if (std::isnan(output[index]) or std::isnan(dem[index])) {
        continue;
      }

      if (output[index] > dem[index]) {
        std::cout << "(" << col << ", " << row << "): " << output[index]
                  << " > " << dem[index] << std::endl;
        assert(0);
      }
    }
  }

  return 0;
}

/**
 * @brief Test that the min_filter_square function doesn't result
   in values larger than input cells.

   @note Structuring element is of size 3 by 3 and filled with ones
 *
 * @return int
 */
int test_no_larger_values_square(float *dem, float *output, float *tmp,
                                 uint8_t se_width, ptrdiff_t dims[2]) {
  tt::min_filter_square(output, dem, tmp, se_width, dims);

  for (ptrdiff_t row = 0; row < dims[1]; row++) {
    for (ptrdiff_t col = 0; col < dims[0]; col++) {
      ptrdiff_t index = col + row * dims[0];

      if (std::isnan(output[index]) or std::isnan(dem[index])) {
        continue;
      }

      if (output[index] > dem[index]) {
        std::cout << "(" << col << ", " << row << "): " << output[index]
                  << " > " << dem[index] << std::endl;
        assert(0);
      }
    }
  }
  return 0;
}

/**
 * @brief Test that the max_filter function doesn't result
   in values smaller than input cells.

     @note Structuring element is of size 3 by 3 and filled with ones
 *
 * @return int
 */
int test_no_smaller_values(float *dem, float *output, uint8_t se[9],
                           ptrdiff_t dims[2], ptrdiff_t se_dims[3]) {
  tt::max_filter(output, dem, se, dims, se_dims);

  for (ptrdiff_t row = 0; row < dims[1]; row++) {
    for (ptrdiff_t col = 0; col < dims[0]; col++) {
      ptrdiff_t index = col + row * dims[0];

      if (std::isnan(output[index]) or std::isnan(dem[index])) {
        continue;
      }

      if (output[index] < dem[index]) {
        std::cout << "(" << col << ", " << row << "): " << output[index]
                  << " < " << dem[index] << std::endl;
        assert(0);
      }
    }
  }

  return 0;
}

/**
 * @brief Test that the max_filter_square function doesn't result
   in values smaller than input cells.

     @note Structuring element is of size 3 by 3 and filled with ones
 *
 * @return int
 */
int test_no_smaller_values_square(float *dem, float *output, float *tmp,
                                  uint8_t se_width, ptrdiff_t dims[2]) {
  tt::max_filter_square(output, dem, tmp, se_width, dims);

  for (ptrdiff_t row = 0; row < dims[1]; row++) {
    for (ptrdiff_t col = 0; col < dims[0]; col++) {
      ptrdiff_t index = col + row * dims[0];

      if (std::isnan(output[index]) or std::isnan(dem[index])) {
        continue;
      }

      if (output[index] < dem[index]) {
        std::cout << "(" << col << ", " << row << "): " << output[index]
                  << " < " << dem[index] << std::endl;
        assert(0);
      }
    }
  }
  return 0;
}

/**
 * @brief Test that NAN values persist in min_filter

   @note Structuring element is of size 3 by 3 and filled with ones
 *
 * @return int
 */
int test_nan_values_persist_min_filter(float *dem, float *output, uint8_t se[9],
                                       ptrdiff_t dims[2],
                                       ptrdiff_t se_dims[3]) {
  tt::min_filter(output, dem, se, dims, se_dims);

  for (ptrdiff_t row = 0; row < dims[1]; row++) {
    for (ptrdiff_t col = 0; col < dims[0]; col++) {
      ptrdiff_t index = col + row * dims[0];

      if (std::isnan(output[index]) != std::isnan(dem[index])) {
        std::cout << "(" << col << ", " << row << "): isnan(" << output[index]
                  << ") != isnan(" << dem[index] << ")" << std::endl;
        assert(0);
      }
    }
  }

  return 0;
}

/**
 * @brief Test that NAN values persist in min_filter_square

   @note Structuring element is of size 3 by 3 and filled with ones
 *
 * @return int
 */
int test_nan_values_persist_min_filter_square(float *dem, float *output,
                                              float *tmp, uint8_t se_width,
                                              ptrdiff_t dims[2]) {
  tt::min_filter_square(output, dem, tmp, se_width, dims);

  for (ptrdiff_t row = 0; row < dims[1]; row++) {
    for (ptrdiff_t col = 0; col < dims[0]; col++) {
      ptrdiff_t index = col + row * dims[0];

      if (std::isnan(output[index]) != std::isnan(dem[index])) {
        std::cout << "(" << col << ", " << row << "): isnan(" << output[index]
                  << ") != isnan(" << dem[index] << ")" << std::endl;
        assert(0);
      }
    }
  }

  return 0;
}

/**
 * @brief Test that NAN values persist in max_filter

   @note Structuring element is of size 3 by 3 and filled with ones
 *
 * @return int
 */
int test_nan_values_persist_max_filter(float *dem, float *output, uint8_t se[9],
                                       ptrdiff_t dims[2],
                                       ptrdiff_t se_dims[3]) {
  tt::max_filter(output, dem, se, dims, se_dims);

  for (ptrdiff_t row = 0; row < dims[1]; row++) {
    for (ptrdiff_t col = 0; col < dims[0]; col++) {
      ptrdiff_t index = col + row * dims[0];

      if (std::isnan(output[index]) != std::isnan(dem[index])) {
        std::cout << "(" << col << ", " << row << "): isnan(" << output[index]
                  << ") != isnan(" << dem[index] << ")" << std::endl;
        assert(0);
      }
    }
  }

  return 0;
}

/**
 * @brief Test that NAN values persist in max_filter_square

   @note Structuring element is of size 3 by 3 and filled with ones
 *
 * @return int
 */
int test_nan_values_persist_max_filter_square(float *dem, float *output,
                                              float *tmp, uint8_t se_width,
                                              ptrdiff_t dims[2]) {
  tt::max_filter_square(output, dem, tmp, se_width, dims);

  for (ptrdiff_t row = 0; row < dims[1]; row++) {
    for (ptrdiff_t col = 0; col < dims[0]; col++) {
      ptrdiff_t index = col + row * dims[0];

      if (std::isnan(output[index]) != std::isnan(dem[index])) {
        std::cout << "(" << col << ", " << row << "): isnan(" << output[index]
                  << ") != isnan(" << dem[index] << ")" << std::endl;
        assert(0);
      }
    }
  }

  return 0;
}

/**
 * @brief Test that min_filter and min_filter_square result in same output for
 * identical structuring elements

   @note Structuring element is of size 3 by 3 and filled with ones
 *
 * @return int
 */
int test_min_filter_implementations_agree(float *dem, float *foutput,
                                          float *soutput, float *tmp,
                                          uint8_t se_width, uint8_t se[9],
                                          ptrdiff_t dims[2],
                                          ptrdiff_t se_dims[3]) {
  tt::min_filter(foutput, dem, se, dims, se_dims);

  tt::min_filter_square(soutput, dem, tmp, se_width, dims);

  for (ptrdiff_t row = 0; row < dims[1]; row++) {
    for (ptrdiff_t col = 0; col < dims[0]; col++) {
      ptrdiff_t index = col + row * dims[0];

      if ((std::isnan(foutput[index]) != std::isnan(soutput[index])) ||
          // second isnan check not needed but added for completeness
          (!std::isnan(foutput[index]) && !std::isnan(soutput[index]) &&
           foutput[index] != soutput[index])) {
        std::cout << "max rows: " << dims[1] << ", max cols: " << dims[0]
                  << std::endl;
        std::cout << "(" << col << ", " << row << "): " << foutput[index]
                  << " != " << soutput[index] << std::endl;
        assert(0);
      }
    }
  }
  return 0;
}

/**
 * @brief Test that max_filter and max_filter_square result in same output for
 * identical structuring elements

   @note Structuring element is of size 3 by 3 and filled with ones
 *
 * @return int
 */
int test_max_filter_implementations_agree(float *dem, float *foutput,
                                          float *soutput, float *tmp,
                                          uint8_t se_width, uint8_t se[9],
                                          ptrdiff_t dims[2],
                                          ptrdiff_t se_dims[3]) {
  tt::max_filter(foutput, dem, se, dims, se_dims);

  tt::max_filter_square(soutput, dem, tmp, se_width, dims);

  for (ptrdiff_t row = 0; row < dims[1]; row++) {
    for (ptrdiff_t col = 0; col < dims[0]; col++) {
      ptrdiff_t index = col + row * dims[0];

      if ((std::isnan(foutput[index]) != std::isnan(soutput[index])) ||
          // second isnan check not needed but added for completeness
          (!std::isnan(foutput[index]) && !std::isnan(soutput[index]) &&
           foutput[index] != soutput[index])) {
        std::cout << "(" << col << ", " << row << "): " << foutput[index]
                  << " != " << soutput[index] << std::endl;
        assert(0);
      }
    }
  }
  return 0;
}

/**
 * @brief Test that min_filter and morphological_erosion agree when 0 = NAN and
 * 1 = 0

   @note Structuring element is of size 3 by 3 and filled with ones
 *
 * @return int
 */
int test_min_filter_and_erosion_agree(float *dem, float *foutput,
                                      float *soutput, uint8_t se[9],
                                      ptrdiff_t dims[3], ptrdiff_t se_dims[2],
                                      float *float_se) {
  tt::min_filter(foutput, dem, se, dims, se_dims);

  tt::morphological_erosion(soutput, dem, float_se, dims, se_dims);

  for (ptrdiff_t row = 0; row < dims[1]; row++) {
    for (ptrdiff_t col = 0; col < dims[0]; col++) {
      ptrdiff_t index = col + row * dims[0];

      if ((std::isnan(foutput[index]) != std::isnan(soutput[index])) ||
          // second isnan check not needed but added for completeness
          (!std::isnan(foutput[index]) && !std::isnan(soutput[index]) &&
           foutput[index] != soutput[index])) {
        std::cout << "(" << col << ", " << row << "): " << foutput[index]
                  << " != " << soutput[index] << std::endl;
        assert(0);
      }
    }
  }

  return 0;
}

/**
 * @brief Test that min_filter_square and morphological_erosion agree when 0 =
 * NAN and 1 = 0

   @note Structuring element is of size 3 by 3 and filled with ones
 *
 * @return int
 */
int test_min_square_filter_and_erosion_agree(
    float *dem, float *foutput, float *soutput, float *tmp, uint8_t se_width,
    ptrdiff_t dims[3], ptrdiff_t se_dims[2], float *float_se) {
  tt::min_filter_square(foutput, dem, tmp, se_width, dims);

  tt::morphological_erosion(soutput, dem, float_se, dims, se_dims);

  for (ptrdiff_t row = 0; row < dims[1]; row++) {
    for (ptrdiff_t col = 0; col < dims[0]; col++) {
      ptrdiff_t index = col + row * dims[0];

      if ((std::isnan(foutput[index]) != std::isnan(soutput[index])) ||
          // second isnan check not needed but added for completeness
          (!std::isnan(foutput[index]) && !std::isnan(soutput[index]) &&
           foutput[index] != soutput[index])) {
        std::cout << "(" << col << ", " << row << "): " << foutput[index]
                  << " != " << soutput[index] << std::endl;
        assert(0);
      }
    }
  }

  return 0;
}

/**
 * @brief Test that max_filter and morphological_dilation agree when 0 = NAN and
 * 1 = 0

   @note Structuring element is of size 3 by 3 and filled with ones
 *
 * @return int
 */
int test_max_filter_and_dilation_agree(float *dem, float *foutput,
                                       float *soutput, uint8_t se[9],
                                       ptrdiff_t dims[3], ptrdiff_t se_dims[2],
                                       float *float_se) {
  tt::max_filter(foutput, dem, se, dims, se_dims);

  tt::morphological_dilation(soutput, dem, float_se, dims, se_dims);

  for (ptrdiff_t row = 0; row < dims[1]; row++) {
    for (ptrdiff_t col = 0; col < dims[0]; col++) {
      ptrdiff_t index = col + row * dims[0];

      if ((std::isnan(foutput[index]) != std::isnan(soutput[index])) ||
          // second isnan check not needed but added for completeness
          (!std::isnan(foutput[index]) && !std::isnan(soutput[index]) &&
           foutput[index] != soutput[index])) {
        std::cout << "(" << col << ", " << row << "): " << foutput[index]
                  << " != " << soutput[index] << std::endl;
        assert(0);
      }
    }
  }

  return 0;
}

/**
 * @brief Test that max_filter_square and morphological_dilation agree when 0 =
 * NAN and 1 = 0

   @note Structuring element is of size 3 by 3 and filled with ones
 *
 * @return int
 */
int test_max_square_filter_and_dilation_agree(
    float *dem, float *foutput, float *soutput, float *tmp, uint8_t se_width,
    ptrdiff_t dims[3], ptrdiff_t se_dims[2], float *float_se) {
  tt::max_filter_square(foutput, dem, tmp, se_width, dims);

  tt::morphological_dilation(soutput, dem, float_se, dims, se_dims);

  for (ptrdiff_t row = 0; row < dims[1]; row++) {
    for (ptrdiff_t col = 0; col < dims[0]; col++) {
      ptrdiff_t index = col + row * dims[0];

      if ((std::isnan(foutput[index]) != std::isnan(soutput[index])) ||
          // second isnan check not needed but added for completeness
          (!std::isnan(foutput[index]) && !std::isnan(soutput[index]) &&
           foutput[index] != soutput[index])) {
        std::cout << "(" << col << ", " << row << "): " << foutput[index]
                  << " != " << soutput[index] << std::endl;
        assert(0);
      }
    }
  }

  return 0;
}

int test_on_random_dem(uint32_t seed) {
  // sizes between 1 and 513
  ptrdiff_t dimensions[2] = {rand() % 512 + 1, rand() % 512 + 1};

  float *dem = new float[dimensions[0] * dimensions[1]];
  float *tmp = new float[dimensions[0] * dimensions[1]];
  float *output = new float[dimensions[0] * dimensions[1]];
  float *second_output = new float[dimensions[0] * dimensions[1]];

  // Initialize (pseudo) random DEM
  for (uint32_t col = 0; col < dimensions[1]; col++) {
    for (uint32_t row = 0; row < dimensions[0]; row++) {
      if (rand() % 100 < 5) {
        dem[col * dimensions[0] + row] = NAN;
      } else {
        dem[col * dimensions[0] + row] =
            100.0f * utils::pcg4d(row, col, seed, 1);
      }
    }
  }

  ptrdiff_t se_size[3] = {3, 3, 1};
  uint8_t binary_se[9] = {1, 1, 1, 1, 1, 1, 1, 1, 1};
  float float_se_identical_to_binary[9] = {0.0f, 0.0f, 0.0f, 0.0f, 0.0f,
                                           0.0f, 0.0f, 0.0f, 0.0f};
  uint8_t binary_se_triangle[9] = {1, 0, 0, 1, 1, 0, 1, 1, 1};
  float float_se_triangle[9] = {0.0f, NAN,  NAN,  0.0f, 0.0f,
                                NAN,  0.0f, 0.0f, 0.0f};

  test_no_larger_values(dem, output, binary_se, dimensions, se_size);

  test_no_larger_values_square(dem, output, tmp, std::size(se_size),
                               dimensions);

  test_no_smaller_values(dem, output, binary_se, dimensions, se_size);

  test_no_smaller_values_square(dem, output, tmp, std::size(se_size),
                                dimensions);

  test_nan_values_persist_min_filter(dem, output, binary_se, dimensions,
                                     se_size);

  test_nan_values_persist_min_filter_square(dem, output, tmp,
                                            std::size(se_size), dimensions);

  test_nan_values_persist_max_filter(dem, output, binary_se, dimensions,
                                     se_size);

  test_nan_values_persist_max_filter_square(dem, output, tmp,
                                            std::size(se_size), dimensions);

  test_min_filter_implementations_agree(dem, output, second_output, tmp,
                                        std::size(se_size), binary_se,
                                        dimensions, se_size);

  test_max_filter_implementations_agree(dem, output, second_output, tmp,
                                        std::size(se_size), binary_se,
                                        dimensions, se_size);

  test_min_filter_and_erosion_agree(dem, output, second_output, binary_se,
                                    dimensions, se_size,
                                    float_se_identical_to_binary);

  test_min_filter_and_erosion_agree(dem, output, second_output,
                                    binary_se_triangle, dimensions, se_size,
                                    float_se_triangle);

  test_min_square_filter_and_erosion_agree(
      dem, output, second_output, tmp, std::size(se_size), dimensions, se_size,
      float_se_identical_to_binary);

  test_max_filter_and_dilation_agree(dem, output, second_output, binary_se,
                                     dimensions, se_size,
                                     float_se_identical_to_binary);

  test_max_filter_and_dilation_agree(dem, output, second_output,
                                     binary_se_triangle, dimensions, se_size,
                                     float_se_triangle);

  test_max_square_filter_and_dilation_agree(
      dem, output, second_output, tmp, std::size(se_size), dimensions, se_size,
      float_se_identical_to_binary);

  delete[] dem;
  delete[] tmp;
  delete[] output;
  delete[] second_output;

  return 0;
}

int main(int argc, char *argv[]) {
  for (uint32_t test = 0; test < 100; test++) {
    srand(test);
    test_on_random_dem(test);
  }

  return 0;
}
