#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iostream>

extern "C" {
#include "topotoolbox.h"
#include "utils.h"
}

int32_t random_dem_test(ptrdiff_t nrows, ptrdiff_t ncols, uint32_t seed) {
  // Allocate variables

  // Input DEM
  float *dem = new float[nrows * ncols];

  // Output for fillsinks
  float *filled_dem = new float[nrows * ncols];

  // Output for identifyflats
  int32_t *flats = new int32_t[nrows * ncols];

  // Outputs for compute_costs
  ptrdiff_t *conncomps = new ptrdiff_t[nrows * ncols];
  float *costs = new float[nrows * ncols];

  for (uint32_t col = 0; col < ncols; col++) {
    for (uint32_t row = 0; row < nrows; row++) {
      dem[col * nrows + row] = 100.0f * pcg4d(row, col, seed, 1);
    }
  }

  // Run flow routing algorithms
  fillsinks(filled_dem, dem, nrows, ncols);
  ptrdiff_t count_flats = identifyflats(flats, filled_dem, nrows, ncols);
  gwdt_computecosts(costs, conncomps, flats, dem, filled_dem, nrows, ncols);

  // Number of flats identified in the test
  ptrdiff_t test_count_flats = 0;

  // Test properties of filled DEM and identified flats
  ptrdiff_t col_offset[8] = {-1, -1, -1, 0, 1, 1, 1, 0};
  ptrdiff_t row_offset[8] = {1, 0, -1, -1, -1, 0, 1, 1};

  for (ptrdiff_t col = 0; col < ncols; col++) {
    for (ptrdiff_t row = 0; row < nrows; row++) {
      float z = filled_dem[col * nrows + row];
      int32_t flat = flats[col * nrows + row];
      float cost = costs[col * nrows + row];
      ptrdiff_t label = conncomps[col * nrows + row];

      int32_t current_pixel_on_border =
          row == 0 || row == nrows - 1 || col == 0 || col == ncols - 1;

      // Each pixel of the filled raster should be >= the DEM
      if (z < dem[col * nrows + row]) {
        std::cout << "Pixel (" << row << ", " << col << ") is below the DEM"
                  << std::endl;
        std::cout << "Value: " << z << std::endl;
        std::cout << "DEM: " << dem[col * nrows + row] << std::endl;
        return -1;
      }

      // Test cost computation
      if (flat & 1) {
        if (cost <= 0) {
          std::cout << "The cost at pixel (" << row << ", " << col
                    << ") is nonpositive" << std::endl;
          return -1;
        }

        if (label == 0) {
          std::cout << "Pixel (" << row << ", " << col
                    << ") is a flat but has a connected component label of 0"
                    << std::endl;
          return -1;
        }
      }

      if (!(flat & 1)) {
        if (cost != 0) {
          std::cout << "The cost " << cost << " at nonflat pixel (" << row
                    << ", " << col << ") is nonzero" << std::endl;
          return -1;
        }
      }

      // Number of neighbors lower than the current pixel
      int32_t down_neighbor_count = 0;

      // Number of neighbors higher than the current pixel
      int32_t up_neighbor_count = 0;

      // Number of neighbors that are flats
      int32_t neighboring_flats = 0;

      // Number of neighboring flats that have the same elevation as
      // the current pixel
      int32_t equal_neighboring_flats = 0;

      // Number of neighboring sills that have the same elevation as
      // the current pixel
      int32_t equal_neighboring_sills = 0;

      for (ptrdiff_t neighbor = 0; neighbor < 8; neighbor++) {
        ptrdiff_t neighbor_row = row + row_offset[neighbor];
        ptrdiff_t neighbor_col = col + col_offset[neighbor];

        if (neighbor_row < 0 || neighbor_row >= nrows || neighbor_col < 0 ||
            neighbor_col >= ncols) {
          continue;
        }

        float neighbor_height = filled_dem[neighbor_col * nrows + neighbor_row];
        int32_t neighboring_flat = flats[neighbor_col * nrows + neighbor_row];
        ptrdiff_t neighbor_label =
            conncomps[neighbor_col * nrows + neighbor_row];

        if (z < neighbor_height) {
          up_neighbor_count++;
        }

        if (z > neighbor_height) {
          down_neighbor_count++;
        }

        if (neighboring_flat & 1) {
          neighboring_flats++;

          if (z == neighbor_height) {
            equal_neighboring_flats++;
          }

          // Flats neighboring other flats should all have the same
          // connected component label
          if ((flat & 1) && (neighbor_label != label)) {
            std::cout << "Pixel (" << row << ", " << col
                      << ") is a flat bordering another flat with a different "
                         "connected component label"
                      << std::endl;
            std::cout << "Label: " << label << std::endl;
            std::cout << "Neighbor label: " << neighbor_label << std::endl;
            return -1;
          }
        }

        if (neighboring_flat & 2) {
          if (z == neighbor_height) {
            equal_neighboring_sills++;
          }
        }
      }

      // No pixel of the filled raster should be surrounded by
      // neighbors that are higher than it
      if (up_neighbor_count == 8) {
        std::cout << "Pixel (" << row << ", " << col << ") is a sink"
                  << std::endl;
        return -1;
      }

      if (!current_pixel_on_border && up_neighbor_count < 8 &&
          down_neighbor_count == 0) {
        // This pixel is a flat
        test_count_flats++;
      }

      // Every pixel with no lower neighbors and fewer than 8
      // up_neighbors should be labeled a flat
      if (!current_pixel_on_border && up_neighbor_count < 8 &&
          down_neighbor_count == 0 && !(flat & 1)) {
        std::cout << "Pixel (" << row << ", " << col
                  << ") is has no lower neighbors but is not labeled a flat"
                  << std::endl;
        return -1;
      }

      // Every pixel that has a lower neighbor, borders a flat, and
      // has the same elevation as a flat that it touches should be
      // labeled a sill.
      if (equal_neighboring_flats > 0 && down_neighbor_count > 0 &&
          !(flat & 2)) {
        std::cout << "Pixel (" << row << ", " << col
                  << ") neighbors a flat and has lower neighbors but is not "
                     "labeled a sill"
                  << std::endl;
        return -1;
      }

      // Every pixel that is a flat and borders a sill of the same
      // height should be labeled a presill
      if ((flat & 1) && (equal_neighboring_sills > 0) && !(flat & 4)) {
        std::cout
            << "Pixel (" << row << ", " << col
            << ") is a flat that neighbors a sill but is not labeled a presill"
            << std::endl;
        return -1;
      }

      // No flat should border a lower pixel
      if ((flat & 1) && (down_neighbor_count > 0)) {
        std::cout << "Pixel (" << row << ", " << col
                  << ") is a flat but has a lower neighbor" << std::endl;
        return -1;
      }

      // Every sill pixel should have at least one neighbor lower than it,
      // unless it is on a border
      if (!current_pixel_on_border && (flat & 2) &&
          (down_neighbor_count == 0)) {
        std::cout << "Pixel (" << row << ", " << col
                  << ") is a sill but has no lower neighbor" << std::endl;
        return -1;
      }

      // Every sill pixel should border a flat
      if ((flat & 2) && (neighboring_flats == 0)) {
        std::cout << "Pixel (" << row << ", " << col
                  << ") is a sill but does not border a flat" << std::endl;
        return -1;
      }
    }
  }

  if (test_count_flats != count_flats) {
    std::cout << "Number of flats identified in the test is not equal to the "
                 "number of flats computed by identifyflats"
              << std::endl;
    return -1;
  }

  delete[] dem;
  delete[] filled_dem;
  delete[] flats;
  delete[] conncomps;
  delete[] costs;

  return 0;
}

int main(int argc, char *argv[]) {
  ptrdiff_t nrows = 100;
  ptrdiff_t ncols = 200;

  for (uint32_t test = 0; test < 100; test++) {
    int32_t result = random_dem_test(nrows, ncols, test);
    if (result < 0) {
      return result;
    }
  }
}
