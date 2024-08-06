#undef NDEBUG
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iostream>

extern "C" {
#include "topotoolbox.h"
#include "utils.h"
}

#define SQRT2f 1.41421356237309504880f

/*
  Each pixel of the filled DEM should be greater than or equal to
  the corresponding pixel in the original DEM.
 */
int32_t test_fillsinks_ge(float *original_dem, float *filled_dem,
                          ptrdiff_t dims[2]) {
  for (ptrdiff_t j = 0; j < dims[1]; j++) {
    for (ptrdiff_t i = 0; i < dims[0]; i++) {
      assert(filled_dem[i + dims[0] * j] >= original_dem[i + dims[0] * j]);
    }
  }
  return 0;
}

/*
  No pixel in the filled DEM should be completely surrounded by pixels higher
  than it.
 */
int32_t test_fillsinks_filled(float *filled_dem, ptrdiff_t dims[2]) {
  ptrdiff_t j_offset[8] = {-1, -1, -1, 0, 1, 1, 1, 0};
  ptrdiff_t i_offset[8] = {1, 0, -1, -1, -1, 0, 1, 1};

  for (ptrdiff_t j = 0; j < dims[1]; j++) {
    for (ptrdiff_t i = 0; i < dims[0]; i++) {
      float z = filled_dem[i + dims[0] * j];
      ptrdiff_t up_neighbor_count = 0;
      for (int32_t neighbor = 0; neighbor < 8; neighbor++) {
        ptrdiff_t neighbor_i = i + i_offset[neighbor];
        ptrdiff_t neighbor_j = j + j_offset[neighbor];

        if (neighbor_i < 0 || neighbor_i >= dims[0] || neighbor_j < 0 ||
            neighbor_j >= dims[1]) {
          continue;
        }

        if (filled_dem[neighbor_i + dims[0] * neighbor_j] > z) {
          up_neighbor_count++;
        }
      }
      assert(up_neighbor_count != 8);
    }
  }
  return 0;
}

/*
  Every pixel not on the boundary with no lower neighbors and fewer than
  8 higher neighbors should be labeled a flat. Likewise, every flat
  should have no lower neighbors and fewer than 8 higher neighbors.
 */
int32_t test_identifyflats_flats(int32_t *flats, float *dem,
                                 ptrdiff_t dims[2]) {
  ptrdiff_t col_offset[8] = {-1, -1, -1, 0, 1, 1, 1, 0};
  ptrdiff_t row_offset[8] = {1, 0, -1, -1, -1, 0, 1, 1};

  for (ptrdiff_t j = 0; j < dims[1]; j++) {
    for (ptrdiff_t i = 0; i < dims[0]; i++) {
      float z = dem[i + dims[0] * j];
      int32_t flat = flats[i + dims[0] * j];

      int32_t up_neighbor_count = 0;
      int32_t down_neighbor_count = 0;
      for (int32_t neighbor = 0; neighbor < 8; neighbor++) {
        ptrdiff_t neighbor_i = i + row_offset[neighbor];
        ptrdiff_t neighbor_j = j + col_offset[neighbor];

        if (neighbor_i < 0 || neighbor_i >= dims[0] || neighbor_j < 0 ||
            neighbor_j >= dims[1]) {
          continue;
        }

        float neighbor_height = dem[neighbor_i + dims[0] * neighbor_j];

        if (neighbor_height > z) {
          up_neighbor_count++;
        } else if (neighbor_height < z) {
          down_neighbor_count++;
        }
      }
      int32_t current_pixel_on_boundary =
          i == 0 || i == dims[0] - 1 || j == 0 || j == dims[1] - 1;

      assert((!current_pixel_on_boundary && (down_neighbor_count == 0) &&
              (up_neighbor_count < 8)) == ((flat & 1) == 1));
    }
  }

  return 0;
}

/*
  Every pixel that has a lower neighbor, borders a flat, and has the
  same elevation as a flat that it touches should be labeled a sill.
*/
int32_t test_identifyflats_sills(int32_t *flats, float *dem,
                                 ptrdiff_t dims[2]) {
  ptrdiff_t j_offset[8] = {-1, -1, -1, 0, 1, 1, 1, 0};
  ptrdiff_t i_offset[8] = {1, 0, -1, -1, -1, 0, 1, 1};

  for (ptrdiff_t j = 0; j < dims[1]; j++) {
    for (ptrdiff_t i = 0; i < dims[0]; i++) {
      float z = dem[i + j * dims[0]];
      int32_t flat = flats[i + j * dims[0]];

      int32_t down_neighbor_count = 0;
      int32_t equal_neighbor_flats = 0;

      for (int32_t neighbor = 0; neighbor < 8; neighbor++) {
        ptrdiff_t neighbor_i = i + i_offset[neighbor];
        ptrdiff_t neighbor_j = j + j_offset[neighbor];

        if (neighbor_i < 0 || neighbor_i >= dims[0] || neighbor_j < 0 ||
            neighbor_j >= dims[1]) {
          // Count boundary pixels as down neighbors so sills can
          // drain off the map.
          down_neighbor_count++;
          continue;
        }

        float neighbor_height = dem[neighbor_i + dims[0] * neighbor_j];
        int32_t neighbor_flat = flats[neighbor_i + dims[0] * neighbor_j];

        if (neighbor_height < z) {
          down_neighbor_count++;
        }

        if ((neighbor_flat & 1) && (neighbor_height == z)) {
          equal_neighbor_flats++;
        }
      }
      assert(((flat & 2) > 0) ==
             ((down_neighbor_count > 0) && (equal_neighbor_flats > 0)));
    }
  }

  return 0;
}

/*
  Every pixel that is a flat and borders a sill of the same
  height should be labeled a presill
*/
int32_t test_identifyflats_presills(int32_t *flats, float *dem,
                                    ptrdiff_t dims[2]) {
  ptrdiff_t j_offset[8] = {-1, -1, -1, 0, 1, 1, 1, 0};
  ptrdiff_t i_offset[8] = {1, 0, -1, -1, -1, 0, 1, 1};

  for (ptrdiff_t j = 0; j < dims[1]; j++) {
    for (ptrdiff_t i = 0; i < dims[0]; i++) {
      float z = dem[i + j * dims[0]];
      int32_t flat = flats[i + j * dims[0]];

      int32_t equal_neighbor_sills = 0;

      for (int32_t neighbor = 0; neighbor < 8; neighbor++) {
        ptrdiff_t neighbor_i = i + i_offset[neighbor];
        ptrdiff_t neighbor_j = j + j_offset[neighbor];

        if (neighbor_i < 0 || neighbor_i >= dims[0] || neighbor_j < 0 ||
            neighbor_j >= dims[1]) {
          continue;
        }

        float neighbor_height = dem[neighbor_i + dims[0] * neighbor_j];
        int32_t neighbor_flat = flats[neighbor_i + dims[0] * neighbor_j];

        if (((neighbor_flat & 2) > 0) && (neighbor_height == z)) {
          equal_neighbor_sills++;
        }
      }
      assert(((flat & 1) > 0 && (equal_neighbor_sills > 0)) ==
             ((flat & 4) > 0));
    }
  }

  return 0;
}

/*
  Costs should be zero on nonflats and positive on flats.
 */
int32_t test_gwdt_costs(float *costs, int32_t *flats, ptrdiff_t dims[2]) {
  for (ptrdiff_t j = 0; j < dims[1]; j++) {
    for (ptrdiff_t i = 0; i < dims[0]; i++) {
      float cost = costs[i + dims[0] * j];
      int32_t flat = flats[i + dims[0] * j];

      assert((cost >= 0) && (((flat & 1) > 0) == (cost > 0)));
    }
  }

  return 0;
}

/*
  Flats neighboring other flats should all have the same nonzero
  connected component label
 */
int32_t test_gwdt_conncomps(ptrdiff_t *conncomps, int32_t *flats,
                            ptrdiff_t dims[2]) {
  ptrdiff_t j_offset[8] = {-1, -1, -1, 0, 1, 1, 1, 0};
  ptrdiff_t i_offset[8] = {1, 0, -1, -1, -1, 0, 1, 1};

  for (ptrdiff_t j = 0; j < dims[1]; j++) {
    for (ptrdiff_t i = 0; i < dims[0]; i++) {
      ptrdiff_t label = conncomps[i + dims[0] * j];
      int32_t flat = flats[i + dims[0] * j];

      if (!(flat & 1)) {
        continue;
      }

      for (int32_t neighbor = 0; neighbor < 8; neighbor++) {
        ptrdiff_t neighbor_i = i + i_offset[neighbor];
        ptrdiff_t neighbor_j = j + j_offset[neighbor];

        if (neighbor_i < 0 || neighbor_i >= dims[0] || neighbor_j < 0 ||
            neighbor_j >= dims[1]) {
          continue;
        }

        int32_t neighbor_flat = flats[neighbor_i + dims[0] * neighbor_j];
        ptrdiff_t neighbor_label = conncomps[neighbor_i + dims[0] * neighbor_j];
        if (neighbor_flat & 1) {
          assert(neighbor_label == label);
        }
      }
    }
  }
  return 0;
}

/*
  Gray-weighted distances should be

  1. non-negative for all pixels
  2. 1 for presills
  3. >1 for other flats
  4. equal to the distance to the parent pixel plus the geodesic distance
  between the parent and the current pixel.
 */
int32_t test_gwdt(float *dist, ptrdiff_t *prev, float *costs, int32_t *flats,
                  ptrdiff_t dims[2]) {
  for (ptrdiff_t j = 0; j < dims[1]; j++) {
    for (ptrdiff_t i = 0; i < dims[0]; i++) {
      float d = dist[i + j * dims[0]];
      int32_t flat = flats[i + j * dims[0]];

      assert(d >= 0);

      if (flat & 4) {
        assert(d == 1.0);
      } else if (flat & 1) {
        assert(d > 1.0);

        ptrdiff_t parent = prev[i + j * dims[0]];
        ptrdiff_t parent_j = parent / dims[0];
        ptrdiff_t parent_i = parent % dims[0];
        float chamfer = (parent_j != j) && (parent_i != i) ? SQRT2f : 1.0f;
        float proposed_dist =
            dist[parent] +
            chamfer * (costs[i + j * dims[0]] + costs[parent]) / 2;
        assert(proposed_dist == d);
      }
    }
  }

  return 0;
}

/*
  Flow direction should point downstream or across flats
 */
int32_t test_routeflowd8_direction(uint8_t *direction, float *filled_dem,
                                   float *dist, int32_t *flats,
                                   ptrdiff_t dims[2]) {
  // 1<<5 1<<6  1<<7
  // 1<<4    0  1<<0
  // 1<<3 1<<2  1<<1
  ptrdiff_t i_offset[8] = {0, 1, 1, 1, 0, -1, -1, -1};
  ptrdiff_t j_offset[8] = {1, 1, 0, -1, -1, -1, 0, 1};

  for (ptrdiff_t j = 0; j < dims[1]; j++) {
    for (ptrdiff_t i = 0; i < dims[0]; i++) {
      uint8_t flowdir = direction[j * dims[0] + i];

      if (flowdir == 0) {
        // Pixel is a sink
        continue;
      }

      float z = filled_dem[j * dims[0] + i];
      uint8_t v = flowdir;
      uint8_t r = 0;
      while (v >>= 1) {
        r++;
      }

      ptrdiff_t neighbor_i = i + i_offset[r];
      ptrdiff_t neighbor_j = j + j_offset[r];

      // Check that the downstream neighbor is within the array bounds
      assert(0 <= neighbor_i);
      assert(neighbor_i < dims[0]);
      assert(0 <= neighbor_j);
      assert(neighbor_j < dims[1]);

      // Neighbor elevation should be less than or equal to the
      // current elevation
      assert(filled_dem[neighbor_j * dims[0] + neighbor_i] <= z);
    }
  }
  return 0;
}

/*
  source should be a topological sort of the
  graph defined by direction.

  This is O(N^2) and very slow.
 */
int32_t test_routeflowd8_tsort(uint8_t *marks, ptrdiff_t *source,
                               uint8_t *direction, ptrdiff_t dims[2]) {
  // 1<<5 1<<6  1<<7
  // 1<<4    0  1<<0
  // 1<<3 1<<2  1<<1
  ptrdiff_t i_offset[8] = {0, 1, 1, 1, 0, -1, -1, -1};
  ptrdiff_t j_offset[8] = {1, 1, 0, -1, -1, -1, 0, 1};

  for (ptrdiff_t j = 0; j < dims[1]; j++) {
    for (ptrdiff_t i = 0; i < dims[0]; i++) {
      ptrdiff_t src = source[j * dims[0] + i];
      ptrdiff_t src_i = src % dims[0];
      ptrdiff_t src_j = src / dims[0];

      uint8_t flowdir = direction[src];

      if (flowdir == 0) {
        // src is a sink
        continue;
      }

      uint8_t v = flowdir;
      uint8_t r = 0;
      while (v >>= 1) {
        r++;
      }
      ptrdiff_t neighbor_i = src_i + i_offset[r];
      ptrdiff_t neighbor_j = src_j + j_offset[r];

      ptrdiff_t neighbor_src = neighbor_j * dims[0] + neighbor_i;
      // Check that we haven't seen neighbor in the topological sort yet
      assert(marks[neighbor_src] != 0xff);

      marks[src] = 0xff;
    }
  }
  return 0;
}

int32_t test_flowaccumulation_max(float *acc, ptrdiff_t dims[2]) {
  for (ptrdiff_t j = 0; j < dims[1]; j++) {
    for (ptrdiff_t i = 0; i < dims[0]; i++) {
      assert(acc[j * dims[0] + i] < (dims[0] * dims[1]));
    }
  }
  return 0;
}

int32_t random_dem_test(ptrdiff_t dims[2], uint32_t seed) {
  // Allocate variables

  // Input DEM
  float *dem = new float[dims[0] * dims[1]];
  assert(dem != NULL);

  // Output for fillsinks
  float *filled_dem = new float[dims[0] * dims[1]];
  assert(filled_dem != NULL);

  // Output for identifyflats
  int32_t *flats = new int32_t[dims[0] * dims[1]];
  assert(flats != NULL);

  // Outputs for compute_costs
  ptrdiff_t *conncomps = new ptrdiff_t[dims[0] * dims[1]];
  assert(conncomps != NULL);
  float *costs = new float[dims[0] * dims[1]];
  assert(costs != NULL);

  // Outputs and intermediate needs for gwdt
  float *dist = new float[dims[0] * dims[1]];
  assert(dist != NULL);
  ptrdiff_t *heap = new ptrdiff_t[dims[0] * dims[1]];
  assert(heap != NULL);
  ptrdiff_t *back = new ptrdiff_t[dims[0] * dims[1]];
  assert(back != NULL);
  ptrdiff_t *prev = new ptrdiff_t[dims[0] * dims[1]];
  assert(prev != NULL);

  // Outputs for routeflowd8
  ptrdiff_t *source = new ptrdiff_t[dims[0] * dims[1]];
  assert(source != NULL);
  uint8_t *direction = new uint8_t[dims[0] * dims[1]];
  assert(direction != NULL);
  // marks is only used for testing
  uint8_t *marks =
      new uint8_t[dims[0] * dims[1]]();  // Initialize marks to zero
  assert(marks != NULL);

  float *accum = new float[dims[0] * dims[1]]();
  assert(accum != NULL);

  for (uint32_t col = 0; col < dims[1]; col++) {
    for (uint32_t row = 0; row < dims[0]; row++) {
      dem[col * dims[0] + row] = 100.0f * pcg4d(row, col, seed, 1);
    }
  }

  // Run flow routing algorithms
  fillsinks(filled_dem, dem, dims);

  test_fillsinks_ge(dem, filled_dem, dims);
  test_fillsinks_filled(filled_dem, dims);

  identifyflats(flats, filled_dem, dims);

  test_identifyflats_flats(flats, filled_dem, dims);
  test_identifyflats_sills(flats, filled_dem, dims);
  test_identifyflats_presills(flats, filled_dem, dims);

  gwdt_computecosts(costs, conncomps, flats, dem, filled_dem, dims);

  test_gwdt_costs(costs, flats, dims);
  test_gwdt_conncomps(conncomps, flats, dims);

  gwdt(dist, prev, costs, flats, heap, back, dims);
  test_gwdt(dist, prev, costs, flats, dims);

  flow_routing_d8_carve(source, direction, filled_dem, dist, flats, dims);
  test_routeflowd8_direction(direction, filled_dem, dist, flats, dims);
  test_routeflowd8_tsort(marks, source, direction, dims);

  flow_accumulation(accum, source, direction, NULL, dims);
  test_flowaccumulation_max(accum, dims);

  delete[] dem;
  delete[] filled_dem;
  delete[] flats;
  delete[] conncomps;
  delete[] costs;
  delete[] dist;
  delete[] heap;
  delete[] back;
  delete[] prev;
  delete[] source;
  delete[] direction;
  delete[] marks;
  delete[] accum;

  return 0;
}

int main(int argc, char *argv[]) {
  ptrdiff_t dims[2] = {100, 200};

  for (uint32_t test = 0; test < 100; test++) {
    int32_t result = random_dem_test(dims, test);
    if (result < 0) {
      return result;
    }
  }
}
