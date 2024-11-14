#undef NDEBUG
#include <array>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <ctime>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

// Include topotoolbox.h in its own namespace to help prevent naming
// conflicts in the global scope.
namespace tt {
extern "C" {
#include "topotoolbox.h"
}
}  // namespace tt

extern "C" {
#include "utils.h"
}

#include "profiler.h"
Profiler prof;  // Global Profiler that our test functions can access

#define SQRT2f 1.41421356237309504880f
#define SQRT2 1.41421356237309504880f

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
  For each cell, the saved gradient should be the steepest gradient of all 8
  neighboring cells.
*/
int32_t test_gradient8(float *gradient, float *dem, float cellsize,
                       ptrdiff_t dims[2]) {
  ptrdiff_t i_offset[8] = {0, 1, 1, 1, 0, -1, -1, -1};
  ptrdiff_t j_offset[8] = {1, 1, 0, -1, -1, -1, 0, 1};
  for (ptrdiff_t j = 0; j < dims[1]; j++) {
    for (ptrdiff_t i = 0; i < dims[0]; i++) {
      float max_gradient = 0;

      // Iterate over all 8 neighboring cells
      for (int k = 0; k < 8; k++) {
        ptrdiff_t neighbor_i = i + i_offset[k];
        ptrdiff_t neighbor_j = j + j_offset[k];

        // Check if cells in bounds
        if (neighbor_i >= 0 && neighbor_i < dims[0] && neighbor_j >= 0 &&
            neighbor_j < dims[1]) {
          float horizontal_dist;
          float vertical_dist;
          float local_gradient;

          if (neighbor_i != i && neighbor_j != j) {
            horizontal_dist = SQRT2 * cellsize;
          } else {
            horizontal_dist = cellsize;
          }
          vertical_dist = fabsf(dem[neighbor_j * dims[0] + neighbor_i] -
                                dem[j * dims[0] + i]);

          local_gradient = vertical_dist / horizontal_dist;
          if (local_gradient > max_gradient) {
            max_gradient = local_gradient;
          }
        }
      }
      assert(max_gradient == gradient[j * dims[0] + i]);
    }
  }
  return 0;
}
/*
  Computing the gradient8 using OpenMP should yield the same result as without.
*/
int32_t test_gradient8_mp(float *gradient, float *gradient_mp,
                          ptrdiff_t dims[2]) {
  for (ptrdiff_t j = 0; j < dims[1]; j++) {
    for (ptrdiff_t i = 0; i < dims[0]; i++) {
      ptrdiff_t position = j * dims[0] + i;
      assert(gradient[position] == gradient_mp[position]);
    }
  }
  return 0;
}
/*

*/
int32_t test_acv() { return 0; }
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

int32_t test_flow_accumulation_max(float *acc, ptrdiff_t dims[2]) {
  for (ptrdiff_t j = 0; j < dims[1]; j++) {
    for (ptrdiff_t i = 0; i < dims[0]; i++) {
      assert(acc[j * dims[0] + i] < (dims[0] * dims[1]));
    }
  }
  return 0;
}

int32_t test_flow_accumulation_multimethod(float *acc1, float *acc2,
                                           ptrdiff_t dims[2]) {
  for (ptrdiff_t j = 0; j < dims[1]; j++) {
    for (ptrdiff_t i = 0; i < dims[0]; i++) {
      assert(acc1[j * dims[0] + i] == acc2[j * dims[0] + i]);
    }
  }
  return 0;
}

/*
  The target corresponding to source should be equal to the source
  plus the correct offset in the flow direction.
 */
int32_t test_flow_routing_targets(ptrdiff_t *target, ptrdiff_t *source,
                                  uint8_t *direction, ptrdiff_t dims[2]) {
  ptrdiff_t i_offset[8] = {0, 1, 1, 1, 0, -1, -1, -1};
  ptrdiff_t j_offset[8] = {1, 1, 0, -1, -1, -1, 0, 1};

  for (ptrdiff_t j = 0; j < dims[1]; j++) {
    for (ptrdiff_t i = 0; i < dims[0]; i++) {
      ptrdiff_t u = source[j * dims[0] + i];
      ptrdiff_t u_i = u % dims[0];
      ptrdiff_t u_j = u / dims[0];

      assert(u_i >= 0 && u_i < dims[0]);
      assert(u_j >= 0 && u_j < dims[1]);

      uint8_t flowdir = direction[u];

      if (flowdir == 0) {
        assert(target[j * dims[0] + i] == -1);
      } else {
        uint8_t r = 0;
        while (flowdir >>= 1) {
          r++;
        }

        ptrdiff_t neighbor_i = u_i + i_offset[r];
        ptrdiff_t neighbor_j = u_j + j_offset[r];

        assert(neighbor_i >= 0 && neighbor_i < dims[0]);
        assert(neighbor_j >= 0 && neighbor_j < dims[1]);

        ptrdiff_t v = neighbor_j * dims[0] + neighbor_i;

        assert(v == target[j * dims[0] + i]);
      }
    }
  }
  return 0;
}

/*
  Compute the upstream distance by computing shortest paths in the
  flow network and compare them to the distances computed by streamquad_trapz
  and recorded in node_distance.
 */
int32_t test_stream_distance(float *node_distance, ptrdiff_t *stream_grid,
                             float *distance, ptrdiff_t *source,
                             ptrdiff_t *target, float cellsize,
                             ptrdiff_t dims[2]) {
  for (ptrdiff_t edge = dims[0] * dims[1] - 1; edge >= 0; edge--) {
    ptrdiff_t u = source[edge];
    ptrdiff_t v = target[edge];
    if (v >= 0) {
      float w = llabs(u - v) == 1 || llabs(u - v) == dims[0] ? 1.0f : SQRT2f;
      w *= cellsize;

      distance[u] = distance[v] + w;
    }
  }

  for (ptrdiff_t j = 0; j < dims[1]; j++) {
    for (ptrdiff_t i = 0; i < dims[0]; i++) {
      ptrdiff_t node = stream_grid[j * dims[0] + i];
      if (node > 0) {
        assert(node_distance[node - 1] == distance[j * dims[0] + i]);
      }
    }
  }
  return 0;
}

struct FlowRoutingData {
  std::array<ptrdiff_t, 2> dims;
  float cellsize;

  // input data
  std::vector<float> dem;
  std::vector<uint8_t> bc;
  std::vector<float> fraction;

  // fillsinks
  std::vector<float> filled_dem;
  // fillsinks_hybrid
  std::vector<ptrdiff_t> queue;

  // identifyflats
  std::vector<int32_t> flats;

  // gwdt_compute_costs
  std::vector<ptrdiff_t> conncomps;
  std::vector<float> costs;

  // gwdt
  std::vector<float> dist;
  std::vector<ptrdiff_t> heap;
  std::vector<ptrdiff_t> back;
  std::vector<ptrdiff_t> prev;

  // gradient8
  std::vector<float> gradient;
  std::vector<float> gradient_mp;

  // flow_routing_d8_carve
  std::vector<ptrdiff_t> source;
  std::vector<uint8_t> direction;
  std::vector<uint8_t> marks;

  // flow_routing_targets
  std::vector<ptrdiff_t> target;

  // flow_accumulation
  std::vector<float> accum;
  std::vector<float> accum2;

  // stream network
  std::vector<std::tuple<ptrdiff_t, ptrdiff_t, float>> stream_edges;
  std::vector<ptrdiff_t> stream_source;
  std::vector<ptrdiff_t> stream_target;
  std::vector<float> stream_weight;
  std::vector<ptrdiff_t> stream_grid;
  ptrdiff_t stream_node_count;

  FlowRoutingData(ptrdiff_t input_dims[2], float cs, uint32_t seed)
      : dims({input_dims[0], input_dims[1]}),
        cellsize(cs),
        dem(dims[0] * dims[1]),
        bc(dims[0] * dims[1]),
        fraction(dims[0] * dims[1]),
        filled_dem(dims[0] * dims[1]),
        queue(dims[0] * dims[1]),
        flats(dims[0] * dims[1]),
        conncomps(dims[0] * dims[1]),
        costs(dims[0] * dims[1]),
        dist(dims[0] * dims[1]),
        heap(dims[0] * dims[1]),
        back(dims[0] * dims[1]),
        prev(dims[0] * dims[1]),
        gradient(dims[0] * dims[1]),
        gradient_mp(dims[0] * dims[1]),
        source(dims[0] * dims[1]),
        direction(dims[0] * dims[1]),
        marks(dims[0] * dims[1]),
        target(dims[0] * dims[1]),
        accum(dims[0] * dims[1]),
        accum2(dims[0] * dims[1]),
        stream_grid(dims[0] * dims[1]) {
    // Initialize DEM, boundary conditions for fillsinks and fraction
    for (uint32_t col = 0; col < dims[1]; col++) {
      for (uint32_t row = 0; row < dims[0]; row++) {
        dem[col * dims[0] + row] = 100.0f * pcg4d(row, col, seed, 1);

        if (row == 0 || row == dims[0] - 1 || col == 0 || col == dims[1] - 1) {
          bc[col * dims[0] + row] = 1;
        } else {
          bc[col * dims[0] + row] = 0;
        }

        fraction[col * dims[0] + row] = 1.0f;
      }
    }
  }

  void route_flow() {
    ProfileFunction(prof);

    tt::fillsinks(filled_dem.data(), dem.data(), bc.data(), dims.data());

    tt::identifyflats(flats.data(), filled_dem.data(), dims.data());

    tt::gwdt_computecosts(costs.data(), conncomps.data(), flats.data(),
                          dem.data(), filled_dem.data(), dims.data());

    tt::gwdt(dist.data(), prev.data(), costs.data(), flats.data(), heap.data(),
             back.data(), dims.data());

    tt::flow_routing_d8_carve(source.data(), direction.data(),
                              filled_dem.data(), dist.data(), flats.data(),
                              dims.data());

    tt::flow_routing_targets(target.data(), source.data(), direction.data(),
                             dims.data());

    tt::flow_accumulation_edgelist(accum2.data(), source.data(), target.data(),
                                   fraction.data(), NULL, dims[0] * dims[1],
                                   dims.data());

    // Close timer
  }

  void route_flow_hybrid() {
    ProfileFunction(prof);

    tt::fillsinks_hybrid(filled_dem.data(), queue.data(), dem.data(), bc.data(),
                         dims.data());

    tt::identifyflats(flats.data(), filled_dem.data(), dims.data());

    tt::gwdt_computecosts(costs.data(), conncomps.data(), flats.data(),
                          dem.data(), filled_dem.data(), dims.data());

    tt::gwdt(dist.data(), prev.data(), costs.data(), flats.data(), heap.data(),
             back.data(), dims.data());

    tt::flow_routing_d8_carve(source.data(), direction.data(),
                              filled_dem.data(), dist.data(), flats.data(),
                              dims.data());

    tt::flow_routing_targets(target.data(), source.data(), direction.data(),
                             dims.data());

    tt::flow_accumulation_edgelist(accum2.data(), source.data(), target.data(),
                                   fraction.data(), NULL, dims[0] * dims[1],
                                   dims.data());
  }

  void gradient8() {
    ProfileFunction(prof);

    tt::gradient8(gradient.data(), dem.data(), cellsize, 0, dims.data());
  }

  void gradient8_mp() {
    ProfileFunction(prof);

    tt::gradient8(gradient_mp.data(), dem.data(), cellsize, 1, dims.data());
  }

  void streamnetwork(float threshold) {
    // Process the flow direction and accumulation data to create a
    // stream network
    ptrdiff_t node_index = 1;

    for (ptrdiff_t e = (dims[0] * dims[1] - 1); e >= 0; e--) {
      ptrdiff_t u = source[e];
      ptrdiff_t v = target[e];

      if (accum2[u] >= threshold) {
        // Pixel u is in the stream network

        // Map pixel u to the node index
        stream_grid[u] = node_index++;

        if (v >= 0) {
          // u is not a sink/outlet

          // Compute distance between source and target pixels.
          // If the difference between u is 1 or dims[0], then they
          // are 4 neighbors, otherwise they are 8 neighbors.
          float w = llabs(u - v) == 1 || llabs(u - v) == dims[0]
                        ? 1.0f * cellsize
                        : SQRT2f * cellsize;

          // Add (u,v) to the edge list
          stream_edges.push_back(std::tuple(stream_grid[u], stream_grid[v], w));
        }
      }
    }

    // Fill the source and target vector by reversing stream_edges
    stream_source.reserve(stream_edges.size());
    stream_target.reserve(stream_edges.size());
    for (auto e = stream_edges.rbegin(); e != stream_edges.rend(); e++) {
      stream_source.push_back(std::get<0>(*e));
      stream_target.push_back(std::get<1>(*e));
      stream_weight.push_back(std::get<2>(*e));
    }

    stream_node_count = node_index - 1;
  }

  void runtests(bool hybrid) {
    if (hybrid) {
      route_flow_hybrid();
    } else {
      route_flow();
    }

    test_fillsinks_ge(dem.data(), filled_dem.data(), dims.data());
    test_fillsinks_filled(filled_dem.data(), dims.data());

    test_identifyflats_flats(flats.data(), filled_dem.data(), dims.data());
    test_identifyflats_sills(flats.data(), filled_dem.data(), dims.data());
    test_identifyflats_presills(flats.data(), filled_dem.data(), dims.data());

    test_gwdt_costs(costs.data(), flats.data(), dims.data());
    test_gwdt_conncomps(conncomps.data(), flats.data(), dims.data());

    // route_flow and route_flow_hybrid do not compute gradients
    gradient8();
    gradient8_mp();

    test_gradient8(gradient.data(), dem.data(), cellsize, dims.data());
    test_gradient8_mp(gradient.data(), gradient_mp.data(), dims.data());

    test_routeflowd8_direction(direction.data(), filled_dem.data(), dist.data(),
                               flats.data(), dims.data());
    test_routeflowd8_tsort(marks.data(), source.data(), direction.data(),
                           dims.data());

    test_flow_routing_targets(target.data(), source.data(), direction.data(),
                              dims.data());

    // route_flow and route_flow_hybrid only run the edgelist variant
    // of flow_accumulation, so we must run it explicitly.
    tt::flow_accumulation(accum.data(), source.data(), direction.data(), NULL,
                          dims.data());

    test_flow_accumulation_max(accum.data(), dims.data());

    test_flow_accumulation_multimethod(accum.data(), accum2.data(),
                                       dims.data());

    // Generate stream network
    streamnetwork(dims[0] * dims[1] / 20.0f);

    std::vector<float> integrand(stream_node_count, 1.0f);
    std::vector<float> integral(stream_node_count, 0.0f);
    std::vector<float> distance(dims[0] * dims[1], 0.0f);
    tt::streamquad_trapz_f32(integral.data(), integrand.data(),
                             stream_source.data(), stream_target.data(),
                             stream_weight.data(), stream_source.size());
    test_stream_distance(integral.data(), stream_grid.data(), distance.data(),
                         source.data(), target.data(), cellsize, dims.data());
  }
};

int main(int argc, char *argv[]) {
  // Default size of the DEM
  ptrdiff_t dims[2] = {100, 200};

  if (argc == 3) {
    // If two integers are passed on the command line, use those for
    // the size.

    // This will throw an exception if the inputs cannot be
    // parsed as integers.
    dims[0] = std::stoll(argv[1]);
    dims[1] = std::stoll(argv[2]);
  }

  for (uint32_t test = 0; test < 100; test++) {
    FlowRoutingData frd(dims, 10.0, test);

    if (test & 1) {
      // Run hybrid fillsinks
      frd.runtests(true);
    } else {
      // Run sequential fillsinks
      frd.runtests(false);
    }
  }
  prof.report();
}
