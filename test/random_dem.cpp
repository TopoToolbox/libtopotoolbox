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

// Include topotoolbox in its own namespace to help prevent
// naming conflicts in the global scope.
namespace tt {
extern "C" {
#include "topotoolbox.h"
}
}  // namespace tt

extern "C" {
#include "flow.h"
#include "grid.h"
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
int32_t test_fillsinks_ge(Grid original, Grid filled) {
  Grid sub = GridSubtract(filled, original);
  assert(GridAllNonnegative(sub));
  GridFree(&sub);
  return 0;
}

/*
  No pixel in the filled DEM should be completely surrounded by pixels higher
  than it.
 */
int32_t test_fillsinks_filled(Grid filled) {
  Grid sinks = GridLabelSinks(filled);
  assert(GridAny(sinks) == 0);
  GridFree(&sinks);
  return 0;
}

/*
  Every pixel not on the boundary with no lower neighbors and fewer than
  8 higher neighbors should be labeled a flat. Likewise, every flat
  should have no lower neighbors and fewer than 8 higher neighbors.
 */
int32_t test_identifyflats_flats(Grid flats, Grid dem) {
  ptrdiff_t col_offset[8] = {-1, -1, -1, 0, 1, 1, 1, 0};
  ptrdiff_t row_offset[8] = {1, 0, -1, -1, -1, 0, 1, 1};

  int32_t *flat_data = (int32_t *)flats.data;
  float *dem_data = (float *)dem.data;

  ptrdiff_t dims[2] = {dem.dims[0], dem.dims[1]};

  for (ptrdiff_t j = 0; j < dims[1]; j++) {
    for (ptrdiff_t i = 0; i < dims[0]; i++) {
      float z = dem_data[i + dims[0] * j];
      int32_t flat = flat_data[i + dims[0] * j];

      int32_t up_neighbor_count = 0;
      int32_t down_neighbor_count = 0;
      for (int32_t neighbor = 0; neighbor < 8; neighbor++) {
        ptrdiff_t neighbor_i = i + row_offset[neighbor];
        ptrdiff_t neighbor_j = j + col_offset[neighbor];

        if (neighbor_i < 0 || neighbor_i >= dims[0] || neighbor_j < 0 ||
            neighbor_j >= dims[1]) {
          continue;
        }

        float neighbor_height = dem_data[neighbor_i + dims[0] * neighbor_j];

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
int32_t test_identifyflats_sills(Grid flats, Grid dem) {
  ptrdiff_t j_offset[8] = {-1, -1, -1, 0, 1, 1, 1, 0};
  ptrdiff_t i_offset[8] = {1, 0, -1, -1, -1, 0, 1, 1};

  int32_t *flat_data = (int32_t *)flats.data;
  float *dem_data = (float *)dem.data;

  ptrdiff_t dims[2] = {dem.dims[0], dem.dims[1]};

  for (ptrdiff_t j = 0; j < dims[1]; j++) {
    for (ptrdiff_t i = 0; i < dims[0]; i++) {
      float z = dem_data[i + j * dims[0]];
      int32_t flat = flat_data[i + j * dims[0]];

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

        float neighbor_height = dem_data[neighbor_i + dims[0] * neighbor_j];
        int32_t neighbor_flat = flat_data[neighbor_i + dims[0] * neighbor_j];

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
int32_t test_identifyflats_presills(Grid flats, Grid dem) {
  ptrdiff_t j_offset[8] = {-1, -1, -1, 0, 1, 1, 1, 0};
  ptrdiff_t i_offset[8] = {1, 0, -1, -1, -1, 0, 1, 1};

  int32_t *flat_data = (int32_t *)flats.data;
  float *dem_data = (float *)dem.data;

  ptrdiff_t dims[2] = {dem.dims[0], dem.dims[1]};

  for (ptrdiff_t j = 0; j < dims[1]; j++) {
    for (ptrdiff_t i = 0; i < dims[0]; i++) {
      float z = dem_data[i + j * dims[0]];
      int32_t flat = flat_data[i + j * dims[0]];

      int32_t equal_neighbor_sills = 0;

      for (int32_t neighbor = 0; neighbor < 8; neighbor++) {
        ptrdiff_t neighbor_i = i + i_offset[neighbor];
        ptrdiff_t neighbor_j = j + j_offset[neighbor];

        if (neighbor_i < 0 || neighbor_i >= dims[0] || neighbor_j < 0 ||
            neighbor_j >= dims[1]) {
          continue;
        }

        float neighbor_height = dem_data[neighbor_i + dims[0] * neighbor_j];
        int32_t neighbor_flat = flat_data[neighbor_i + dims[0] * neighbor_j];

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
int32_t test_gwdt_costs(Grid costs, Grid flats) {
  for (ptrdiff_t j = 0; j < flats.dims[1]; j++) {
    for (ptrdiff_t i = 0; i < flats.dims[0]; i++) {
      float cost = ((float *)costs.data)[i + costs.dims[0] * j];
      int32_t flat = ((int32_t *)flats.data)[i + flats.dims[0] * j];

      assert((cost >= 0) && (((flat & 1) > 0) == (cost > 0)));
    }
  }

  return 0;
}

/*
  Flats neighboring other flats should all have the same nonzero
  connected component label
 */
int32_t test_gwdt_conncomps(Grid conncomps, Grid flats) {
  ptrdiff_t j_offset[8] = {-1, -1, -1, 0, 1, 1, 1, 0};
  ptrdiff_t i_offset[8] = {1, 0, -1, -1, -1, 0, 1, 1};

  for (ptrdiff_t j = 0; j < flats.dims[1]; j++) {
    for (ptrdiff_t i = 0; i < flats.dims[0]; i++) {
      ptrdiff_t label =
          ((ptrdiff_t *)conncomps.data)[i + conncomps.dims[0] * j];
      int32_t flat = ((int32_t *)flats.data)[i + flats.dims[0] * j];

      if (!(flat & 1)) {
        continue;
      }

      for (int32_t neighbor = 0; neighbor < 8; neighbor++) {
        ptrdiff_t neighbor_i = i + i_offset[neighbor];
        ptrdiff_t neighbor_j = j + j_offset[neighbor];

        if (neighbor_i < 0 || neighbor_i >= flats.dims[0] || neighbor_j < 0 ||
            neighbor_j >= flats.dims[1]) {
          continue;
        }

        int32_t neighbor_flat =
            ((int32_t *)flats.data)[neighbor_i + flats.dims[0] * neighbor_j];
        ptrdiff_t neighbor_label =
            ((ptrdiff_t *)
                 conncomps.data)[neighbor_i + conncomps.dims[0] * neighbor_j];
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
int32_t test_gwdt(Grid dist, Grid prev, Grid costs, Grid flats) {
  ptrdiff_t dims[2] = {dist.dims[0], dist.dims[1]};

  for (ptrdiff_t j = 0; j < dims[1]; j++) {
    for (ptrdiff_t i = 0; i < dims[0]; i++) {
      float d = ((float *)dist.data)[j * dims[0] + i];
      int32_t flat = ((int32_t *)flats.data)[j * dims[0] + i];

      assert(d >= 0);

      if (flat & 4) {
        assert(d == 1.0);
      } else if (flat & 1) {
        assert(d > 1.0);

        ptrdiff_t parent = ((ptrdiff_t *)prev.data)[j * dims[0] + i];
        ptrdiff_t parent_j = parent / dims[0];
        ptrdiff_t parent_i = parent % dims[0];
        float chamfer = (parent_j != j) && (parent_i != i) ? SQRT2f : 1.0f;
        float proposed_dist = ((float *)dist.data)[parent] +
                              chamfer *
                                  (((float *)costs.data)[j * dims[0] + i] +
                                   ((float *)costs.data)[parent]) /
                                  2;
        assert(proposed_dist == d);
      }
    }
  }

  return 0;
}
/*
  For each cell, the saved gradient should be greater than or equal to
  the signed gradient to all 8 neighboring cells.
*/
int32_t test_gradient8(Grid gradient, Grid dem) {
  ptrdiff_t i_offset[8] = {0, 1, 1, 1, 0, -1, -1, -1};
  ptrdiff_t j_offset[8] = {1, 1, 0, -1, -1, -1, 0, 1};

  ptrdiff_t dims[2] = {dem.dims[0], dem.dims[1]};

  float *gradient_data = (float *)gradient.data;
  float *dem_data = (float *)dem.data;

  for (ptrdiff_t j = 0; j < dims[1]; j++) {
    for (ptrdiff_t i = 0; i < dims[0]; i++) {
      float max_gradient = gradient_data[j * dims[0] + i];

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

          horizontal_dist = dem.cellsize;
          horizontal_dist *= neighbor_i != i && neighbor_j != j ? SQRT2 : 1.0f;

          vertical_dist = dem_data[j * dims[0] + i] -
                          dem_data[neighbor_j * dims[0] + neighbor_i];

          local_gradient = vertical_dist / horizontal_dist;
          assert(max_gradient >= local_gradient);
        }
      }
    }
  }
  return 0;
}
/*
  Computing the gradient8 using OpenMP should yield the same result as without.
*/
int32_t test_gradient8_mp(Grid gradient, Grid gradient_mp) {
  assert(GridEq(gradient, gradient_mp));
  return 0;
}
/*
  Flow direction should point downstream or across flats
 */
int32_t test_routeflowd8_direction(Grid direction, Grid filled_dem) {
  // 1<<5 1<<6  1<<7
  // 1<<4    0  1<<0
  // 1<<3 1<<2  1<<1
  ptrdiff_t i_offset[8] = {0, 1, 1, 1, 0, -1, -1, -1};
  ptrdiff_t j_offset[8] = {1, 1, 0, -1, -1, -1, 0, 1};

  ptrdiff_t dims[2] = {direction.dims[0], direction.dims[1]};

  uint8_t *direction_data = (uint8_t *)direction.data;
  float *filled_data = (float *)filled_dem.data;

  for (ptrdiff_t j = 0; j < dims[1]; j++) {
    for (ptrdiff_t i = 0; i < dims[0]; i++) {
      uint8_t flowdir = direction_data[j * dims[0] + i];

      if (flowdir == 0) {
        // Pixel is a sink
        continue;
      }

      float z = filled_data[j * dims[0] + i];
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
      assert(filled_data[neighbor_j * dims[0] + neighbor_i] <= z);
    }
  }
  return 0;
}

int32_t test_tsort(Flow fd) {
  // source and target arrays should be topologically sorted. All
  // incoming edges to a node must occur before any outgoing edges.

  int32_t flag = 0;

  Grid outdegree = GridCreate(GridU8, NULL, fd.cellsize, fd.dims);
  for (ptrdiff_t e = 0; e < fd.count; e++) {
    ptrdiff_t u = fd.source[e];
    ptrdiff_t v = fd.target[e];

    ((uint8_t *)outdegree.data)[u] += 1;

    if (((uint8_t *)outdegree.data)[v] > 0) {
      flag = 1;
    }
  }
  assert(flag == 0);
  GridFree(&outdegree);
  return flag;
}

int32_t test_flow_accumulation_max(Grid acc) {
  Grid max = GridFill(acc.dims[0] * acc.dims[1], acc.cellsize, acc.dims);
  Grid sub = GridSubtract(max, acc);
  assert(GridAllNonnegative(sub));

  GridFree(&sub);
  GridFree(&max);
  return 0;
}

/*
  The target corresponding to source should be equal to the source
  plus the correct offset in the flow direction.
 */
int32_t test_flow_routing_targets(ptrdiff_t *target, ptrdiff_t *source,
                                  uint8_t *direction, ptrdiff_t edge_count,
                                  ptrdiff_t dims[2]) {
  ptrdiff_t i_offset[8] = {0, 1, 1, 1, 0, -1, -1, -1};
  ptrdiff_t j_offset[8] = {1, 1, 0, -1, -1, -1, 0, 1};

  for (ptrdiff_t e = 0; e < edge_count; e++) {
    ptrdiff_t u = source[e];
    ptrdiff_t u_i = u % dims[0];
    ptrdiff_t u_j = u / dims[0];

    assert(u_i >= 0 && u_i < dims[0]);
    assert(u_j >= 0 && u_j < dims[1]);

    uint8_t flowdir = direction[u];

    if (flowdir != 0) {
      uint8_t r = 0;
      while (flowdir >>= 1) {
        r++;
      }

      ptrdiff_t neighbor_i = u_i + i_offset[r];
      ptrdiff_t neighbor_j = u_j + j_offset[r];

      assert(neighbor_i >= 0 && neighbor_i < dims[0]);
      assert(neighbor_j >= 0 && neighbor_j < dims[1]);

      ptrdiff_t v = neighbor_j * dims[0] + neighbor_i;

      assert(v == target[e]);
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
                             ptrdiff_t edge_count, ptrdiff_t dims[2]) {
  // Compute the upstream distance manually for all edges
  for (ptrdiff_t edge = edge_count - 1; edge >= 0; edge--) {
    ptrdiff_t u = source[edge];
    ptrdiff_t v = target[edge];

    float w = llabs(u - v) == 1 || llabs(u - v) == dims[0] ? 1.0f * cellsize
                                                           : SQRT2f * cellsize;
    distance[u] = distance[v] + w;
  }

  for (ptrdiff_t j = 0; j < dims[1]; j++) {
    for (ptrdiff_t i = 0; i < dims[0]; i++) {
      ptrdiff_t node = stream_grid[j * dims[0] + i];
      if (node > 0) {
        assert(node_distance[node] == distance[j * dims[0] + i]);
      }
    }
  }
  return 0;
}

/*
  A source pixel and its downstream neighbor should be part of the
  same basin, which should have an index greater than 0.
 */
int32_t test_drainagebasins(uint32_t *basins, ptrdiff_t *source,
                            ptrdiff_t *target, ptrdiff_t edge_count,
                            ptrdiff_t dims[2]) {
  for (ptrdiff_t e = 0; e < edge_count; e++) {
    ptrdiff_t src = source[e];
    ptrdiff_t tgt = target[e];

    assert(src < dims[0] * dims[1]);
    assert(src >= 0);
    assert(tgt < dims[0] * dims[1]);
    assert(tgt >= 0);
    assert(basins[src] == basins[tgt]);
    assert(basins[src] > 0);
  }
  return 0;
}

int32_t test_propagatevalues(float *pvF32, double *pvF64, uint8_t *pvU8,
                             uint32_t *pvU32, uint64_t *pvU64, int8_t *pvI8,
                             int32_t *pvI32, int64_t *pvI64,
                             ptrdiff_t node_count) {
  for (ptrdiff_t v = 0; v < node_count; v++) {
    assert(pvF32[v] == pvF64[v] && pvF64[v] == pvU8[v] && pvU32[v] == pvU8[v] &&
           pvU64[v] == pvU32[v] && (uint64_t)pvI8[v] == pvU64[v] &&
           pvI32[v] == pvI8[v] && pvI64[v] == pvI32[v]);
  }
  return 0;
}

int32_t test_lowerenv_convex(float *g, float *z, float *d, uint8_t *knickpoints,
                             ptrdiff_t *source, ptrdiff_t *target,
                             ptrdiff_t edge_count, ptrdiff_t node_count) {
  for (ptrdiff_t e = 0; e < edge_count; e++) {
    ptrdiff_t u = source[e];
    ptrdiff_t v = target[e];

    g[u] = (z[u] - z[v]) / (d[u] - d[v]);
  }

  for (ptrdiff_t e = 0; e < edge_count; e++) {
    ptrdiff_t u = source[e];
    ptrdiff_t v = target[e];

    assert(((g[u] - g[v]) >= -1e-6) || knickpoints[v]);
  }
  return 0;
}

// drainagebasins should return the same drainage basins computed
// using an appropriately initialized traverse_up_u32_or_and.
int32_t test_db_traverse(uint32_t *basins, ptrdiff_t *source, ptrdiff_t *target,
                         ptrdiff_t edge_count, ptrdiff_t dims[2]) {
  std::vector<uint32_t> ones(edge_count, 0xffffffff);           // Edge weights
  std::vector<uint32_t> traverse_basins(dims[0] * dims[1], 0);  // Basin labels
  std::vector<uint8_t> indegree(dims[0] * dims[1], 0);
  std::vector<uint8_t> outdegree(dims[0] * dims[1], 0);

  // Identify outlets
  tt::edgelist_degree(indegree.data(), outdegree.data(), source, target,
                      dims[0] * dims[1], edge_count);

  for (ptrdiff_t j = 0; j < dims[1]; j++) {
    for (ptrdiff_t i = 0; i < dims[0]; i++) {
      ptrdiff_t p = j * dims[0] + i;
      if (indegree[p] > 0 && outdegree[p] == 0) {
        traverse_basins[p] = (uint32_t)basins[p];
      }
    }
  }
  tt::traverse_up_u32_or_and(traverse_basins.data(), ones.data(), source,
                             target, edge_count);

  for (ptrdiff_t j = 0; j < dims[1]; j++) {
    for (ptrdiff_t i = 0; i < dims[0]; i++) {
      assert(basins[j * dims[0] + i] == traverse_basins[j * dims[0] + i]);
    }
  }

  return 0;
}

struct FlowRoutingData {
  std::array<ptrdiff_t, 2> dims;
  float cellsize;

  // input data
  Grid dem;
  Grid bc;
  Grid fraction;

  // fillsinks
  Grid filled_dem;

  // identifyflats
  Grid flats;

  // gwdt_compute_costs
  Grid costs;
  Grid conncomps;

  // gwdt
  Grid dist;
  Grid heap;
  Grid back;
  Grid prev;

  // gradient8
  Grid gradient;
  Grid gradient_mp;

  // flow_routing_d8_carve
  Grid direction;

  // flow_routing_targets
  Flow fd;

  // flow_accumulation
  Grid accum;

  // stream network
  std::vector<ptrdiff_t> stream_source;
  std::vector<ptrdiff_t> stream_target;
  std::vector<float> stream_weight;
  std::vector<ptrdiff_t> stream_grid;
  ptrdiff_t stream_node_count;

  Grid basins;

  FlowRoutingData(ptrdiff_t input_dims[2], float cs, uint32_t seed)
      : dims({input_dims[0], input_dims[1]}),
        cellsize(cs),
        stream_grid(dims[0] * dims[1]) {
    // Initialize DEM, boundary conditions for fillsinks and fraction
    dem = GridRandom(seed, cellsize, dims.data());
    bc = GridCreate(GridU8, NULL, cellsize, dims.data());
    fraction = GridFill(1.0f, cellsize, dims.data());

    for (ptrdiff_t col = 0; col < dims[1]; col++) {
      for (ptrdiff_t row = 0; row < dims[0]; row++) {
        if (row == 0 || row == dims[0] - 1 || col == 0 || col == dims[1] - 1) {
          ((uint8_t *)bc.data)[col * dims[0] + row] = 1;
        } else {
          ((uint8_t *)bc.data)[col * dims[0] + row] = 0;
        }
      }
    }
  }

  ~FlowRoutingData() {
    GridFree(&dem);
    GridFree(&bc);
    GridFree(&fraction);
    GridFree(&filled_dem);
    GridFree(&flats);
    GridFree(&costs);
    GridFree(&conncomps);
    GridFree(&dist);
    GridFree(&prev);
    GridFree(&back);
    GridFree(&heap);
    GridFree(&gradient);
    GridFree(&gradient_mp);
    GridFree(&direction);
    GridFree(&accum);
    GridFree(&basins);
    FlowFree(&fd);
  }

  void fillsinks_hybrid() {
    ProfileFunction(prof);
    filled_dem = GridFillsinks(dem);
  }

  void fillsinks() {
    ProfileFunction(prof);
    filled_dem = GridCreate(GridF32, NULL, dem.cellsize, dem.dims);

    tt::fillsinks((float *)filled_dem.data, (float *)dem.data,
                  (uint8_t *)bc.data, dims.data());
  }

  void route_flow(bool hybrid) {
    ProfileFunction(prof);

    if (hybrid) {
      fillsinks_hybrid();
    } else {
      fillsinks();
    }

    flats = GridIdentifyFlats(filled_dem);

    costs = GridZerosLike(dem);
    conncomps = GridCreate(GridIdx, NULL, dem.cellsize, dem.dims);

    tt::gwdt_computecosts((float *)costs.data, (ptrdiff_t *)conncomps.data,
                          (int32_t *)flats.data, (float *)dem.data,
                          (float *)filled_dem.data, dims.data());

    dist = GridZerosLike(costs);
    prev = GridCreate(GridIdx, NULL, dist.cellsize, dist.dims);
    heap = GridCreate(GridIdx, NULL, dist.cellsize, dist.dims);
    back = GridCreate(GridIdx, NULL, dist.cellsize, dist.dims);
    tt::gwdt((float *)dist.data, (ptrdiff_t *)prev.data, (float *)costs.data,
             (int32_t *)flats.data, (ptrdiff_t *)heap.data,
             (ptrdiff_t *)back.data, dims.data());

    fd = FlowInit(dem);

    direction = GridCreate(GridU8, NULL, dem.cellsize, dem.dims);
    tt::flow_routing_d8_carve((ptrdiff_t *)fd.stream, (uint8_t *)direction.data,
                              (float *)filled_dem.data, (float *)dist.data,
                              (int32_t *)flats.data, dims.data(), 0);

    fd.count = tt::flow_routing_d8_edgelist(
        (ptrdiff_t *)fd.source, (ptrdiff_t *)fd.target, (ptrdiff_t *)fd.stream,
        (uint8_t *)direction.data, dims.data(), 0);

    for (ptrdiff_t e = 0; e < fd.count; e++) {
      fd.fraction[e] = 1.0;
    }

    accum = FlowAccumulation(fd);
  }

  void gradient8() {
    ProfileFunction(prof);
    gradient = GridGradient8(dem, 0);
  }

  void gradient8_mp() {
    ProfileFunction(prof);
    gradient_mp = GridGradient8(dem, 1);
  }

  void streamnetwork(float threshold) {
    // Process the flow direction and accumulation data to create a
    // stream network
    ptrdiff_t node_index = 0;

    // Label all of the nodes in the stream network with indices into
    // the node attribute list.
    for (ptrdiff_t j = 0; j < dims[1]; j++) {
      for (ptrdiff_t i = 0; i < dims[0]; i++) {
        ptrdiff_t u = ((ptrdiff_t *)fd.stream)[j * dims[0] + i];
        if (((float *)accum.data)[u] >= threshold) {
          // Pixel u is in the stream network
          stream_grid[u] = node_index++;
        } else {
          stream_grid[u] = -1;
        }
      }
    }

    // Find edges that are part of the stream network
    for (ptrdiff_t e = 0; e < fd.count; e++) {
      ptrdiff_t u = fd.source[e];
      ptrdiff_t v = fd.target[e];

      // Compute distance between source and target pixels. If the
      // difference between u and v is 1 or dims[0], then they are
      // cardinal neighbors, otherwise they are diagonal neighbors.
      float w = llabs(u - v) == 1 || llabs(u - v) == dims[0]
                    ? 1.0f * cellsize
                    : SQRT2f * cellsize;

      if (((float *)accum.data)[u] >= threshold) {
        // Pixel u is in the stream network

        // Add this edge to the stream network, but map u and v to
        // their node attribute indices.
        stream_source.push_back(stream_grid[u]);
        stream_target.push_back(stream_grid[v]);
        stream_weight.push_back(w);
      }
    }
    stream_node_count = node_index;
  }

  void drainagebasins() {
    ProfileFunction(prof);

    basins = FlowDrainageBasins(fd);
  }

  void runtests(bool hybrid) {
    route_flow(hybrid);

    test_fillsinks_ge(dem, filled_dem);
    test_fillsinks_filled(filled_dem);

    test_identifyflats_flats(flats, filled_dem);
    test_identifyflats_sills(flats, filled_dem);
    test_identifyflats_presills(flats, filled_dem);

    test_gwdt_costs(costs, flats);
    test_gwdt_conncomps(conncomps, flats);

    test_gwdt(dist, prev, costs, flats);

    // route_flow and route_flow_hybrid do not compute gradients
    gradient8();
    gradient8_mp();

    test_gradient8(gradient, dem);
    test_gradient8_mp(gradient, gradient_mp);

    test_routeflowd8_direction(direction, filled_dem);

    test_tsort(fd);

    test_flow_routing_targets((ptrdiff_t *)fd.target, (ptrdiff_t *)fd.source,
                              (uint8_t *)direction.data, fd.count, dims.data());

    // Compute and test drainage basins
    drainagebasins();
    test_drainagebasins((uint32_t *)basins.data, (ptrdiff_t *)fd.source,
                        (ptrdiff_t *)fd.target, fd.count, dims.data());
    test_db_traverse((uint32_t *)basins.data, (ptrdiff_t *)fd.source,
                     (ptrdiff_t *)fd.target, fd.count, dims.data());

    test_flow_accumulation_max(accum);

    // Generate stream network
    streamnetwork(dims[0] * dims[1] / 20.0f);

    std::vector<float> integrand(stream_node_count, 1.0f);
    std::vector<float> integral(stream_node_count, 0.0f);
    std::vector<float> distance(dims[0] * dims[1], 0.0f);
    tt::streamquad_trapz_f32(integral.data(), integrand.data(),
                             stream_source.data(), stream_target.data(),
                             stream_weight.data(), stream_source.size());
    test_stream_distance(integral.data(), stream_grid.data(), distance.data(),
                         (ptrdiff_t *)fd.source, (ptrdiff_t *)fd.target,
                         cellsize, fd.count, dims.data());

    std::vector<float> pvF32(stream_node_count, 0.0f);
    std::vector<double> pvF64(stream_node_count, 0.0);
    std::vector<uint8_t> pvU8(stream_node_count, 0);
    std::vector<uint32_t> pvU32(stream_node_count, 0);
    std::vector<uint64_t> pvU64(stream_node_count, 0);
    std::vector<int8_t> pvI8(stream_node_count, 0);
    std::vector<int32_t> pvI32(stream_node_count, 0);
    std::vector<int64_t> pvI64(stream_node_count, 0);

    for (ptrdiff_t i = 0; i < stream_node_count; i++) {
      uint8_t v = (uint64_t)i & 0x7f;
      pvF32[i] = v;
      pvF64[i] = v;
      pvU8[i] = v;
      pvU32[i] = v;
      pvU64[i] = v;
      pvI8[i] = v;
      pvI32[i] = v;
      pvI64[i] = v;
    }

    tt::propagatevaluesupstream_f32(pvF32.data(), stream_source.data(),
                                    stream_target.data(), stream_source.size());
    tt::propagatevaluesupstream_f64(pvF64.data(), stream_source.data(),
                                    stream_target.data(), stream_source.size());
    tt::propagatevaluesupstream_u8(pvU8.data(), stream_source.data(),
                                   stream_target.data(), stream_source.size());
    tt::propagatevaluesupstream_u32(pvU32.data(), stream_source.data(),
                                    stream_target.data(), stream_source.size());
    tt::propagatevaluesupstream_u64(pvU64.data(), stream_source.data(),
                                    stream_target.data(), stream_source.size());
    tt::propagatevaluesupstream_i8(pvI8.data(), stream_source.data(),
                                   stream_target.data(), stream_source.size());
    tt::propagatevaluesupstream_i32(pvI32.data(), stream_source.data(),
                                    stream_target.data(), stream_source.size());
    tt::propagatevaluesupstream_i64(pvI64.data(), stream_source.data(),
                                    stream_target.data(), stream_source.size());
    test_propagatevalues(pvF32.data(), pvF64.data(), pvU8.data(), pvU32.data(),
                         pvU64.data(), pvI8.data(), pvI32.data(), pvI64.data(),
                         stream_node_count);

    std::vector<uint8_t> kn(stream_node_count, 0);
    std::vector<uint8_t> onenvelope(stream_node_count, 0);
    std::vector<ptrdiff_t> ix(stream_node_count, 0);
    std::vector<float> z(stream_node_count, 1.0f);
    std::vector<float> d(stream_node_count, 0.0f);
    std::vector<float> w(stream_source.size(), 0.0f);

    for (ptrdiff_t j = 0; j < dims[1]; j++) {
      for (ptrdiff_t i = 0; i < dims[0]; i++) {
        ptrdiff_t node = j * dims[0] + i;
        if (stream_grid[node] >= 0) {
          z[stream_grid[node]] = ((float *)dem.data)[node];
          d[stream_grid[node]] = distance[node];
        }
      }
    }

    for (ptrdiff_t i = 0; i < stream_node_count; i++) {
      kn[i] = (i % 10) == 0;
    }

    // imposemin
    tt::traverse_down_f32_min_add(z.data(), w.data(), stream_source.data(),
                                  stream_target.data(), stream_source.size());

    tt::lowerenv(z.data(), kn.data(), d.data(), ix.data(), onenvelope.data(),
                 stream_source.data(), stream_target.data(),
                 stream_source.size(), stream_node_count);

    std::vector<float> g(stream_node_count, 0.0f);
    test_lowerenv_convex(g.data(), z.data(), d.data(), kn.data(),
                         stream_source.data(), stream_target.data(),
                         stream_source.size(), stream_node_count);
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
