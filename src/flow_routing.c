#define TOPOTOOLBOX_BUILD

#include <assert.h>
#include <math.h>
#include <stddef.h>
#include <stdint.h>

#include "topotoolbox.h"

// Compute the steepest descent flow direction from the DEM with flow
// routed over flat regions by the auxiliary topography in dist.
uint8_t compute_flowdirection(ptrdiff_t i, ptrdiff_t j, float *dem, float *dist,
                              int32_t *flats, ptrdiff_t dims[2]) {
  int32_t is_flat = flats[j * dims[0] + i] & 1;
  float z = is_flat > 0 ? dist[j * dims[0] + i] : dem[j * dims[0] + i];
  uint8_t direction = 0;
  float max_gradient = 0.0;

  float chamfer[8] = {1.0, sqrtf(2.0), 1.0, sqrtf(2.0),
                      1.0, sqrtf(2.0), 1.0, sqrtf(2.0)};

  ptrdiff_t i_offset[8] = {0, 1, 1, 1, 0, -1, -1, -1};
  ptrdiff_t j_offset[8] = {1, 1, 0, -1, -1, -1, 0, 1};

  for (int32_t neighbor = 0; neighbor < 8; neighbor++) {
    uint8_t query_direction = 1 << neighbor;
    // step_downstream(idxs, query_direction, idx, dims);
    ptrdiff_t neighbor_i = i + i_offset[neighbor];
    ptrdiff_t neighbor_j = j + j_offset[neighbor];

    if (neighbor_i < 0 || neighbor_i >= dims[0] || neighbor_j < 0 ||
        neighbor_j >= dims[1]) {
      continue;
    }

    float neighbor_z = dem[neighbor_j * dims[0] + neighbor_i];
    if (neighbor_z > dem[j * dims[0] + i]) {
      // This will skip any surrounding higher neighbors of flats,
      // including sills that have higher elevations
      continue;
    }

    if (is_flat > 0) {
      neighbor_z = dist[neighbor_j * dims[0] + neighbor_i];
    }

    float g = (z - neighbor_z) / chamfer[neighbor];
    if (g > max_gradient) {
      max_gradient = g;
      direction = query_direction;
    }
  }
  return direction;
}

uint8_t compute_flowdirection_TT2(ptrdiff_t i, ptrdiff_t j, float *dem,
                                  float *dist, int32_t *flats,
                                  ptrdiff_t dims[2]) {
  uint8_t direction = 0;

  float z = dem[j * dims[0] + i];
  float max_gradient = 0.0;
  float Ggwdt = dist[j * dims[0] + i];
  float chamfer[8] = {1.0, sqrtf(2.0), 1.0, sqrtf(2.0),
                      1.0, sqrtf(2.0), 1.0, sqrtf(2.0)};

  int32_t is_flat = flats[j * dims[0] + i] & 1;

  ptrdiff_t i_offset[8] = {0, 1, 1, 1, 0, -1, -1, -1};
  ptrdiff_t j_offset[8] = {1, 1, 0, -1, -1, -1, 0, 1};

  for (int32_t neighbor = 0; neighbor < 8; neighbor++) {
    uint8_t query_direction = 1 << neighbor;

    ptrdiff_t neighbor_i = i + i_offset[neighbor];
    ptrdiff_t neighbor_j = j + j_offset[neighbor];

    if (neighbor_i < 0 || neighbor_i >= dims[0] || neighbor_j < 0 ||
        neighbor_j >= dims[1]) {
      continue;
    }
    float neighbor_z = dem[neighbor_j * dims[0] + neighbor_i];

    if (neighbor_z > z) {
      continue;
    }

    if (is_flat) {
      if (flats[neighbor_j * dims[0] + neighbor_i] & 2) {
        // neighbor is a sill
        direction = query_direction;
        break;
      } else {
        float ggwdt = dist[neighbor_j * dims[0] + neighbor_i];
        if (ggwdt < Ggwdt) {
          direction = query_direction;
          Ggwdt = ggwdt;
        }
      }
    } else {
      // Pixel is not flat
      float g = (z - neighbor_z) / chamfer[neighbor];
      if (g > max_gradient) {
        max_gradient = g;
        direction = query_direction;
      }
    }
  }
  return direction;
}

TOPOTOOLBOX_API
void flow_routing_d8_carve(ptrdiff_t *source, uint8_t *direction, float *dem,
                           float *dist, int32_t *flats, ptrdiff_t dims[2]) {
  // source[i] is the linear index of the i-th source in the
  // topologically sorted graph.
  //
  // direction[idx] is the bitfield-encoded flow direction for the
  // pixel at idx
  //
  // 1<<5 1<<6  1<<7
  // 1<<4    0  1<<0
  // 1<<3 1<<2  1<<1
  //
  // Strategy: we need to conduct a depth-first search of the graph
  // The topological sort is a reversed postordering: when we finish
  // visiting a node, we prepend it to the source array

  // Use source[--next] = idx to push the index onto the source array.
  ptrdiff_t next = dims[0] * dims[1];
  // Use source[stack_top++] = idx to push the stack
  // Use idx = source[--stack_top] top pop the stack
  ptrdiff_t stack_top = 0;

  // Initialize the directions with 0xff, which is our sentinel value
  // for unvisited nodes
  for (ptrdiff_t j = 0; j < dims[1]; j++) {
    for (ptrdiff_t i = 0; i < dims[0]; i++) {
      direction[j * dims[0] + i] = 0xff;
    }
  }

  ptrdiff_t i_offset[8] = {0, 1, 1, 1, 0, -1, -1, -1};
  ptrdiff_t j_offset[8] = {1, 1, 0, -1, -1, -1, 0, 1};

  for (ptrdiff_t j = 0; j < dims[1]; j++) {
    for (ptrdiff_t i = 0; i < dims[0]; i++) {
      if (direction[j * dims[0] + i] == 0xff) {
        // Pixel is unvisited, push it onto the stack and do a
        // depth-first search.
        // The stack should always be empty by the time we get here
        assert(stack_top == 0);
        source[stack_top++] = j * dims[0] + i;

        while (stack_top > 0) {
          // Pop a node from the stack
          ptrdiff_t node = source[--stack_top];

          if (direction[node] == 0xff) {
            // node has not been discovered yet

            // Compute its Cartesian indices
            ptrdiff_t node_i = node % dims[0];
            ptrdiff_t node_j = node / dims[0];

            uint8_t flowdir =
                compute_flowdirection(node_i, node_j, dem, dist, flats, dims);

            if (flowdir == 0) {
              // This node is a sink/outlet
              direction[node] = flowdir;
              // Prepend the node to the edge list.
              assert(next > 0);
              assert(next > stack_top);
              source[--next] = node;
              continue;
            }
            // flowdir is not an index, but 1<<index
            // Compute the index
            uint8_t v = flowdir;
            uint8_t r = 0;
            while (v >>= 1) {
              r++;
            }
            ptrdiff_t neighbor_i = node_i + i_offset[r];
            ptrdiff_t neighbor_j = node_j + j_offset[r];

            if (neighbor_i < 0 || neighbor_i >= dims[0] || neighbor_j < 0 ||
                neighbor_j >= dims[1]) {
              continue;
            }

            if (direction[neighbor_j * dims[0] + neighbor_i] == 0xff) {
              // The neighbor has not been visited yet
              // Push the node back on the stack, then push the neighbor
              assert(stack_top < dims[0] * dims[1]);
              source[stack_top++] = node;

              assert(stack_top < dims[0] * dims[1]);
              source[stack_top++] = neighbor_j * dims[0] + neighbor_i;
            } else {
              // The downstream neighbor has been visited, visit the node
              direction[node] = flowdir;
              // Prepend the node to the edge list.
              assert(next > 0);
              assert(next > stack_top);
              source[--next] = node;
            }
          } else {
            // This is a cycle
            assert(0);
          }
        }
      }
    }
  }
}

TOPOTOOLBOX_API
void flow_routing_targets(ptrdiff_t *target, ptrdiff_t *source,
                          uint8_t *direction, ptrdiff_t dims[2]) {
  ptrdiff_t offsets[8] = {dims[0],  dims[0] + 1,  1,  -dims[0] + 1,
                          -dims[0], -dims[0] - 1, -1, dims[0] - 1};

  for (ptrdiff_t j = 0; j < dims[1]; j++) {
    for (ptrdiff_t i = 0; i < dims[0]; i++) {
      ptrdiff_t node = source[j * dims[0] + i];

      uint8_t flowdir = direction[node];

      uint8_t v = flowdir;
      uint8_t r = 0;
      while (v >>= 1) {
        r++;
      }

      target[j * dims[0] + i] = node + offsets[r];
    }
  }
}
