#define TOPOTOOLBOX_BUILD

#include <assert.h>
#include <math.h>
#include <stddef.h>
#include <stdint.h>

#include "helpers/priority_queue.h"
#include "topotoolbox.h"

#define SQRT2 1.41421356237309504880f

// Compute the steepest descent flow direction from the DEM with flow
// routed over flat regions by the auxiliary topography in dist.
uint8_t compute_flowdirection(ptrdiff_t i, ptrdiff_t j, float *dem, float *dist,
                              int32_t *flats, ptrdiff_t dims[2],
                              unsigned int order) {
  int32_t is_flat = flats[j * dims[0] + i] & 1;
  float z = is_flat > 0 ? dist[j * dims[0] + i] : dem[j * dims[0] + i];
  uint8_t direction = 0;
  float max_gradient = 0.0;

  /*
    We need the index offsets so we can easily check the bounds of the
    array. To index in the correct order, we need to reverse the
    direction of indexing when working with row-major arrays. This
    amounts to swapping the i_offsets and j_offsets when the input
    arrays are row-major (order == 1). ijoffsets[order] is the correct
    index offset in the first dimension and ijoffsets[order ^ 1] is
    the corect index offset in the second dimension.

    for (int neighbor = 0; neighbor < 8; neighbor++) {
      ptrdiff_t i_offset = ijoffsets[order & 1][neighbor];
      ptrdiff_t j_offset = ijoffsets[(order ^ 1) & 1][neighbor];
      // ...
    }
   */

  ptrdiff_t ij_offsets[2][8] = {{0, 1, 1, 1, 0, -1, -1, -1},
                                {1, 1, 0, -1, -1, -1, 0, 1}};

  float chamfer[8] = {1.0, sqrtf(2.0), 1.0, sqrtf(2.0),
                      1.0, sqrtf(2.0), 1.0, sqrtf(2.0)};

  for (int32_t neighbor = 0; neighbor < 8; neighbor++) {
    uint8_t query_direction = 1 << neighbor;
    // step_downstream(idxs, query_direction, idx, dims);
    ptrdiff_t neighbor_i = i + ij_offsets[order & 1][neighbor];
    ptrdiff_t neighbor_j = j + ij_offsets[(order ^ 1) & 1][neighbor];

    if (neighbor_i < 0 || neighbor_i >= dims[0] || neighbor_j < 0 ||
        neighbor_j >= dims[1]) {
      // This check is important, because we are not yet sure whether
      // it is in bounds or not.
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
void flow_routing_d8_carve(ptrdiff_t *node, uint8_t *direction, float *dem,
                           float *dist, int32_t *flats, ptrdiff_t dims[2],
                           unsigned int order) {
  // node contains an array of dims[0] * dims[1] linear pixel indices
  // into dem. These indices are sorted topologically, so that if
  // there is an edge from node u to node v, u comes before v in the
  // array.
  //
  // direction[i] is the bitfield-encoded flow direction for the
  // pixel at i
  //
  // 1<<5 1<<6  1<<7
  // 1<<4    0  1<<0
  // 1<<3 1<<2  1<<1
  //
  // To construct a topological ordering of the flow graph, we conduct
  // a depth first traversal of the flow graph. The topological order
  // is given by a reversed postorder of the nodes encountered during
  // the traversal. As a result, the node array fills up from the
  // bottom.
  //
  // Use node[--next] = u to append vertex u onto the node list.
  ptrdiff_t next = dims[0] * dims[1];

  // To implement the depth first traversal, we need a stack. The
  // stack can only hold up to as many vertices as have not yet been
  // assigned to the node list.  We can therefore using the top of the
  // node array to hold our stack.
  //
  // Use source[stack_top++] = u to push vertex u onto the stack
  // Use idx = source[--stack_top] top pop the stack
  ptrdiff_t stack_top = 0;

  // Initialize the directions with 0xff, which is our sentinel value
  // for unvisited nodes
  for (ptrdiff_t j = 0; j < dims[1]; j++) {
    for (ptrdiff_t i = 0; i < dims[0]; i++) {
      direction[j * dims[0] + i] = 0xff;
    }
  }

  ptrdiff_t strides[2] = {0};
  if (order & 1) {
    // row-major
    strides[0] = dims[0];
    strides[1] = 1;
  } else {
    strides[0] = 1;
    strides[1] = dims[0];
  }

  ptrdiff_t offsets[8] = {strides[1],  strides[0] + strides[1],
                          strides[0],  strides[0] - strides[1],
                          -strides[1], -strides[0] - strides[1],
                          -strides[0], -strides[0] + strides[1]};

  for (ptrdiff_t j = 0; j < dims[1]; j++) {
    for (ptrdiff_t i = 0; i < dims[0]; i++) {
      if (direction[j * dims[0] + i] == 0xff) {
        // Pixel is unvisited, push it onto the stack and do a
        // depth-first search.
        // The stack should always be empty by the time we get here
        assert(stack_top == 0);
        node[stack_top++] = j * dims[0] + i;

        while (stack_top > 0) {
          // Pop a node from the stack
          ptrdiff_t u = node[--stack_top];

          if (direction[u] == 0xff) {
            // node has not been discovered yet

            // Compute its Cartesian indices
            ptrdiff_t u_i = u % dims[0];
            ptrdiff_t u_j = u / dims[0];

            uint8_t flowdir =
                compute_flowdirection(u_i, u_j, dem, dist, flats, dims, order);

            if (flowdir == 0) {
              // This node is a sink/outlet
              direction[u] = flowdir;
              // Prepend the node to the edge list.
              assert(next > 0);
              assert(next > stack_top);
              node[--next] = u;
              continue;
            }
            // flowdir is not an index, but 1<<index
            // Compute the index
            uint8_t d = flowdir;
            uint8_t r = 0;
            while (d >>= 1) {
              r++;
            }
            ptrdiff_t v = u + offsets[r];

            if (direction[v] == 0xff) {
              // The neighbor has not been visited yet
              // Push the node back on the stack, then push the neighbor
              assert(stack_top < next);
              node[stack_top++] = u;

              assert(stack_top < next);
              node[stack_top++] = v;
            } else {
              // The downstream neighbor has been visited, visit the node
              direction[u] = flowdir;
              // Prepend the node to the edge list.
              assert(next > 0);
              assert(next > stack_top);
              node[--next] = u;
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
ptrdiff_t flow_routing_d8_edgelist(ptrdiff_t *source, ptrdiff_t *target,
                                   ptrdiff_t *node, uint8_t *direction,
                                   ptrdiff_t dims[2], unsigned int order) {
  // For an array of size {m,n}, strides is {n,1} if row-major, {1, m} if
  // column-major Note dims = {m,n} for column-major, {n,m} for row-major
  ptrdiff_t strides[2] = {0};
  if (order & 1) {
    // row-major
    strides[0] = dims[0];
    strides[1] = 1;
  } else {
    strides[0] = 1;
    strides[1] = dims[0];
  }

  /*
    These are the offsets of the linear addresses for each neighbor
    of the central cell, regardless of the memory order,

    ---------------------------------------
    | -s[0] - s[1] | -s[0] | -s[0] + s[1] |
    |       - s[1] |   0   |         s[1] |
    |  s[0] - s[1] |  s[0] |  s[0] + s[1] |
    ---------------------------------------

    assuming that s[0] is the stride between rows and s[1] is the
    stride between columns.
   */
  ptrdiff_t offsets[8] = {strides[1],  strides[0] + strides[1],
                          strides[0],  strides[0] - strides[1],
                          -strides[1], -strides[0] - strides[1],
                          -strides[0], -strides[0] + strides[1]};

  ptrdiff_t edge_count = 0;
  for (ptrdiff_t j = 0; j < dims[1]; j++) {
    for (ptrdiff_t i = 0; i < dims[0]; i++) {
      ptrdiff_t u = node[j * dims[0] + i];

      uint8_t flowdir = direction[u];

      if (flowdir != 0) {
        uint8_t v = flowdir;
        uint8_t r = 0;
        while (v >>= 1) {
          r++;
        }

        source[edge_count] = u;
        target[edge_count++] = u + offsets[r];
      }
    }
  }
  return edge_count;
}

///////////////////////////////////
// Topological sorting of edge sets

static uint8_t next_neighbor(uint8_t d) {
  uint8_t n = 0;
  while (d >>= 1) {
    n++;
  }
  return n;
}

TOPOTOOLBOX_API
void flow_routing_tsort(ptrdiff_t *stream, ptrdiff_t *source, ptrdiff_t *target,
                        float *sorted_weight, ptrdiff_t *stack,
                        uint8_t *stackdir, uint8_t *direction, float *weight,
                        ptrdiff_t *weightscan, uint8_t *visited,
                        ptrdiff_t edge_count, ptrdiff_t dims[2], int order) {
  ptrdiff_t e[2][8] = {{0, -1, -1, -1, 0, 1, 1, 1},
                       {1, 1, 0, -1, -1, -1, 0, 1}};

  for (ptrdiff_t j = 0; j < dims[1]; j++) {
    for (ptrdiff_t i = 0; i < dims[0]; i++) {
      visited[j * dims[0] + i] = 0;
    }
  }

  ptrdiff_t stack_top = 0;
  ptrdiff_t next_node = dims[0] * dims[1];

  for (ptrdiff_t j = 0; j < dims[1]; j++) {
    for (ptrdiff_t i = 0; i < dims[0]; i++) {
      ptrdiff_t idx = j * dims[0] + i;

      if (visited[idx]) continue;

      visited[idx] = 1;

      stack[stack_top] = idx;
      stackdir[stack_top++] = direction[idx];

      while (stack_top > 0) {
        while (stackdir[stack_top - 1]) {
          // We still have neighbors left to check
          ptrdiff_t node = stack[stack_top - 1];
          ptrdiff_t iu = node % dims[0];
          ptrdiff_t ju = node / dims[0];

          uint8_t n = next_neighbor(stackdir[stack_top - 1]);
          stackdir[stack_top - 1] ^= (1 << n);

          ptrdiff_t in = iu + e[order & 1][n];
          ptrdiff_t jn = ju + e[(order ^ 1) & 1][n];

          if (!visited[jn * dims[0] + in]) {
            visited[jn * dims[0] + in] = 1;

            // Push the neighbor and its outgoing edges on the stack
            stack[stack_top] = jn * dims[0] + in;
            stackdir[stack_top++] = direction[jn * dims[0] + in];
          } else if (visited[jn * dims[0] + in] == 1) {
            // This is a cycle. The neighbor is already on the stack
            assert(0 && "flow_routing_tsort detected a cycle. This is a bug.");
          }
        }
        // Pop the stack.
        ptrdiff_t node = stack[--stack_top];
        ptrdiff_t iu = node % dims[0];
        ptrdiff_t ju = node / dims[0];

        visited[node] = 2;

        if (next_node <= 0) {
          assert(0 && "flow_routing_tsort ran out of nodes. This is a bug.");
        }
        stream[--next_node] = node;

        ptrdiff_t offset = weightscan[node];

        // Add its edges to source and target
        for (int n = 0; n < 8; n++) {
          if (direction[node] & (1 << n)) {
            ptrdiff_t in = iu + e[order & 1][n];
            ptrdiff_t jn = ju + e[(order ^ 1) & 1][n];

            if (edge_count <= 0) {
              assert(0 &&
                     "flow_routing_tsort ran out of nodes. This is a bug.");
            }
            source[--edge_count] = node;
            target[edge_count] = jn * dims[0] + in;

            sorted_weight[edge_count] = weight[offset++];
          }
        }
      }
    }
  }
}

static uint8_t get_direction(ptrdiff_t i, ptrdiff_t j, int order) {
  uint8_t e[2][9] = {{8, 16, 32, 4, 0, 64, 2, 1, 128},
                     {8, 4, 2, 16, 0, 1, 32, 64, 128}};

  return e[order][(j + 1) * 3 + (i + 1)];
}

TOPOTOOLBOX_API
void resolve_flats_shortest_path(uint8_t *direction, uint8_t *resolved,
                                 float *path_distance, ptrdiff_t *heap,
                                 ptrdiff_t *back, float *dem, ptrdiff_t dims[2],
                                 int order) {
  ptrdiff_t e[2][8] = {{0, -1, -1, -1, 0, 1, 1, 1},
                       {1, 1, 0, -1, -1, -1, 0, 1}};
  float chamfer[8] = {1.0f, SQRT2, 1.0f, SQRT2, 1.0f, SQRT2, 1.0f, SQRT2};

  PriorityQueue pq = pq_create(dims[0] * dims[1], heap, back, path_distance, 0);

  // Initialize distances based on the resolved array.
  for (ptrdiff_t j = 0; j < dims[1]; j++) {
    for (ptrdiff_t i = 0; i < dims[0]; i++) {
      ptrdiff_t idx = j * dims[0] + i;
      if (resolved[idx] == 0) {
        pq_insert(&pq, idx, INFINITY);
      } else {
        path_distance[idx] = -INFINITY;
      }
    }
  }
  // Identify presills, flats that border a resolved pixel with the same
  // elevation
  for (ptrdiff_t j = 0; j < dims[1]; j++) {
    for (ptrdiff_t i = 0; i < dims[0]; i++) {
      ptrdiff_t idx = j * dims[0] + i;
      if (resolved[idx]) {
        continue;
      }
      float z = dem[idx];
      for (ptrdiff_t n = 0; n < 8; n++) {
        ptrdiff_t in = i + e[order & 1][n];
        ptrdiff_t jn = j + e[(order ^ 1) & 1][n];

        if (in < 0 || in >= dims[0] || jn < 0 || jn >= dims[1]) {
          // Pixel is on the boundary, call it a presill
          pq_decrease_key(&pq, idx, 0.0);
          break;
        }

        ptrdiff_t idxn = jn * dims[0] + in;

        if (resolved[idxn] && z == dem[idxn]) {
          // These are the presills, so their distance is 0
          pq_decrease_key(&pq, idx, 0.0);

          // Flow the presill towards the neighbor
          direction[idx] = get_direction(in - i, jn - j, order);
          break;
        }
      }
    }
  }

  while (!pq_isempty(&pq)) {
    ptrdiff_t idx = pq_deletemin(&pq);
    float d = pq_get_priority(&pq, idx);

    ptrdiff_t j = idx / dims[0];
    ptrdiff_t i = idx % dims[0];

    float z = dem[idx];

    for (ptrdiff_t n = 0; n < 8; n++) {
      ptrdiff_t in = i + e[order & 1][n];
      ptrdiff_t jn = j + e[(order ^ 1) & 1][n];

      // Skip any neighboring pixels that are outside the domain, are
      // resolved or are not equal in elevation. The weird elevation
      // comparison is in case the neighboring pixel is a NaN, in
      // which case !(neighbor_z <= z) is true and we skip the pixel.
      if (in < 0 || in >= dims[0] || jn < 0 || jn >= dims[1] ||
          resolved[jn * dims[0] + in] || !(dem[jn * dims[0] + in] == z)) {
        continue;
      }
      float proposal = d + chamfer[n];

      if (proposal < pq_get_priority(&pq, jn * dims[0] + in)) {
        pq_decrease_key(&pq, jn * dims[0] + in, proposal);

        // Save the flow direction to the current pixel as the neighbor's flow
        // direction
        direction[jn * dims[0] + in] = get_direction(i - in, j - jn, order);
      }
    }
  }
}

// Weights are 1 for shortest path flat resolution
TOPOTOOLBOX_API
void resolve_flats_shortest_path_weights(float *weight, ptrdiff_t count) {
  for (ptrdiff_t e = 0; e < count; e++) {
    weight[e] = 1.0;
  }
}
