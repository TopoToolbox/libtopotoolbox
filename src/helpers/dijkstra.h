#pragma once

#include <stddef.h>
#include <stdint.h>

#include "priority_queue.h"

/*
  Generic 8-connected grid Dijkstra
  ===================================

  Caller provides a cost function and optionally an on-update callback.
  Common cost functions (e.g. D8 step) are provided here; domain-specific
  ones (e.g. perpendicular track distance) live in their own modules.

  Usage pattern:
    GridDijkstra gd;
    grid_dijkstra_init(&gd, dist_array, dims);
    grid_dijkstra_seed(&gd, idx, priority, on_update, ctx);  // repeat as needed
    grid_dijkstra_run(&gd, mask, cost_fn, on_update, ctx);
    grid_dijkstra_free(&gd);
*/

// Cost function: returns candidate priority for neighbor pixel 'to'.
//   from      - settled source pixel index
//   to        - neighbor pixel index
//   from_dist - current priority of 'from'
//   dir       - neighbor direction 0-7 (0,2,5,7 diagonal; 1,3,4,6 cardinal)
//   ctx       - user context
// Return FLT_MAX to skip the neighbor.
typedef float (*DijkstraCostFn)(ptrdiff_t from, ptrdiff_t to, float from_dist,
                                int dir, void *ctx);

// Update callback: called when pixel 'idx' is assigned or updated with
// 'new_dist', both during seeding and expansion. May be NULL.
typedef void (*DijkstraUpdateFn)(ptrdiff_t idx, float new_dist, void *ctx);

typedef struct {
  float *dist;         // priority array (size dims[0]*dims[1]), caller-owned
  ptrdiff_t *pq_heap;  // internally allocated
  ptrdiff_t *pq_back;  // internally allocated
  PriorityQueue q;
  ptrdiff_t dims[2];
} GridDijkstra;

// Allocate internal PQ arrays and initialise dist to FLT_MAX.
void grid_dijkstra_init(GridDijkstra *gd, float *dist, ptrdiff_t dims[2]);

// Insert or update a seed pixel.
// Only accepted if priority < current dist[idx]. Handles duplicates.
void grid_dijkstra_seed(GridDijkstra *gd, ptrdiff_t idx, float priority,
                        DijkstraUpdateFn on_update, void *ctx);

// Run expansion until the queue is empty.
//   mask      - skip neighbor if mask[nidx] == 0; NULL = all valid
//   cost      - cost function
//   on_update - called on insert or priority update; NULL to skip
//   ctx       - passed to cost and on_update
void grid_dijkstra_run(GridDijkstra *gd, const int8_t *mask,
                       DijkstraCostFn cost, DijkstraUpdateFn on_update,
                       void *ctx);

// Free internally allocated PQ arrays.
void grid_dijkstra_free(GridDijkstra *gd);

// Pre-coded cost: accumulated D8 step distance.
// sqrt(2) for diagonal neighbors, 1.0 for cardinal. ctx: unused.
float dijkstra_cost_d8(ptrdiff_t from, ptrdiff_t to, float from_dist, int dir,
                       void *ctx);
