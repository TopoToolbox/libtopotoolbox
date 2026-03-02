#include "dijkstra.h"

#include <float.h>
#include <stdlib.h>

static const int k_di8[8] = {-1, -1, -1, 0, 0, 1, 1, 1};
static const int k_dj8[8] = {-1, 0, 1, -1, 1, -1, 0, 1};

float dijkstra_cost_d8(ptrdiff_t from, ptrdiff_t to, float from_dist, int dir,
                       void *ctx) {
  (void)from;
  (void)to;
  (void)ctx;
  static const float dd8[8] = {1.41421356f, 1.0f, 1.41421356f, 1.0f, 1.0f,
                               1.41421356f, 1.0f, 1.41421356f};
  return from_dist + dd8[dir];
}

void grid_dijkstra_init(GridDijkstra *gd, float *dist, ptrdiff_t dims[2]) {
  ptrdiff_t total = dims[0] * dims[1];
  gd->dist = dist;
  gd->dims[0] = dims[0];
  gd->dims[1] = dims[1];
  for (ptrdiff_t i = 0; i < total; i++) dist[i] = FLT_MAX;
  gd->pq_heap = (ptrdiff_t *)malloc(total * sizeof(ptrdiff_t));
  gd->pq_back = (ptrdiff_t *)malloc(total * sizeof(ptrdiff_t));
  for (ptrdiff_t i = 0; i < total; i++) gd->pq_back[i] = -1;
  gd->q = pq_create(total, gd->pq_heap, gd->pq_back, dist, 0);
}

void grid_dijkstra_seed(GridDijkstra *gd, ptrdiff_t idx, float priority,
                        DijkstraUpdateFn on_update, void *ctx) {
  if (priority < gd->dist[idx]) {
    if (on_update) on_update(idx, priority, ctx);
    if (gd->q.back[idx] >= 0)
      pq_decrease_key(&gd->q, idx, priority);
    else
      pq_insert(&gd->q, idx, priority);
  }
}

void grid_dijkstra_run(GridDijkstra *gd, const int8_t *mask,
                       DijkstraCostFn cost, DijkstraUpdateFn on_update,
                       void *ctx) {
  while (!pq_isempty(&gd->q)) {
    ptrdiff_t idx = pq_deletemin(&gd->q);
    ptrdiff_t ci = idx % gd->dims[0];
    ptrdiff_t cj = idx / gd->dims[0];

    for (int n = 0; n < 8; n++) {
      ptrdiff_t ni = ci + k_di8[n];
      ptrdiff_t nj = cj + k_dj8[n];
      if (ni < 0 || ni >= gd->dims[0] || nj < 0 || nj >= gd->dims[1]) continue;
      ptrdiff_t nidx = nj * gd->dims[0] + ni;
      if (mask && !mask[nidx]) continue;

      float new_dist = cost(idx, nidx, gd->dist[idx], n, ctx);
      if (new_dist < gd->dist[nidx]) {
        if (gd->q.back[nidx] >= 0) {
          if (on_update) on_update(nidx, new_dist, ctx);
          pq_decrease_key(&gd->q, nidx, new_dist);
        } else if (gd->dist[nidx] >= FLT_MAX) {
          if (on_update) on_update(nidx, new_dist, ctx);
          pq_insert(&gd->q, nidx, new_dist);
        }
        // settled pixels (back < 0, dist < FLT_MAX) are skipped
      }
    }
  }
}

void grid_dijkstra_free(GridDijkstra *gd) {
  free(gd->pq_heap);
  free(gd->pq_back);
  gd->pq_heap = NULL;
  gd->pq_back = NULL;
}
