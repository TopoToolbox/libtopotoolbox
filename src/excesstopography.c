#define TOPOTOOLBOX_BUILD

#include <math.h>
#include <stddef.h>
#include <stdint.h>

#include "priority_queue.h"
#include "topotoolbox.h"

/*
  Propose an updated value for pixel (i,j) based on the upwind
  discretization of the eikonal equation.

  fi is cellsize/f[i,j] or cellsize * threshold_slope[i,j]
 */
static float eikonal_solver(float *solution, float fi, ptrdiff_t i, ptrdiff_t j,
                            ptrdiff_t dims[2]) {
  // Find the vertical neighbors.
  //
  // Points outside the DEM are given infinite values. Since the
  // minimum value of these is taken, gradients at boundary pixels are
  // always computed from their neighboring interior pixels.
  float south = i < dims[0] - 1 ? solution[j * dims[0] + i + 1] : INFINITY;
  float north = i > 0 ? solution[j * dims[0] + i - 1] : INFINITY;

  float u1 = fminf(south, north);

  // Find the horizontal neighbors
  float east = j < dims[1] - 1 ? solution[(j + 1) * dims[0] + i] : INFINITY;
  float west = j > 0 ? solution[(j - 1) * dims[0] + i] : INFINITY;

  float u2 = fminf(east, west);

  // Solve the discretized eikonal equation for the central pixel
  // using its two upwind derivatives.
  if (fabsf(u1 - u2) < fi) {
    return (u1 + u2) / 2 + sqrtf(-(u1 - u2) * (u1 - u2) + 2 * (fi * fi)) / 2;
  } else {
    return fminf(u1 + fi, u2 + fi);
  }
}

/*
  Compute the two dimensional excess topography by solving the eikonal
  equation with the fast sweeping method.
 */
TOPOTOOLBOX_API
void excesstopography_fsm2d(float *excess, float *dem, float *threshold_slopes,
                            float cellsize, ptrdiff_t dims[2]) {
  // Initialize excess == dem
  for (ptrdiff_t j = 0; j < dims[1]; j++) {
    for (ptrdiff_t i = 0; i < dims[0]; i++) {
      excess[j * dims[0] + i] = dem[j * dims[0] + i];
    }
  }
  ptrdiff_t count = dims[0] * dims[1];

  while (count > 0) {
    count = 0;

    // Perform four eikonal_solver sweeps in alternating directions

    // TODO(wsk): Do we still have to skip the boundary pixels?

    // Sweep 1
    for (ptrdiff_t j = 1; j < dims[1] - 1; j++) {
      for (ptrdiff_t i = 1; i < dims[0] - 1; i++) {
        float fi = cellsize * threshold_slopes[j * dims[0] + i];
        float proposal = eikonal_solver(excess, fi, i, j, dims);
        if (proposal < excess[j * dims[0] + i]) {
          excess[j * dims[0] + i] = proposal;
          count += 1;
        }
      }
    }

    // Sweep 2
    for (ptrdiff_t j = dims[1] - 2; j > 0; j--) {
      for (ptrdiff_t i = 1; i < dims[0] - 1; i++) {
        float fi = cellsize * threshold_slopes[j * dims[0] + i];
        float proposal = eikonal_solver(excess, fi, i, j, dims);
        if (proposal < excess[j * dims[0] + i]) {
          excess[j * dims[0] + i] = proposal;
          count += 1;
        }
      }
    }

    // Sweep 3
    for (ptrdiff_t j = dims[1] - 2; j > 0; j--) {
      for (ptrdiff_t i = dims[0] - 2; i > 0; i--) {
        float fi = cellsize * threshold_slopes[j * dims[0] + i];
        float proposal = eikonal_solver(excess, fi, i, j, dims);
        if (proposal < excess[j * dims[0] + i]) {
          excess[j * dims[0] + i] = proposal;
          count += 1;
        }
      }
    }

    // Sweep 4
    for (ptrdiff_t j = 1; j < dims[1] - 1; j++) {
      for (ptrdiff_t i = dims[0] - 2; i > 0; i--) {
        float fi = cellsize * threshold_slopes[j * dims[0] + i];
        float proposal = eikonal_solver(excess, fi, i, j, dims);
        if (proposal < excess[j * dims[0] + i]) {
          excess[j * dims[0] + i] = proposal;
          count += 1;
        }
      }
    }
  }
}

/*
  Compute the two-dimensional excess topography by solving the eikonal
  equation with the fast marching method.
 */
TOPOTOOLBOX_API
void excesstopography_fmm2d(float *excess, ptrdiff_t *heap, ptrdiff_t *back,
                            float *dem, float *threshold_slopes, float cellsize,
                            ptrdiff_t dims[2]) {
  // Initialize the arrays
  // Pixels start with the elevation given by the DEM
  for (ptrdiff_t j = 0; j < dims[1]; j++) {
    for (ptrdiff_t i = 0; i < dims[0]; i++) {
      ptrdiff_t idx = j * dims[0] + i;
      back[idx] = idx;
      excess[idx] = dem[idx];
      heap[idx] = idx;
    }
  }

  // Initialize a priority queue.
  // See priority_queue.h for more details
  ptrdiff_t count = dims[1] * dims[0];
  PriorityQueue q = pq_create(count, heap, back, excess, 1);

  while (!pq_isempty(&q)) {
    ptrdiff_t trial = pq_deletemin(&q);
    float trial_elevation = pq_get_priority(&q, trial);

    ptrdiff_t j = trial / dims[0];
    ptrdiff_t i = trial % dims[0];

    // South neighbor
    if (i < dims[0] - 1 && excess[j * dims[0] + i + 1] >= trial_elevation) {
      float fi = cellsize * threshold_slopes[j * dims[0] + i + 1];
      float proposal = eikonal_solver(q.priorities, fi, i + 1, j, dims);
      pq_decrease_key(&q, j * dims[0] + i + 1, proposal);
    }

    // North neighbor
    if (i > 0 && excess[j * dims[0] + i - 1] >= trial_elevation) {
      float fi = cellsize * threshold_slopes[j * dims[0] + i - 1];
      float proposal = eikonal_solver(q.priorities, fi, i - 1, j, dims);
      pq_decrease_key(&q, j * dims[0] + i - 1, proposal);
    }

    // East neighbor
    if (j < dims[1] - 1 && excess[(j + 1) * dims[0] + i] >= trial_elevation) {
      float fi = cellsize * threshold_slopes[(j + 1) * dims[0] + i];
      float proposal = eikonal_solver(q.priorities, fi, i, j + 1, dims);
      pq_decrease_key(&q, (j + 1) * dims[0] + i, proposal);
    }

    // West neighbor
    if (j > 0 && excess[(j - 1) * dims[0] + i] >= trial_elevation) {
      float fi = cellsize * threshold_slopes[(j - 1) * dims[0] + i];
      float proposal = eikonal_solver(q.priorities, fi, i, j - 1, dims);
      pq_decrease_key(&q, (j - 1) * dims[0] + i, proposal);
    }
  }
}

/*
  Compute the excess topography using three-dimensional lithological
  variability and the fast marching method.
 */
TOPOTOOLBOX_API
void excesstopography_fmm3d(float *excess, ptrdiff_t *heap, ptrdiff_t *back,
                            float *dem, float *lithstack,
                            float *threshold_slopes, float cellsize,
                            ptrdiff_t dims[2], ptrdiff_t nlayers) {
  // Initialize the arrays
  // Pixels start with the elevation given by the DEM
  for (ptrdiff_t j = 0; j < dims[1]; j++) {
    for (ptrdiff_t i = 0; i < dims[0]; i++) {
      ptrdiff_t idx = j * dims[0] + i;
      back[idx] = idx;
      excess[idx] = dem[idx];
      heap[idx] = idx;
    }
  }
  ptrdiff_t count = dims[1] * dims[0];
  PriorityQueue q = pq_create(count, heap, back, excess, 1);

  while (!pq_isempty(&q)) {
    ptrdiff_t trial = pq_deletemin(&q);
    float trial_elevation = pq_get_priority(&q, trial);

    ptrdiff_t j = trial / dims[0];
    ptrdiff_t i = trial % dims[0];

    // South neighbor
    if (i < dims[0] - 1 &&
        pq_get_priority(&q, j * dims[0] + i + 1) >= trial_elevation) {
      float proposal = pq_get_priority(&q, j * dims[0] + i + 1);

      // Choose the lowest layer that produces a proposal solution
      // that is below the upper surface of that layer.
      for (int layer = 0; layer < nlayers; layer++) {
        float fi = cellsize * threshold_slopes[layer];
        proposal = eikonal_solver(q.priorities, fi, i + 1, j, dims);
        if (proposal < lithstack[(j * dims[0] + i + 1) * nlayers + layer]) {
          pq_decrease_key(&q, j * dims[0] + i + 1, proposal);
          break;
        }
      }
    }

    // North neighbor
    if (i > 0 && pq_get_priority(&q, j * dims[0] + i - 1) >= trial_elevation) {
      float proposal = pq_get_priority(&q, j * dims[0] + i - 1);
      for (int layer = 0; layer < nlayers; layer++) {
        float fi = cellsize * threshold_slopes[layer];
        proposal = eikonal_solver(q.priorities, fi, i - 1, j, dims);
        if (proposal < lithstack[(j * dims[0] + i - 1) * nlayers + layer]) {
          pq_decrease_key(&q, j * dims[0] + i - 1, proposal);
          break;
        }
      }
    }

    // East neighbor
    if (j < dims[1] - 1 &&
        pq_get_priority(&q, (j + 1) * dims[0] + i) >= trial_elevation) {
      float proposal = pq_get_priority(&q, (j + 1) * dims[0] + i);

      for (int layer = 0; layer < nlayers; layer++) {
        float fi = cellsize * threshold_slopes[layer];
        proposal = eikonal_solver(q.priorities, fi, i, j + 1, dims);
        if (proposal < lithstack[((j + 1) * dims[0] + i) * nlayers + layer]) {
          pq_decrease_key(&q, (j + 1) * dims[0] + i, proposal);
          break;
        }
      }
    }

    // West neighbor
    if (j > 0 &&
        pq_get_priority(&q, (j - 1) * dims[0] + i) >= trial_elevation) {
      float proposal = pq_get_priority(&q, (j - 1) * dims[0] + i);
      for (int layer = 0; layer < nlayers; layer++) {
        float fi = cellsize * threshold_slopes[layer];
        proposal = eikonal_solver(q.priorities, fi, i, j - 1, dims);
        if (proposal < lithstack[((j - 1) * dims[0] + i) * nlayers + layer]) {
          pq_decrease_key(&q, (j - 1) * dims[0] + i, proposal);
          break;
        }
      }
    }
  }
}
