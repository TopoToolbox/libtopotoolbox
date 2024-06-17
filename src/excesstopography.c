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
                            ptrdiff_t nrows, ptrdiff_t ncols) {
  // Find the vertical neighbors.
  //
  // Points outside the DEM are given infinite values. Since the
  // minimum value of these is taken, gradients at boundary pixels are
  // always computed from their neighboring interior pixels.
  float south = i < nrows - 1 ? solution[j * nrows + i + 1] : INFINITY;
  float north = i > 0 ? solution[j * nrows + i - 1] : INFINITY;

  float u1 = fminf(south, north);

  // Find the horizontal neighbors
  float east = j < ncols - 1 ? solution[(j + 1) * nrows + i] : INFINITY;
  float west = j > 0 ? solution[(j - 1) * nrows + i] : INFINITY;

  float u2 = fminf(east, west);

  // Solve the discretized eikonal equation for the central pixel
  // using its two upwind derivatives.
  if (fabsf(u1 - u2) < fi) {
    return (u1 + u2) / 2 +
           sqrtf((u1 + u2) * (u1 + u2) - 2 * (u1 * u1 + u2 * u2 - fi * fi)) / 2;
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
                            float cellsize, ptrdiff_t nrows, ptrdiff_t ncols) {
  // Initialize excess == dem
  for (ptrdiff_t col = 0; col < ncols; col++) {
    for (ptrdiff_t row = 0; row < nrows; row++) {
      excess[col * nrows + row] = dem[col * nrows + row];
    }
  }
  ptrdiff_t count = nrows * ncols;

  while (count > 0) {
    count = 0;

    // Perform four eikonal_solver sweeps in alternating directions

    // TODO(wsk): Do we still have to skip the boundary pixels?

    // Sweep 1
    for (ptrdiff_t col = 1; col < ncols - 1; col++) {
      for (ptrdiff_t row = 1; row < nrows - 1; row++) {
        float fi = cellsize * threshold_slopes[col * nrows + row];
        float proposal = eikonal_solver(excess, fi, row, col, nrows, ncols);
        if (proposal < excess[col * nrows + row]) {
          excess[col * nrows + row] = proposal;
          count += 1;
        }
      }
    }

    // Sweep 2
    for (ptrdiff_t col = ncols - 2; col > 0; col--) {
      for (ptrdiff_t row = 1; row < nrows - 1; row++) {
        float fi = cellsize * threshold_slopes[col * nrows + row];
        float proposal = eikonal_solver(excess, fi, row, col, nrows, ncols);
        if (proposal < excess[col * nrows + row]) {
          excess[col * nrows + row] = proposal;
          count += 1;
        }
      }
    }

    // Sweep 3
    for (ptrdiff_t col = ncols - 2; col > 0; col--) {
      for (ptrdiff_t row = nrows - 2; row > 0; row--) {
        float fi = cellsize * threshold_slopes[col * nrows + row];
        float proposal = eikonal_solver(excess, fi, row, col, nrows, ncols);
        if (proposal < excess[col * nrows + row]) {
          excess[col * nrows + row] = proposal;
          count += 1;
        }
      }
    }

    // Sweep 4
    for (ptrdiff_t col = 1; col < ncols - 1; col++) {
      for (ptrdiff_t row = nrows - 2; row > 0; row--) {
        float fi = cellsize * threshold_slopes[col * nrows + row];
        float proposal = eikonal_solver(excess, fi, row, col, nrows, ncols);
        if (proposal < excess[col * nrows + row]) {
          excess[col * nrows + row] = proposal;
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
                            ptrdiff_t nrows, ptrdiff_t ncols) {
  // Initialize the arrays
  // Pixels start with the elevation given by the DEM
  for (ptrdiff_t j = 0; j < ncols; j++) {
    for (ptrdiff_t i = 0; i < nrows; i++) {
      ptrdiff_t idx = j * nrows + i;
      back[idx] = idx;
      excess[idx] = dem[idx];
      heap[idx] = idx;
    }
  }

  // Initialize a priority queue.
  // See priority_queue.h for more details
  ptrdiff_t count = ncols * nrows;
  PriorityQueue q = pq_create(count, heap, back, excess, 1);

  while (!pq_isempty(&q)) {
    ptrdiff_t trial = pq_deletemin(&q);
    float trial_elevation = pq_get_priority(&q, trial);

    ptrdiff_t col = trial / nrows;
    ptrdiff_t row = trial % nrows;

    // South neighbor
    if (row < nrows - 1 && excess[col * nrows + row + 1] >= trial_elevation) {
      float fi = cellsize * threshold_slopes[col * nrows + row + 1];
      float proposal =
          eikonal_solver(q.priorities, fi, row + 1, col, nrows, ncols);
      pq_decrease_key(&q, col * nrows + row + 1, proposal);
    }

    // North neighbor
    if (row > 0 && excess[col * nrows + row - 1] >= trial_elevation) {
      float fi = cellsize * threshold_slopes[col * nrows + row - 1];
      float proposal =
          eikonal_solver(q.priorities, fi, row - 1, col, nrows, ncols);
      pq_decrease_key(&q, col * nrows + row - 1, proposal);
    }

    // East neighbor
    if (col < ncols - 1 && excess[(col + 1) * nrows + row] >= trial_elevation) {
      float fi = cellsize * threshold_slopes[(col + 1) * nrows + row];
      float proposal =
          eikonal_solver(q.priorities, fi, row, col + 1, nrows, ncols);
      pq_decrease_key(&q, (col + 1) * nrows + row, proposal);
    }

    // West neighbor
    if (col > 0 && excess[(col - 1) * nrows + row] >= trial_elevation) {
      float fi = cellsize * threshold_slopes[(col - 1) * nrows + row];
      float proposal =
          eikonal_solver(q.priorities, fi, row, col - 1, nrows, ncols);
      pq_decrease_key(&q, (col - 1) * nrows + row, proposal);
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
                            ptrdiff_t nrows, ptrdiff_t ncols,
                            ptrdiff_t nlayers) {
  // Initialize the arrays
  // Pixels start with the elevation given by the DEM
  for (ptrdiff_t j = 0; j < ncols; j++) {
    for (ptrdiff_t i = 0; i < nrows; i++) {
      ptrdiff_t idx = j * nrows + i;
      back[idx] = idx;
      excess[idx] = dem[idx];
      heap[idx] = idx;
    }
  }
  ptrdiff_t count = ncols * nrows;
  PriorityQueue q = pq_create(count, heap, back, excess, 1);

  while (!pq_isempty(&q)) {
    ptrdiff_t trial = pq_deletemin(&q);
    float trial_elevation = pq_get_priority(&q, trial);

    ptrdiff_t col = trial / nrows;
    ptrdiff_t row = trial % nrows;

    // South neighbor
    if (row < nrows - 1 &&
        pq_get_priority(&q, col * nrows + row + 1) >= trial_elevation) {
      float proposal = pq_get_priority(&q, col * nrows + row + 1);

      // Choose the lowest layer that produces a proposal solution
      // that is below the upper surface of that layer.
      for (int layer = 0; layer < nlayers; layer++) {
        float fi = cellsize * threshold_slopes[layer];
        proposal = eikonal_solver(q.priorities, fi, row + 1, col, nrows, ncols);
        if (proposal < lithstack[(col * nrows + row + 1) * nlayers + layer]) {
          pq_decrease_key(&q, col * nrows + row + 1, proposal);
          break;
        }
      }
    }

    // North neighbor
    if (row > 0 &&
        pq_get_priority(&q, col * nrows + row - 1) >= trial_elevation) {
      float proposal = pq_get_priority(&q, col * nrows + row - 1);
      for (int layer = 0; layer < nlayers; layer++) {
        float fi = cellsize * threshold_slopes[layer];
        proposal = eikonal_solver(q.priorities, fi, row - 1, col, nrows, ncols);
        if (proposal < lithstack[(col * nrows + row - 1) * nlayers + layer]) {
          pq_decrease_key(&q, col * nrows + row - 1, proposal);
          break;
        }
      }
    }

    // East neighbor
    if (col < ncols - 1 &&
        pq_get_priority(&q, (col + 1) * nrows + row) >= trial_elevation) {
      float proposal = pq_get_priority(&q, (col + 1) * nrows + row);

      for (int layer = 0; layer < nlayers; layer++) {
        float fi = cellsize * threshold_slopes[layer];
        proposal = eikonal_solver(q.priorities, fi, row, col + 1, nrows, ncols);
        if (proposal < lithstack[((col + 1) * nrows + row) * nlayers + layer]) {
          pq_decrease_key(&q, (col + 1) * nrows + row, proposal);
          break;
        }
      }
    }

    // West neighbor
    if (col > 0 &&
        pq_get_priority(&q, (col - 1) * nrows + row) >= trial_elevation) {
      float proposal = pq_get_priority(&q, (col - 1) * nrows + row);
      for (int layer = 0; layer < nlayers; layer++) {
        float fi = cellsize * threshold_slopes[layer];
        proposal = eikonal_solver(q.priorities, fi, row, col - 1, nrows, ncols);
        if (proposal < lithstack[((col - 1) * nrows + row) * nlayers + layer]) {
          pq_decrease_key(&q, (col - 1) * nrows + row, proposal);
          break;
        }
      }
    }
  }
}
