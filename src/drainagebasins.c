#define TOPOTOOLBOX_BUILD

#include "topotoolbox.h"

TOPOTOOLBOX_API
void drainagebasins(ptrdiff_t *basins, ptrdiff_t *source, ptrdiff_t *target,
                    ptrdiff_t dims[2]) {
  // Initialize the basins array
  for (ptrdiff_t j = 0; j < dims[1]; j++) {
    for (ptrdiff_t i = 0; i < dims[0]; i++) {
      basins[j * dims[0] + i] = 0;
    }
  }

  ptrdiff_t basin_count = 1;  // Start the basin labels at one.
  for (ptrdiff_t j = dims[1] - 1; j >= 0; j--) {
    for (ptrdiff_t i = dims[0] - 1; i >= 0; i--) {
      ptrdiff_t src = source[j * dims[0] + i];
      ptrdiff_t tgt = target[j * dims[0] + i];

      if (tgt < 0) {
        // src has no downstream neighbors and is the root of the
        // drainage basin tree
        basins[src] = basin_count++;
      } else {
        basins[src] = basins[tgt];
      }
    }
  }
}
