#define TOPOTOOLBOX_BUILD

#include <math.h>
#include <stddef.h>
#include <stdint.h>

#include "morphology/reconstruct.h"
#include "topotoolbox.h"

/*
  Fill sinks in a digital elevation model.

  The arrays pointed to by `output` and `dem` should represent a two
  dimensional array of size (nrows, ncols). The filled DEM is written
  to the `output` array.

  `dem` is modified during the processing, but the modifications are
  reverted before the function returns. Be careful when accessing
  `dem` concurrently with `fillsinks`.

  Sinks are filled using grayscale morphological reconstruction.
*/
TOPOTOOLBOX_API
void fillsinks(float *output, float *dem, ptrdiff_t dims[2]) {
  for (ptrdiff_t j = 0; j < dims[1]; j++) {
    for (ptrdiff_t i = 0; i < dims[0]; i++) {
      // Invert the DEM
      dem[j * dims[0] + i] *= -1.0;

      // Set the boundary pixels of the output equal to the DEM and
      // the interior pixels equal to -INFINITY.
      if ((i == 0 || i == (dims[0] - 1)) || (j == 0 || j == (dims[1] - 1))) {
        output[j * dims[0] + i] = dem[j * dims[0] + i];
      } else {
        output[j * dims[0] + i] = -INFINITY;
      }
    }
  }

  reconstruct(output, dem, dims);

  // Revert the DEM and the output
  for (ptrdiff_t j = 0; j < dims[1]; j++) {
    for (ptrdiff_t i = 0; i < dims[0]; i++) {
      dem[j * dims[0] + i] *= -1.0;
      output[j * dims[0] + i] *= -1.0;
    }
  }
}

TOPOTOOLBOX_API
void fillsinks_hybrid(float *output, ptrdiff_t *queue, float *dem,
                      ptrdiff_t dims[2]) {
  for (ptrdiff_t j = 0; j < dims[1]; j++) {
    for (ptrdiff_t i = 0; i < dims[0]; i++) {
      // Invert the DEM
      dem[j * dims[0] + i] *= -1.0;

      // Set the boundary pixels of the output equal to the DEM and
      // the interior pixels equal to -INFINITY.
      if ((i == 0 || i == (dims[0] - 1)) || (j == 0 || j == (dims[1] - 1))) {
        output[j * dims[0] + i] = dem[j * dims[0] + i];
      } else {
        output[j * dims[0] + i] = -INFINITY;
      }
    }
  }

  reconstruct_hybrid(output, queue, dem, dims);

  // Revert the DEM and the output
  for (ptrdiff_t j = 0; j < dims[1]; j++) {
    for (ptrdiff_t i = 0; i < dims[0]; i++) {
      dem[j * dims[0] + i] *= -1.0;
      output[j * dims[0] + i] *= -1.0;
    }
  }
}
