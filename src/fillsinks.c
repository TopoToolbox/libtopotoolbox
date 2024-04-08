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
void fillsinks(float *output, float *dem, ptrdiff_t nrows, ptrdiff_t ncols) {
  for (ptrdiff_t col = 0; col < ncols; col++) {
    for (ptrdiff_t row = 0; row < nrows; row++) {
      // Invert the DEM
      dem[col * nrows + row] *= -1.0;

      // Set the boundary pixels of the output equal to the DEM and
      // the interior pixels equal to -INFINITY.
      if ((row == 0 || row == (nrows - 1)) ||
          (col == 0 || col == (ncols - 1))) {
        output[col * nrows + row] = dem[col * nrows + row];
      } else {
        output[col * nrows + row] = -INFINITY;
      }
    }
  }

  reconstruct(output, dem, nrows, ncols);

  // Revert the DEM and the output
  for (ptrdiff_t col = 0; col < ncols; col++) {
    for (ptrdiff_t row = 0; row < nrows; row++) {
      dem[col * nrows + row] *= -1.0;
      output[col * nrows + row] *= -1.0;
    }
  }
}
