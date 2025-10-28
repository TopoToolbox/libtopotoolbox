#include <math.h>
#include <stddef.h>
#include <stdint.h>

#include "topotoolbox.h"

TOPOTOOLBOX_API
void erode(float* restrict output, float* restrict dem,
           uint8_t* restrict structuring_element, ptrdiff_t io_dims[2],
           ptrdiff_t se_dims[3]) {
  /* COMPUTE SE CENTER (BIAS TOWARDS TL); SAME AS (se_dims[0] - 1) / 2 */
  ptrdiff_t slow_dim_center = (se_dims[1] + 1) / 2 - 1;
  ptrdiff_t fast_dim_center = (se_dims[0] + 1) / 2 - 1;

  /* GRAYSCALE EROSION */
  for (ptrdiff_t slow_dim_idx = 0; slow_dim_idx < io_dims[1]; slow_dim_idx++) {
    for (ptrdiff_t fast_dim_idx = 0; fast_dim_idx < io_dims[0];
         fast_dim_idx++) {
      if (isnan(dem[fast_dim_idx + slow_dim_idx * io_dims[0]])) {
        output[fast_dim_idx + slow_dim_idx * io_dims[0]] = NAN;
        continue;
      }

      output[fast_dim_idx + slow_dim_idx * io_dims[0]] = INFINITY;

      for (ptrdiff_t se_slice = 0; se_slice < se_dims[2]; se_slice++) {
        for (ptrdiff_t se_slow = 0; se_slow < se_dims[1]; se_slow++) {
          // (se_slow - slow_dim_center): offset in input array to SE cell
          // on input cell; same for other dimensions except slowest
          ptrdiff_t new_slow_dim_index =
              slow_dim_idx + (se_slow - slow_dim_center);

          if (new_slow_dim_index < 0 || new_slow_dim_index >= io_dims[1])
            continue;

          for (ptrdiff_t se_fast = 0; se_fast < se_dims[0]; se_fast++) {
            if (structuring_element[se_fast + se_slow * se_dims[0] +
                                    se_slice * se_dims[0] * se_dims[1]] == 0)
              continue;

            ptrdiff_t new_fast_dim_index =
                fast_dim_idx + (se_fast - fast_dim_center);

            if (new_fast_dim_index < 0 || new_fast_dim_index >= io_dims[0])
              continue;

            if (isnan(
                    dem[new_fast_dim_index + new_slow_dim_index * io_dims[0]]))
              continue;

            if (output[fast_dim_idx + slow_dim_idx * io_dims[0]] <
                dem[new_fast_dim_index + new_slow_dim_index * io_dims[0]])
              continue;

            output[fast_dim_idx + slow_dim_idx * io_dims[0]] =
                dem[new_fast_dim_index + new_slow_dim_index * io_dims[0]];
          }
        }
      }
    }
  }
  return;
}
