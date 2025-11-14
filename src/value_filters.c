#define TOPOTOOLBOX_BUILD

#include <math.h>
#include <stddef.h>
#include <stdint.h>

#include "topotoolbox.h"

TOPOTOOLBOX_API
void min_filter(float* restrict output, float* restrict dem,
                uint8_t* restrict structuring_element, ptrdiff_t io_dims[2],
                ptrdiff_t se_dims[3]) {
  /* COMPUTE SE CENTER (BIAS TOWARDS TL); SAME AS (se_dims[0] - 1) / 2 */
  ptrdiff_t slow_dim_center = (se_dims[1] + 1) / 2 - 1;
  ptrdiff_t fast_dim_center = (se_dims[0] + 1) / 2 - 1;

  for (ptrdiff_t slow_dim_idx = 0; slow_dim_idx < io_dims[1]; slow_dim_idx++) {
    for (ptrdiff_t fast_dim_idx = 0; fast_dim_idx < io_dims[0];
         fast_dim_idx++) {
      ptrdiff_t original_index = fast_dim_idx + slow_dim_idx * io_dims[0];

      if (isnan(dem[original_index])) {
        output[original_index] = NAN;
        continue;
      }

      output[original_index] = INFINITY;

      for (ptrdiff_t se_slice = 0; se_slice < se_dims[2]; se_slice++) {
        for (ptrdiff_t se_slow = 0; se_slow < se_dims[1]; se_slow++) {
          // (se_slow - slow_dim_center): offset in input array to SE cell
          // on input cell; same for other dimensions except slowest
          ptrdiff_t slow_offset = slow_dim_idx + (se_slow - slow_dim_center);

          if (slow_offset < 0 || slow_offset >= io_dims[1]) continue;

          for (ptrdiff_t se_fast = 0; se_fast < se_dims[0]; se_fast++) {
            if (structuring_element[se_fast + se_slow * se_dims[0] +
                                    se_slice * se_dims[0] * se_dims[1]] == 0)
              continue;

            ptrdiff_t fast_offset = fast_dim_idx + (se_fast - fast_dim_center);

            ptrdiff_t offset_index = fast_offset + slow_offset * io_dims[0];

            if (fast_offset < 0 || fast_offset >= io_dims[0]) continue;

            if (isnan(dem[offset_index])) continue;

            if (output[original_index] < dem[offset_index]) continue;

            output[original_index] = dem[offset_index];
          }
        }
      }
    }
  }
  return;
}

TOPOTOOLBOX_API
void max_filter(float* restrict output, float* restrict dem,
                uint8_t* restrict structuring_element, ptrdiff_t io_dims[2],
                ptrdiff_t se_dims[3]) {
  /* COMPUTE SE CENTER (BIAS TOWARDS TL); SAME AS (se_dims[0] - 1) / 2 */
  ptrdiff_t slow_dim_center = (se_dims[1] + 1) / 2 - 1;
  ptrdiff_t fast_dim_center = (se_dims[0] + 1) / 2 - 1;

  for (ptrdiff_t slow_dim_idx = 0; slow_dim_idx < io_dims[1]; slow_dim_idx++) {
    for (ptrdiff_t fast_dim_idx = 0; fast_dim_idx < io_dims[0];
         fast_dim_idx++) {
      ptrdiff_t original_index = fast_dim_idx + slow_dim_idx * io_dims[0];
      if (isnan(dem[original_index])) {
        output[original_index] = NAN;
        continue;
      }

      output[original_index] = -INFINITY;

      for (ptrdiff_t se_slice = 0; se_slice < se_dims[2]; se_slice++) {
        for (ptrdiff_t se_slow = 0; se_slow < se_dims[1]; se_slow++) {
          // (se_slow - slow_dim_center): offset in input array to SE cell
          // on input cell; same for other dimensions except slowest
          ptrdiff_t slow_offset = slow_dim_idx + (se_slow - slow_dim_center);

          if (slow_offset < 0 || slow_offset >= io_dims[1]) continue;

          for (ptrdiff_t se_fast = 0; se_fast < se_dims[0]; se_fast++) {
            if (structuring_element[se_fast + se_slow * se_dims[0] +
                                    se_slice * se_dims[0] * se_dims[1]] == 0)
              continue;

            ptrdiff_t fast_offset = fast_dim_idx + (se_fast - fast_dim_center);

            ptrdiff_t offset_index = fast_offset + slow_offset * io_dims[0];

            if (fast_offset < 0 || fast_offset >= io_dims[0]) continue;

            if (isnan(dem[offset_index])) continue;

            if (output[original_index] > dem[offset_index]) continue;

            output[original_index] = dem[offset_index];
          }
        }
      }
    }
  }
  return;
}

TOPOTOOLBOX_API
void min_filter_square(float* restrict output, float* restrict dem,
                       float* restrict tmp, uint8_t width,
                       ptrdiff_t io_dims[2]) {
  /* COMPUTE SE CENTER (BIAS TOWARDS TL); SAME AS (width - 1) / 2 */
  ptrdiff_t se_center = (width + 1) / 2 - 1;

  /* CHECK ALONG FIRST DIMENSION */
  for (ptrdiff_t second_dim_idx = 0; second_dim_idx < io_dims[1];
       second_dim_idx++) {
    for (ptrdiff_t first_dim_idx = 0; first_dim_idx < io_dims[0];
         first_dim_idx++) {
      ptrdiff_t original_index = first_dim_idx + second_dim_idx * io_dims[0];

      // no skipping of NANs in temporary storage

      tmp[original_index] = INFINITY;

      for (ptrdiff_t se_idx = 0; se_idx < width; se_idx++) {
        ptrdiff_t axis_offset = first_dim_idx + (se_idx - se_center);

        ptrdiff_t offset_index = axis_offset + second_dim_idx * io_dims[0];

        if (axis_offset < 0 || axis_offset >= io_dims[0]) continue;

        if (isnan(dem[offset_index])) continue;

        if (tmp[original_index] < dem[offset_index]) continue;

        tmp[original_index] = dem[offset_index];
      }
    }
  }

  /* CHECK ALONG SECOND DIMENSION */
  for (ptrdiff_t second_dim_idx = 0; second_dim_idx < io_dims[1];
       second_dim_idx++) {
    for (ptrdiff_t first_dim_idx = 0; first_dim_idx < io_dims[0];
         first_dim_idx++) {
      ptrdiff_t original_index = first_dim_idx + second_dim_idx * io_dims[0];

      // copy over NANs from input to output
      if (isnan(dem[original_index])) {
        output[original_index] = NAN;
        continue;
      }

      output[original_index] = INFINITY;

      for (ptrdiff_t se_idx = 0; se_idx < width; se_idx++) {
        ptrdiff_t axis_offset = second_dim_idx + (se_idx - se_center);

        ptrdiff_t offset_index = first_dim_idx + axis_offset * io_dims[0];

        if (axis_offset < 0 || axis_offset >= io_dims[1]) continue;

        if (isnan(tmp[offset_index])) continue;

        if (output[original_index] < tmp[offset_index]) continue;

        output[original_index] = tmp[offset_index];
      }
    }
  }

  return;
}

TOPOTOOLBOX_API
void max_filter_square(float* restrict output, float* restrict dem,
                       float* restrict tmp, uint8_t width,
                       ptrdiff_t io_dims[2]) {
  /* COMPUTE SE CENTER (BIAS TOWARDS TL); SAME AS (width - 1) / 2 */
  ptrdiff_t se_center = (width + 1) / 2 - 1;

  /* CHECK ALONG FIRST DIMENSION */
  for (ptrdiff_t second_dim_idx = 0; second_dim_idx < io_dims[1];
       second_dim_idx++) {
    for (ptrdiff_t first_dim_idx = 0; first_dim_idx < io_dims[0];
         first_dim_idx++) {
      ptrdiff_t original_index = first_dim_idx + second_dim_idx * io_dims[0];

      // no skipping of NANs in temporary storage

      tmp[original_index] = -INFINITY;

      for (ptrdiff_t se_idx = 0; se_idx < width; se_idx++) {
        ptrdiff_t axis_offset = first_dim_idx + (se_idx - se_center);

        ptrdiff_t offset_index = axis_offset + second_dim_idx * io_dims[0];

        if (axis_offset < 0 || axis_offset >= io_dims[0]) continue;

        if (isnan(dem[offset_index])) continue;

        if (tmp[original_index] > dem[offset_index]) continue;

        tmp[original_index] = dem[offset_index];
      }
    }
  }

  /* CHECK ALONG SECOND DIMENSION */
  for (ptrdiff_t second_dim_idx = 0; second_dim_idx < io_dims[1];
       second_dim_idx++) {
    for (ptrdiff_t first_dim_idx = 0; first_dim_idx < io_dims[0];
         first_dim_idx++) {
      ptrdiff_t original_index = first_dim_idx + second_dim_idx * io_dims[0];

      // copy over NANs from input to output
      if (isnan(dem[original_index])) {
        output[original_index] = NAN;
        continue;
      }

      output[original_index] = -INFINITY;

      for (ptrdiff_t se_idx = 0; se_idx < width; se_idx++) {
        ptrdiff_t axis_offset = second_dim_idx + (se_idx - se_center);

        ptrdiff_t offset_index = first_dim_idx + axis_offset * io_dims[0];

        if (axis_offset < 0 || axis_offset >= io_dims[1]) continue;

        if (isnan(tmp[offset_index])) continue;

        if (output[original_index] > tmp[offset_index]) continue;

        output[original_index] = tmp[offset_index];
      }
    }
  }

  return;
}
