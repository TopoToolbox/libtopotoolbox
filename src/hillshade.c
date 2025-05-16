#define TOPOTOOLBOX_BUILD

#include <math.h>

#include "topotoolbox.h"

#define PI 3.14159265358979323846f
#define PI_2 1.57079632679489661923f

// Compute hillshade value from gradient and illumination vector
static inline float compute_hillshade(float dx, float dy, float sx, float sy,
                                      float sz) {
  float inorm = 1.0f / sqrtf(dx * dx + dy * dy + 1.0f);

  float nx = -dx * inorm;
  float ny = -dy * inorm;
  float nz = inorm;

  return nx * sx + ny * sy + nz * sz;
}

TOPOTOOLBOX_API
void gradient_secondorder(float *restrict p0, float *restrict p1, float *dem,
                          float cellsize, ptrdiff_t dims[2]) {
  // Compute second order
  ptrdiff_t m = dims[0];
  ptrdiff_t n = dims[1];

  // Boundaries use a stencil of [-3 4 -1]/(2*cs) and [1 -4 3]/(2*cs)
  // Interior points use a stencil of [-1 0 1]/(2*cs)

  // Compute gradient in the first dimension
  for (ptrdiff_t j = 0; j < n; j++) {
    ptrdiff_t i = 0;

    p0[j * m + i] =
        -dem[j * m + i + 2] + 4 * dem[j * m + i + 1] - 3 * dem[j * m + i];
    p0[j * m + i] /= 2 * cellsize;

    for (i = 1; i < m - 1; i++) {
      p0[j * m + i] =
          (dem[j * m + i + 1] - dem[j * m + i - 1]) / (2 * cellsize);
    }

    i = m - 1;
    p0[j * m + i] =
        dem[j * m + i - 2] - 4 * dem[j * m + i - 1] + 3 * dem[j * m + i];
    p0[j * m + i] /= 2 * cellsize;
  }

  // Compute gradient in the second dimension
  ptrdiff_t j = 0;
  for (ptrdiff_t i = 0; i < m; i++) {
    p1[j * m + i] = -dem[(j + 2) * m + i] + 4 * dem[(j + 1) * m + i] -
                    3 * dem[(j + 0) * m + i];
    p1[j * m + i] /= 2 * cellsize;
  }

  for (j = 1; j < n - 1; j++) {
    for (ptrdiff_t i = 0; i < m; i++) {
      p1[j * m + i] =
          (dem[(j + 1) * m + i] - dem[(j - 1) * m + i]) / (2 * cellsize);
    }
  }

  j = n - 1;
  for (ptrdiff_t i = 0; i < m; i++) {
    p1[j * m + i] =
        dem[(j - 2) * m + i] - 4 * dem[(j - 1) * m + i] + 3 * dem[j * m + i];
    p1[j * m + i] /= 2 * cellsize;
  }
}

TOPOTOOLBOX_API
void hillshade(float *output, float *dx, float *dy, float *dem, float azimuth,
               float altitude, float cellsize, ptrdiff_t dims[2]) {
  float sx = sinf(PI_2 - altitude) * cosf(azimuth);
  float sy = sinf(PI_2 - altitude) * sinf(azimuth);
  float sz = cosf(PI_2 - altitude);

  gradient_secondorder(dx, dy, dem, cellsize, dims);

  for (ptrdiff_t j = 0; j < dims[1]; j++) {
    for (ptrdiff_t i = 0; i < dims[0]; i++) {
      float dxij = dx[j * dims[0] + i];
      float dyij = dy[j * dims[0] + i];

      output[j * dims[0] + i] = compute_hillshade(dxij, dyij, sx, sy, sz);
    }
  }
}

////////////////////
// A lower memory version of hillshade fuses the loops together so
// intermediate arrays can be avoided

TOPOTOOLBOX_API
void hillshade_fused(float *output, float *dem, float azimuth, float altitude,
                     float cellsize, ptrdiff_t dims[2]) {
  float sx = sinf(PI_2 - altitude) * cosf(azimuth);
  float sy = sinf(PI_2 - altitude) * sinf(azimuth);
  float sz = cosf(PI_2 - altitude);

  ptrdiff_t m = dims[0];
  ptrdiff_t n = dims[1];

  ptrdiff_t i, j;

  j = 0;
  {
    i = 0;
    float dx =
        (-dem[j * m + i + 2] + 4 * dem[j * m + i + 1] - 3 * dem[j * m + i]) /
        (2 * cellsize);
    float dy = (-dem[(j + 2) * m + i] + 4 * dem[(j + 1) * m + i] -
                3 * dem[(j + 0) * m + i]) /
               (2 * cellsize);
    output[j * m + i] = compute_hillshade(dx, dy, sx, sy, sz);
  }
  for (i = 1; i < (m - 1); i++) {
    float dx = (dem[j * m + i + 1] - dem[j * m + i - 1]) / (2 * cellsize);
    float dy = (-dem[(j + 2) * m + i] + 4 * dem[(j + 1) * m + i] -
                3 * dem[(j + 0) * m + i]) /
               (2 * cellsize);

    output[j * m + i] = compute_hillshade(dx, dy, sx, sy, sz);
  }
  {
    i = m - 1;
    float dx =
        (dem[j * m + i - 2] - 4 * dem[j * m + i - 1] + 3 * dem[j * m + i]) /
        (2 * cellsize);
    float dy = (-dem[(j + 2) * m + i] + 4 * dem[(j + 1) * m + i] -
                3 * dem[(j + 0) * m + i]) /
               (2 * cellsize);
    output[j * m + i] = compute_hillshade(dx, dy, sx, sy, sz);
  }

  for (j = 1; j < (n - 1); j++) {
    {
      i = 0;
      float dx =
          (-dem[j * m + i + 2] + 4 * dem[j * m + i + 1] - 3 * dem[j * m + i]) /
          (2 * cellsize);
      float dy = (dem[(j + 1) * m + i] - dem[(j - 1) * m + i]) / (2 * cellsize);
      output[j * m + i] = compute_hillshade(dx, dy, sx, sy, sz);
    }
    for (i = 1; i < (m - 1); i++) {
      float dx = (dem[j * m + i + 1] - dem[j * m + i - 1]) / (2 * cellsize);
      float dy = (dem[(j + 1) * m + i] - dem[(j - 1) * m + i]) / (2 * cellsize);

      output[j * m + i] = compute_hillshade(dx, dy, sx, sy, sz);
    }
    {
      i = m - 1;
      float dx =
          (dem[j * m + i - 2] - 4 * dem[j * m + i - 1] + 3 * dem[j * m + i]) /
          (2 * cellsize);
      float dy = (dem[(j + 1) * m + i] - dem[(j - 1) * m + i]) / (2 * cellsize);
      output[j * m + i] = compute_hillshade(dx, dy, sx, sy, sz);
    }
  }

  j = n - 1;
  {
    i = 0;
    float dx =
        (-dem[j * m + i + 2] + 4 * dem[j * m + i + 1] - 3 * dem[j * m + i]) /
        (2 * cellsize);
    float dy =
        (dem[(j - 2) * m + i] - 4 * dem[(j - 1) * m + i] + 3 * dem[j * m + i]) /
        (2 * cellsize);
    output[j * m + i] = compute_hillshade(dx, dy, sx, sy, sz);
  }
  for (i = 1; i < (m - 1); i++) {
    float dx = (dem[j * m + i + 1] - dem[j * m + i - 1]) / (2 * cellsize);
    float dy =
        (dem[(j - 2) * m + i] - 4 * dem[(j - 1) * m + i] + 3 * dem[j * m + i]) /
        (2 * cellsize);

    output[j * m + i] = compute_hillshade(dx, dy, sx, sy, sz);
  }
  {
    i = m - 1;
    float dx =
        (dem[j * m + i - 2] - 4 * dem[j * m + i - 1] + 3 * dem[j * m + i]) /
        (2 * cellsize);
    float dy =
        (dem[(j - 2) * m + i] - 4 * dem[(j - 1) * m + i] + 3 * dem[j * m + i]) /
        (2 * cellsize);
    output[j * m + i] = compute_hillshade(dx, dy, sx, sy, sz);
  }
}
