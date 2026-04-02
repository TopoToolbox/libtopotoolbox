#include "grid.h"

#include <assert.h>
#include <gdal.h>
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>

#include "../include/topotoolbox.h"
#include "utils.h"

void GridFree(Grid *g) {
  free(g->data);
  g->data = NULL;
  g->dims[0] = 0;
  g->dims[0] = 0;
  g->cellsize = 0.0;
}

Grid GridCreate(GridDataType dtype, void *buffer, float cellsize,
                ptrdiff_t dims[2]) {
  Grid res = {0};

  res.cellsize = cellsize;
  res.dtype = dtype;
  res.dims[0] = dims[0];
  res.dims[1] = dims[1];

  if (res.dims[0] > 0 && (res.dims[1] > PTRDIFF_MAX / res.dims[0])) {
    return res;
  }

  if (buffer == NULL) {
    buffer = calloc(GridElementCount(res), GridElementSize(res));
  }
  res.data = buffer;

  return res;
}

Grid GridFromFile(const char *path) {
  Grid res = {0};

  const GDALAccess eAccess = GA_ReadOnly;

  GDALDatasetH hDataset = GDALOpen(path, eAccess);

  if (hDataset == NULL) {
    return res;
  }

  GDALRasterBandH hBand = GDALGetRasterBand(hDataset, 1);

  ptrdiff_t nx = GDALGetRasterBandXSize(hBand);
  ptrdiff_t ny = GDALGetRasterBandYSize(hBand);

  res.dims[0] = nx;
  res.dims[1] = ny;

  if (nx > 0 && (ny > PTRDIFF_MAX / nx)) {
    return res;
  }

  double adfTransform[6];
  GDALGetGeoTransform(hDataset, adfTransform);

  if (fabs(fabs(adfTransform[1]) - fabs(adfTransform[5])) > 1e-6) {
    return res;
  }

  res.cellsize = adfTransform[1];
  res.dtype = GridF32;

  float *buffer = calloc(GridElementCount(res), GridElementSize(res));

  if (buffer == NULL) {
    return res;
  }

  CPLErr err = GDALRasterIO(hBand, GF_Read, 0, 0, nx, ny, buffer, nx, ny,
                            GDT_Float32, 0, 0);
  if (err) {
    return res;
  }

  res.data = buffer;

  GDALClose(hDataset);

  return res;
}

GDALDataType GridGDALDataType(Grid g) {
  switch (g.dtype) {
    case GridF32:
      return GDT_Float32;
    case GridU8:
      return GDT_Byte;
    case GridU32:
      return GDT_UInt32;
    case GridI32:
      return GDT_Int32;
    case GridIdx:
      return GDT_Int64;
    default:
      return GDT_Unknown;
  }
}

int GridToFile(Grid g, const char *dst) {
  GDALDriverH driver = GDALGetDriverByName("GTiff");
  if (driver == NULL) return -1;

  GDALDataType dtype = GridGDALDataType(g);
  GDALDatasetH ds =
      GDALCreate(driver, dst, g.dims[0], g.dims[1], 1, dtype, NULL);

  if (ds == NULL) return -1;

  double gt[6] = {0, g.cellsize, 0, 0, 0, -g.cellsize};
  if (GDALSetGeoTransform(ds, gt)) return -1;

  GDALRasterBandH band = GDALGetRasterBand(ds, 1);

  if (band == NULL) return -1;

  CPLErr err = GDALRasterIO(band, GF_Write, 0, 0, g.dims[0], g.dims[1], g.data,
                            g.dims[0], g.dims[1], dtype, 0, 0);
  if (err) return -1;

  GDALClose(ds);

  return 0;
}

Grid GridRandom(uint64_t seed, float cellsize, ptrdiff_t dims[2]) {
  Grid res = {0};

  ptrdiff_t nx = dims[0];
  ptrdiff_t ny = dims[1];

  res.dims[0] = nx;
  res.dims[1] = ny;

  if (nx > 0 && (ny > PTRDIFF_MAX / nx)) {
    return res;
  }

  res.dtype = GridF32;

  float *buffer = calloc(GridElementCount(res), GridElementSize(res));

  if (buffer == NULL) {
    return res;
  }

  uint32_t c = seed >> 32;
  uint32_t d = seed & ((1UL << 32) - 1);
  for (ptrdiff_t j = 0; j < dims[1]; j++) {
    for (ptrdiff_t i = 0; i < dims[0]; i++) {
      buffer[j * dims[0] + i] = 100.0f * pcg4d(i, j, c, d);
    }
  }

  res.data = buffer;

  return res;
}

size_t GridElementSize(Grid g) {
  switch (g.dtype) {
    case GridF32:
      return sizeof(float);
    case GridU8:
      return sizeof(uint8_t);
    case GridU32:
      return sizeof(uint32_t);
    case GridI32:
      return sizeof(int32_t);
    case GridIdx:
      return sizeof(ptrdiff_t);
    default:
      return 0;
  }
}

ptrdiff_t GridElementCount(Grid g) { return g.dims[0] * g.dims[1]; }

size_t GridSize(Grid g) { return GridElementSize(g) * GridElementCount(g); }

Grid GridCopy(Grid g) {
  Grid res = {0};

  res.dims[0] = g.dims[0];
  res.dims[1] = g.dims[1];

  res.cellsize = g.cellsize;

  void *buffer = calloc(GridElementCount(g), GridElementSize(g));

  if (buffer == NULL) {
    return res;
  }

  res.data = buffer;
  res.dtype = g.dtype;

  memcpy(res.data, g.data, GridSize(g));

  return res;
}

Grid GridZeros(float cellsize, ptrdiff_t dims[2]) {
  Grid res = {0};

  res.cellsize = cellsize;
  res.dtype = GridF32;
  res.dims[0] = dims[0];
  res.dims[1] = dims[1];

  if (res.dims[0] > 0 && (res.dims[1] > PTRDIFF_MAX / res.dims[0])) {
    return res;
  }

  void *buffer = calloc(GridElementCount(res), GridElementSize(res));

  res.data = buffer;

  return res;
}

Grid GridZerosLike(Grid g) {
  Grid res = {0};

  res.dims[0] = g.dims[0];
  res.dims[1] = g.dims[1];

  res.cellsize = g.cellsize;

  res.dtype = g.dtype;

  res.data = calloc(GridElementCount(g), GridElementSize(g));

  return res;
}

Grid GridFill(float value, float cellsize, ptrdiff_t dims[2]) {
  Grid res = {0};

  res.cellsize = cellsize;
  res.dtype = GridF32;
  res.dims[0] = dims[0];
  res.dims[1] = dims[1];

  if (res.dims[0] > 0 && (res.dims[1] > PTRDIFF_MAX / res.dims[0])) {
    return res;
  }

  float *buffer = calloc(GridElementCount(res), GridElementSize(res));

  if (buffer == NULL) return res;

  for (ptrdiff_t j = 0; j < dims[1]; j++) {
    for (ptrdiff_t i = 0; i < dims[0]; i++) {
      buffer[j * dims[0] + i] = value;
    }
  }

  res.data = buffer;

  return res;
}

Grid GridTranspose(Grid g) {
  Grid res = {0};
  res.cellsize = g.cellsize;
  res.dims[0] = g.dims[1];
  res.dims[1] = g.dims[0];
  res.dtype = g.dtype;

  float *buffer = calloc(GridElementCount(res), GridElementSize(res));

  if (buffer == NULL) return res;

  for (ptrdiff_t j = 0; j < g.dims[1]; j++) {
    for (ptrdiff_t i = 0; i < g.dims[0]; i++) {
      buffer[i * res.dims[0] + j] = ((float *)g.data)[j * g.dims[0] + i];
    }
  }

  res.data = buffer;

  return res;
}

Grid GridSubtract(Grid g1, Grid g2) {
  assert(g1.dims[0] == g2.dims[0]);
  assert(g1.dims[1] == g2.dims[1]);

  // Arithmetic operations are only defined for float Grids
  assert(g1.dtype == g2.dtype);
  assert(g1.dtype == GridF32);

  Grid res = GridZerosLike(g1);

  if (g1.data && g2.data && res.data) {
    float *r = res.data;
    float *d1 = g1.data;
    float *d2 = g2.data;

    for (ptrdiff_t j = 0; j < g1.dims[1]; j++) {
      for (ptrdiff_t i = 0; i < g1.dims[0]; i++) {
        r[j * g1.dims[0] + i] = d1[j * g1.dims[0] + i] - d2[j * g1.dims[0] + i];
      }
    }
  }

  return res;
}

float GridMin(Grid g) {
  assert(g.dtype == GridF32);
  float *data = (float *)g.data;

  float min = INFINITY;
  if (data) {
    for (ptrdiff_t j = 0; j < g.dims[1]; j++) {
      for (ptrdiff_t i = 0; i < g.dims[0]; i++) {
        min = fminf(min, data[j * g.dims[0] + i]);
      }
    }
  }
  return min;
}

float GridMax(Grid g) {
  assert(g.dtype == GridF32);
  float *data = (float *)g.data;

  float max = -INFINITY;
  if (data) {
    for (ptrdiff_t j = 0; j < g.dims[1]; j++) {
      for (ptrdiff_t i = 0; i < g.dims[0]; i++) {
        max = fmaxf(max, data[j * g.dims[0] + i]);
      }
    }
  }
  return max;
}

int GridAllNonnegative(Grid g) {
  assert(g.dtype == GridF32);
  float *data = (float *)g.data;

  if (!data) return 1;

  for (ptrdiff_t j = 0; j < g.dims[1]; j++) {
    for (ptrdiff_t i = 0; i < g.dims[0]; i++) {
      if (data[j * g.dims[0] + i] < 0) {
        return 0;
      }
    }
  }
  return 1;
}

int GridAll(Grid g) {
  assert(g.dtype == GridF32);
  float *data = (float *)g.data;

  if (!data) return 0;

  for (ptrdiff_t j = 0; j < g.dims[1]; j++) {
    for (ptrdiff_t i = 0; i < g.dims[0]; i++) {
      if (data[j * g.dims[0] + i] == 0.0) {
        return 0;
      }
    }
  }
  return 1;
}

int GridAny(Grid g) {
  assert(g.dtype == GridF32);
  float *data = (float *)g.data;

  if (!data) return 0;

  for (ptrdiff_t j = 0; j < g.dims[1]; j++) {
    for (ptrdiff_t i = 0; i < g.dims[0]; i++) {
      if (data[j * g.dims[0] + i] != 0.0) {
        return 1;
      }
    }
  }
  return 0;
}

int GridEq(Grid g1, Grid g2) {
  assert(g1.dtype == GridF32);
  assert(g2.dtype == GridF32);

  if (g1.dims[0] != g2.dims[0] || g1.dims[1] != g2.dims[1]) {
    return 1;
  }

  float *data1 = g1.data;
  float *data2 = g2.data;

  if (!data1 || !data2) return 0;

  for (ptrdiff_t j = 0; j < g1.dims[1]; j++) {
    for (ptrdiff_t i = 0; i < g2.dims[0]; i++) {
      if (data1[j * g1.dims[0] + i] != data2[j * g1.dims[0] + i]) {
        return 0;
      }
    }
  }
  return 1;
}

int GridAllClose(Grid g1, Grid g2, double rtol) {
  assert(g1.dtype == GridF32);
  assert(g2.dtype == GridF32);

  if (g1.dims[0] != g2.dims[0] || g1.dims[1] != g2.dims[1]) {
    return 1;
  }

  float *data1 = g1.data;
  float *data2 = g2.data;

  if (!data1 || !data2) return 0;

  for (ptrdiff_t j = 0; j < g1.dims[1]; j++) {
    for (ptrdiff_t i = 0; i < g2.dims[0]; i++) {
      float z1 = data1[j * g1.dims[0] + i];
      float z2 = data2[j * g1.dims[0] + i];
      if (fabsf(z1 - z2) > fabsf(z2) * rtol) {
        return 0;
      }
    }
  }

  return 1;
}

Grid GridInvert(Grid g) {
  Grid res = GridCopy(g);

  assert(res.dtype == GridF32);

  if (res.data) {
    float *data = (float *)res.data;
    for (ptrdiff_t j = 0; j < g.dims[1]; j++) {
      for (ptrdiff_t i = 0; i < g.dims[0]; i++) {
        data[j * g.dims[0] + i] *= -1;
      }
    }
  }
  return res;
}

Grid GridNot(Grid g) {
  Grid res = GridCopy(g);

  assert(res.dtype == GridF32);

  if (res.data) {
    float *data = (float *)res.data;
    for (ptrdiff_t j = 0; j < g.dims[1]; j++) {
      for (ptrdiff_t i = 0; i < g.dims[0]; i++) {
        data[j * g.dims[0] + i] = 1 - data[j * g.dims[0] + i];
      }
    }
  }
  return res;
}

Grid GridIsNan(Grid g) {
  Grid res = GridCopy(g);

  if (res.data) {
    float *data = res.data;
    for (ptrdiff_t j = 0; j < g.dims[1]; j++) {
      for (ptrdiff_t i = 0; i < g.dims[0]; i++) {
        data[j * g.dims[0] + i] = isnan(data[j * g.dims[0] + i]);
      }
    }
  }

  return res;
}

Grid GridFillsinks(Grid g) {
  assert(g.dtype == GridF32);
  float *data = g.data;

  Grid g2 = GridZerosLike(g);

  if (g2.data) {
    Grid bc = GridCreate(GridU8, NULL, g.cellsize, g.dims);

    if (bc.data) {
      for (ptrdiff_t j = 0; j < g.dims[1]; j++) {
        for (ptrdiff_t i = 0; i < g.dims[0]; i++) {
          if (i == 0 || i == g.dims[0] - 1 || j == 0 || j == g.dims[1] - 1) {
            ((uint8_t *)bc.data)[j * g.dims[0] + i] = 1;
          }
          if (isnan(data[j * g.dims[0] + i])) {
            ((uint8_t *)bc.data)[j * g.dims[0] + i] = 1;
            data[j * g.dims[0] + i] = -INFINITY;
          }
        }
      }

      Grid queue = GridCreate(GridIdx, NULL, g.cellsize, g.dims);
      if (queue.data) {
        fillsinks_hybrid(g2.data, queue.data, data, bc.data, g.dims);

        // Restore NaNs
        for (ptrdiff_t j = 0; j < g.dims[1]; j++) {
          for (ptrdiff_t i = 0; i < g.dims[0]; i++) {
            if (data[j * g.dims[0] + i] == -INFINITY) {
              data[j * g.dims[0] + i] = NAN;
              ((float *)g2.data)[j * g.dims[0] + i] = NAN;
            }
          }
        }

        GridFree(&queue);
      }
      GridFree(&bc);
    }
  }
  return g2;
}

Grid GridLabelSinks(Grid g) {
  assert(g.dtype == GridF32);

  Grid res = GridZerosLike(g);

  if (!res.data || !g.data) return res;

  float *r = res.data;
  float *data = g.data;

  ptrdiff_t j_offset[8] = {-1, -1, -1, 0, 1, 1, 1, 0};
  ptrdiff_t i_offset[8] = {1, 0, -1, -1, -1, 0, 1, 1};

  for (ptrdiff_t j = 1; j < g.dims[1] - 1; j++) {
    for (ptrdiff_t i = 1; i < g.dims[0] - 1; i++) {
      float z = data[i + g.dims[0] * j];
      ptrdiff_t up_neighbor_count = 0;
      for (int32_t neighbor = 0; neighbor < 8; neighbor++) {
        ptrdiff_t neighbor_i = i + i_offset[neighbor];
        ptrdiff_t neighbor_j = j + j_offset[neighbor];

        if (data[g.dims[0] * neighbor_j + neighbor_i] > z) {
          up_neighbor_count++;
        }
      }
      r[j * g.dims[0] + i] = up_neighbor_count == 8 ? 1.0 : 0.0;
    }
  }

  return res;
}

Grid GridIdentifyFlats(Grid g) {
  assert(g.dtype == GridF32);

  Grid res = GridCreate(GridI32, NULL, g.cellsize, g.dims);

  if (!g.data || !res.data) return res;

  identifyflats(res.data, g.data, g.dims);

  return res;
}

Grid GridGWDTComputeCosts(Grid flats, Grid dem, Grid demf) {
  assert(flats.dtype == GridI32);
  assert(dem.dtype == GridF32);
  assert(dem.dtype == GridF32);

  Grid costs = GridZerosLike(dem);
  if (!flats.data || !dem.data || !costs.data) return costs;

  Grid conncomps = GridCreate(GridIdx, NULL, dem.cellsize, dem.dims);
  if (!conncomps.data) return costs;

  gwdt_computecosts(costs.data, conncomps.data, flats.data, dem.data, demf.data,
                    dem.dims);

  GridFree(&conncomps);

  return costs;
}

Grid GridGWDT(Grid costs, Grid flats) {
  assert(costs.dtype == GridF32);
  assert(flats.dtype == GridI32);

  Grid dist = GridZerosLike(costs);
  if (dist.data && costs.data && flats.data) {
    Grid prev = GridCreate(GridIdx, NULL, dist.cellsize, dist.dims);
    if (prev.data) {
      Grid heap = GridCreate(GridIdx, NULL, dist.cellsize, dist.dims);

      if (heap.data) {
        Grid back = GridCreate(GridIdx, NULL, dist.cellsize, dist.dims);

        if (back.data) {
          gwdt(dist.data, prev.data, costs.data, flats.data, heap.data,
               back.data, costs.dims);

          GridFree(&back);
        }
        GridFree(&heap);
      }
      GridFree(&prev);
    }
  }

  return dist;
}

Grid GridGradient8(Grid dem, int use_mp) {
  Grid res = GridZerosLike(dem);
  if (res.data) {
    gradient8(res.data, dem.data, dem.cellsize, use_mp, dem.dims);
  }
  return res;
}
