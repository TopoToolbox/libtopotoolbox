#ifndef TT3_TEST_GRID_H
#define TT3_TEST_GRID_H

#include <stddef.h>
#include <stdint.h>

typedef enum GridDataType {
  GridF32,
  GridU8,
  GridU32,
  GridI32,
  GridIdx
} GridDataType;

typedef struct Grid {
  void *data;
  ptrdiff_t dims[2];
  float cellsize;
  GridDataType dtype;
} Grid;

void GridFree(Grid *g);

// Creation
Grid GridCreate(GridDataType, void *buffer, float cellsize, ptrdiff_t dims[2]);
Grid GridFromFile(const char *path);
Grid GridRandom(uint64_t seed, float cellsize, ptrdiff_t dims[2]);
Grid GridCopy(Grid g);
Grid GridZeros(float cellsize, ptrdiff_t dims[2]);
Grid GridZerosLike(Grid g);
Grid GridFill(float value, float cellsize, ptrdiff_t dims[2]);
Grid GridTranspose(Grid g);

// Output
int GridToFile(Grid g, const char *dst);

// Details
size_t GridElementSize(Grid g);
ptrdiff_t GridElementCount(Grid g);
size_t GridSize(Grid g);

// Arithmetic
Grid GridSubtract(Grid g1, Grid g2);
int GridAll(Grid g);
int GridAny(Grid g);
int GridAllNonnegative(Grid g);
int GridEq(Grid g1, Grid g2);
int GridAllClose(Grid g1, Grid g2, double rtol);
Grid GridInvert(Grid g);
Grid GridNot(Grid g);
Grid GridIsNan(Grid g);

float GridMin(Grid g);
float GridMax(Grid g);

// TopoToolbox
Grid GridFillsinks(Grid g);
Grid GridLabelSinks(Grid g);
Grid GridIdentifyFlats(Grid g);
Grid GridGWDTComputeCosts(Grid flats, Grid dem, Grid demf);
Grid GridGWDT(Grid costs, Grid flats);
Grid GridGradient8(Grid dem, int use_mp);

#endif
