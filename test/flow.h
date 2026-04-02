#ifndef TT3_TEST_FLOW_H
#define TT3_TEST_FLOW_H

#include <stddef.h>

#include "grid.h"

typedef struct Flow {
  ptrdiff_t *source;
  ptrdiff_t *target;
  ptrdiff_t *stream;
  float *fraction;
  ptrdiff_t count;
  ptrdiff_t dims[2];
  float cellsize;
} Flow;

// Destructor
void FlowFree(Flow *f);

// Creation
Flow FlowInit(Grid g);

// Flow routing
Flow FlowRouteD8Carve(Grid g);

// TopoToolbox
Grid FlowAccumulation(Flow f);
Grid FlowDrainageBasins(Flow f);

#endif
