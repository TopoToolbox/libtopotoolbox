#include "flow.h"

#include <stdlib.h>

#include "../include/topotoolbox.h"
#include "grid.h"

void FlowFree(Flow *f) {
  free(f->source);
  free(f->target);
  free(f->stream);
  free(f->fraction);

  f->source = NULL;
  f->target = NULL;
  f->stream = NULL;
  f->fraction = NULL;
  f->dims[0] = 0;
  f->dims[1] = 0;
  f->count = 0;
  f->cellsize = 0.0;
}

Flow FlowInit(Grid dem) {
  Flow fd = {0};
  fd.cellsize = dem.cellsize;
  fd.dims[0] = dem.dims[0];
  fd.dims[1] = dem.dims[1];

  fd.stream = calloc(GridElementCount(dem), sizeof *fd.stream);
  fd.source = calloc(GridElementCount(dem), sizeof *fd.source);
  fd.target = calloc(GridElementCount(dem), sizeof *fd.target);
  fd.fraction = calloc(GridElementCount(dem), sizeof *fd.fraction);

  return fd;
}

Flow FlowRouteD8Carve(Grid dem) {
  Flow fd = FlowInit(dem);

  Grid demf = GridFillsinks(dem);
  Grid flats = GridIdentifyFlats(demf);
  Grid costs = GridGWDTComputeCosts(flats, dem, demf);
  Grid dist = GridGWDT(costs, flats);

  if (demf.data && flats.data && costs.data && dist.data) {
    Grid direction = GridCreate(GridU8, NULL, dem.cellsize, dem.dims);
    if (direction.data && fd.stream) {
      flow_routing_d8_carve(fd.stream, direction.data, demf.data, dist.data,
                            flats.data, dem.dims, 1);

      if (fd.source && fd.target && fd.fraction) {
        fd.count = flow_routing_d8_edgelist(fd.source, fd.target, fd.stream,
                                            direction.data, dem.dims, 1);
        for (ptrdiff_t e = 0; e < fd.count; e++) {
          fd.fraction[e] = 1.0;
        }
      }
    }
    GridFree(&direction);
  }

  GridFree(&demf);
  GridFree(&flats);
  GridFree(&costs);
  GridFree(&dist);
  return fd;
}

Grid FlowAccumulation(Flow f) {
  Grid a = GridFill(1.0, f.cellsize, f.dims);

  if (a.data && f.count > 0) {
    traverse_down_f32_add_mul(a.data, f.fraction, f.source, f.target, f.count);
  }

  return a;
}

Grid Outlets(Flow f) {
  Grid outlets = GridCreate(GridU8, NULL, f.cellsize, f.dims);
  if (outlets.data && f.count > 0) {
    uint8_t *outletdata = outlets.data;

    Grid outdegree = GridCreate(GridU8, NULL, f.cellsize, f.dims);
    if (outdegree.data) {
      uint8_t *outdegreedata = outdegree.data;

      ptrdiff_t dims[2] = {f.dims[0], f.dims[1]};

      edgelist_degree(outletdata, outdegreedata, f.source, f.target,
                      GridElementCount(outlets), f.count);

      for (ptrdiff_t j = 0; j < outlets.dims[1]; j++) {
        for (ptrdiff_t i = 0; i < outlets.dims[0]; i++) {
          outletdata[j * dims[0] + i] = (outletdata[j * dims[0] + i] > 0) &
                                        (outdegreedata[j * dims[0] + i] == 0);
        }
      }

      GridFree(&outdegree);
    }
  }
  return outlets;
}

Grid FlowDrainageBasins(Flow f) {
  Grid basins = GridCreate(GridU32, NULL, f.cellsize, f.dims);

  if (f.count > 0 && basins.data) {
    uint32_t *basinsdata = basins.data;

    uint32_t *weights = calloc(f.count, sizeof *weights);

    if (weights) {
      for (ptrdiff_t e = 0; e < f.count; e++) {
        weights[e] = 0xffffffff;
      }

      Grid outlets = Outlets(f);
      if (outlets.data) {
        uint8_t *outdata = outlets.data;

        uint32_t basincount = 1;
        for (ptrdiff_t j = 0; j < basins.dims[1]; j++) {
          for (ptrdiff_t i = 0; i < basins.dims[0]; i++) {
            if (outdata[j * basins.dims[0] + i]) {
              basinsdata[j * basins.dims[0] + i] = basincount;
              basincount++;
            }
          }
        }

        traverse_up_u32_or_and(basinsdata, weights, f.source, f.target,
                               f.count);

        GridFree(&outlets);
      }
      free(weights);
    }
  }
  return basins;
}
