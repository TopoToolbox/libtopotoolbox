/*
This file contains routine to accumulate flow downstream, a way or another

*/

#define TOPOTOOLBOX_BUILD

#include <math.h>
#include <stddef.h>
#include <stdint.h>

#include "gf_utils.h"
#include "graphflood/define_types.h"
#include "topotoolbox.h"

/*
Calculate the drainage area from a topological order
*/
TOPOTOOLBOX_API
void compute_drainage_area_single_flow(GF_FLOAT* output, GF_UINT* Sreceivers,
                                       GF_UINT* Stack, GF_UINT* dim,
                                       GF_FLOAT dx) {
  GF_UINT tnxy = nxy(dim);
  const GF_FLOAT cell_area = dx * dx;

  for (GF_UINT i = 0; i < tnxy; ++i) output[i] = 0.;

  for (GF_UINT i = 0; i < tnxy; ++i) {
    GF_UINT ri = tnxy - 1 - i;
    GF_UINT node = Stack[ri];
    if (node == Sreceivers[node]) continue;
    output[node] += cell_area;
    output[Sreceivers[node]] += output[node];
  }
}

/*
Calculate the drainage area from a topological order
*/
TOPOTOOLBOX_API
void compute_weighted_drainage_area_single_flow(GF_FLOAT* output,
                                                GF_FLOAT* weights,
                                                GF_UINT* Sreceivers,
                                                GF_UINT* Stack, GF_UINT* dim,
                                                GF_FLOAT dx) {
  GF_UINT tnxy = nxy(dim);
  const GF_FLOAT cell_area = dx * dx;

  for (GF_UINT i = 0; i < tnxy; ++i) output[i] = 0.;

  for (GF_UINT i = 0; i < tnxy; ++i) {
    GF_UINT ri = tnxy - 1 - i;
    GF_UINT node = Stack[ri];
    if (node == Sreceivers[node]) continue;
    output[node] += cell_area * weights[node];
    output[Sreceivers[node]] += output[node];
  }
}
