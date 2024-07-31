#define TOPOTOOLBOX_BUILD

#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>

#include "graphflood/define_types.h"
#include "gf_utils.h"

#include "../morphology/reconstruct.h"
#include "topotoolbox.h"

/*
  todo
*/
TOPOTOOLBOX_API
void graphflood_full(GF_FLOAT* Z, GF_FLOAT* hw, uint8_t* BCs, GF_FLOAT* Precipitations, GF_FLOAT* manning, GF_UINT* dim, GF_FLOAT dt, GF_FLOAT dx, bool SFD, bool D8, int N_iterations) {

  if(SFD == false){
    printf("Multiple Flow Direction is work in progress and won't work yet\n");
    return;
  }

  // Creating an array of Zw (hydraulic surface = Z + hw)
  GF_FLOAT* Zw = (GF_FLOAT*)malloc(sizeof(GF_FLOAT) * nxy(dim));
  for(GF_UINT i=0; i < nxy(dim); ++i)
    Zw[i] = Z[i] + hw[i];


  // Init the graph structure locally
  GF_UINT* Sreceivers = (GF_UINT*)malloc(sizeof(GF_UINT) * nxy(dim));
  GF_FLOAT* distToReceivers = (GF_FLOAT*)malloc(sizeof(GF_FLOAT) * nxy(dim));
  GF_UINT* Sdonors = (GF_UINT*)malloc(sizeof(GF_UINT) * nxy(dim) * (D8 ? 8:4) );
  uint8_t* NSdonors = (uint8_t*)malloc(sizeof(uint8_t) * nxy(dim));
  GF_UINT* Stack = (GF_UINT*)malloc(sizeof(GF_UINT) * nxy(dim));
  GF_FLOAT* Qwin = (GF_FLOAT*)malloc(sizeof(GF_FLOAT) * nxy(dim));

  GF_FLOAT cell_area = dx*dx;


  for(GF_UINT iteration = 0; iteration<N_neighbour; ++iteration){

    // At each iteration I update the graph while filling every depressions (*in the hydraulic surface) with water
    compute_sfgraph_priority_flood(Zw, Sreceivers, distToReceivers, Sdonors, NSdonors, Stack, BCs, dim, dx, D8);

    // From the graph hence created I accumulate the flow (steady conditions)
    compute_weighted_drainage_area_single_flow(Qwin, Precipitations, Sreceivers, Stack, dim, dx);

    for(GF_UINT i=0; i<nxy(dim);++i){

      // Traversing the stack in reverse, super important because it allows us to update the Zw on the go
      // (it ensures receivers are never processed before the donors and therefor the hydraulic slope remains explicit even if we update a donor)
      GF_UINT node = Stack[nxy(dim) - i - 1];
      GF_UINT rec = Sreceivers[node];

      // Checking if the node needs to be processed 
      // Note that a lot of the checks are actually already done by the graph calculation
      if(rec == node) continue;
      // Additional check: if no water and no input, no need to calculate
      if(Zw[node] == Z[node] && Qwin[node] == 0) continue;

      // Calculating the hydraulic slope
      GF_FLOAT tSw = min(Zw[node] - Zw[rec], (GF_FLOAT)1e-6)/distToReceivers[node];
      
      // Calculating the Volumetric discharge based on Manning's friction equation
      GF_FLOAT tQwout = distToReceivers[node]/manning[node] * pow(Zw[node] - Z[node], 5./3.) * sqrt(tSw);

      // Applying the divergence
      Zw[node] = max(Z[node], Zw[node] + dt*(Qwin[node] - tQwout)/cell_area);
    }
  }


  // back translate Zw into hw
  for(GF_UINT i=0; i < nxy(dim); ++i)
    hw[i] = -Z[i] + Zw[i];

  // don't forget to free memory
  free(Zw);
  free(Qwin);
  free(Sreceivers);
  free(distToReceivers);
  free(Sdonors);
  free(NSdonors);
  free(Stack);
  
}
