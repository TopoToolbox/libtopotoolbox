#define TOPOTOOLBOX_BUILD

#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>

#include "gf_utils.h"
#include "graphflood/define_types.h"
#include "topotoolbox.h"

// Helper functions ton mimic c++ std::min and std::max
inline float max_float(float a, float b) { return (a > b) ? a : b; }
inline float min_float(float a, float b) { return (a < b) ? a : b; }

/*
Internal function running graphflood in its full vanilla version in single flow
direction
*/
void _graphflood_full_sfd(GF_FLOAT* Z, GF_FLOAT* hw, uint8_t* BCs,
                          GF_FLOAT* Precipitations, GF_FLOAT* manning,
                          GF_UINT* dim, GF_FLOAT dt, GF_FLOAT dx, bool SFD,
                          bool D8, GF_UINT N_iterations, GF_FLOAT step) {
  // Creating an array of Zw (hydraulic surface = Z + hw)
  GF_FLOAT* Zw = (GF_FLOAT*)malloc(sizeof(GF_FLOAT) * nxy(dim));
  for (GF_UINT i = 0; i < nxy(dim); ++i) Zw[i] = Z[i] + hw[i];

  // Init the graph structure locally
  GF_UINT* Sreceivers = (GF_UINT*)malloc(sizeof(GF_UINT) * nxy(dim));
  GF_FLOAT* distToReceivers = (GF_FLOAT*)malloc(sizeof(GF_FLOAT) * nxy(dim));
  GF_UINT* Sdonors =
      (GF_UINT*)malloc(sizeof(GF_UINT) * nxy(dim) * (D8 ? 8 : 4));
  uint8_t* NSdonors = (uint8_t*)malloc(sizeof(uint8_t) * nxy(dim));
  GF_UINT* Stack = (GF_UINT*)malloc(sizeof(GF_UINT) * nxy(dim));
  GF_FLOAT* Qwin = (GF_FLOAT*)malloc(sizeof(GF_FLOAT) * nxy(dim));

  GF_FLOAT cell_area = dx * dx;

  for (GF_UINT iteration = 0; iteration < N_iterations; ++iteration) {
    // At each iteration I update the graph while filling every depressions (*in
    // the hydraulic surface) with water
    compute_sfgraph_priority_flood(Zw, Sreceivers, distToReceivers, Sdonors,
                                   NSdonors, Stack, BCs, dim, dx, D8, step);

    // From the graph hence created I accumulate the flow (steady conditions)
    compute_weighted_drainage_area_single_flow(Qwin, Precipitations, Sreceivers,
                                               Stack, dim, dx);

    for (GF_UINT i = 0; i < nxy(dim); ++i) {
      // Traversing the stack in reverse, super important because it allows us
      // to update the Zw on the go (it ensures receivers are never processed
      // before the donors and therefor the hydraulic slope remains explicit
      // even if we update a donor)
      GF_UINT node = Stack[nxy(dim) - i - 1];
      GF_UINT rec = Sreceivers[node];

      // Checking if the node needs to be processed
      // Note that a lot of the checks are actually already done by the graph
      // calculation
      if (rec == node) continue;
      // Boundary condition: if the flow can out I do not touch hw
      if (can_out(node, BCs)) continue;
      // Additional check: if no water and no input, no need to calculate
      if (Zw[node] == Z[node] && Qwin[node] == 0) continue;

      // Calculating the hydraulic slope
      GF_FLOAT tSw =
          min_float(Zw[node] - Zw[rec], (GF_FLOAT)1e-6) / distToReceivers[node];

      // Calculating the Volumetric discharge based on Manning's friction
      // equation
      GF_FLOAT tQwout =
          (GF_FLOAT)(distToReceivers[node] / manning[node] *
                     pow(Zw[node] - Z[node], 5. / 3.) * sqrt(tSw));

      // Applying the divergence
      Zw[node] =
          max_float(Z[node], Zw[node] + dt * (Qwin[node] - tQwout) / cell_area);
    }
  }

  // back translate Zw into hw
  for (GF_UINT i = 0; i < nxy(dim); ++i) hw[i] = -Z[i] + Zw[i];

  // Don't forget to free memory
  free(Zw);
  free(Qwin);
  free(Sreceivers);
  free(distToReceivers);
  free(Sdonors);
  free(NSdonors);
  free(Stack);
}

void _graphflood_full_mfd(GF_FLOAT* Z, GF_FLOAT* hw, uint8_t* BCs,
                          GF_FLOAT* Precipitations, GF_FLOAT* manning,
                          GF_UINT* dim, GF_FLOAT dt, GF_FLOAT dx, bool SFD,
                          bool D8, GF_UINT N_iterations, GF_FLOAT step) {
  // Initialising the offset for neighbouring operations
  GF_INT offset[8];
  (D8 == false) ? generate_offset_D4_flat(offset, dim)
                : generate_offset_D8_flat(offset, dim);
  // // Initialising the offset distance for each neighbour
  GF_FLOAT offdx[8];
  (D8 == false) ? generate_offsetdx_D4(offdx, dx)
                : generate_offsetdx_D8(offdx, dx);

  GF_FLOAT dxy = (GF_FLOAT)sqrt(2) * dx;
  GF_FLOAT cell_area = dx * dx;

  // Creating an array of Zw (hydraulic surface = Z + hw)
  GF_FLOAT* Zw = (GF_FLOAT*)malloc(sizeof(GF_FLOAT) * nxy(dim));
  for (GF_UINT i = 0; i < nxy(dim); ++i) Zw[i] = Z[i] + hw[i];

  GF_FLOAT* Qwin = (GF_FLOAT*)malloc(sizeof(GF_FLOAT) * nxy(dim));

  GF_UINT* Stack = (GF_UINT*)malloc(sizeof(GF_UINT) * nxy(dim));

  // reintialising Qw
  for (GF_UINT i = 0; i < nxy(dim); ++i) {
    Qwin[i] = 0.;
    Stack[i] = i;
  }

  GF_FLOAT weights[8];

  for (GF_UINT iteration = 0; iteration < N_iterations; ++iteration) {
    // First priority flooding and calculating stack
    compute_priority_flood_plus_topological_ordering(Zw, Stack, BCs, dim, D8,
                                                     step);

    // reintialising Qw
    for (GF_UINT i = 0; i < nxy(dim); ++i) {
      Qwin[i] = 0.;
    }

    // printf("%u\n", Stack[560250]);

    // processing nodes from top to bottom
    for (GF_UINT i = 0; i < nxy(dim); ++i) {
      // Traversing the stack in reverse, super important because it allows us
      // to update the Zw on the go (it ensures receivers are never processed
      // before the donors and therefor the hydraulic slope remains explicit
      // even if we update a donor)
      GF_UINT node = Stack[nxy(dim) - i - 1];

      // If no data: pass
      if (is_nodata(node, BCs)) continue;
      // Boundary condition: if the flow can out I do not touch hw
      if (can_out(node, BCs)) continue;

      // First, incrementing local Qwin
      Qwin[node] += Precipitations[node] * dx * dx;

      // Now calculating the gradients: local, steepest and weighted
      GF_FLOAT sumslope = 0., maxslope = 0., dxmaxdir = dx;
      for (uint8_t n = 0; n < N_neighbour(D8); ++n) {
        // Checking if the neighbour belongs to the grid
        if (check_bound_neighbour(node, n, dim, BCs, D8) == false) {
          weights[n] = 0;
          continue;
        }

        GF_UINT nnode = node + offset[n];

        if (Zw[nnode] >= Zw[node] || can_receive(nnode, BCs) == false ||
            can_give(node, BCs) == false) {
          weights[n] = 0;
          continue;
        }

        GF_FLOAT tSw =
            max_float((GF_FLOAT)1e-4, (Zw[node] - Zw[nnode]) / offdx[n]);

        weights[n] = tSw * ((dx == offdx[n] || D8 == false) ? dx : dxy);

        sumslope += weights[n];
        if (tSw > maxslope) {
          maxslope = tSw;
          dxmaxdir = offdx[n];
        }
      }

      // Transferring the flux
      // GF_FLOAT sumtransfer = 0.; // DEBUG to keep
      if (sumslope > 0) {
        for (GF_UINT n = 0; n < N_neighbour(D8); ++n) {
          if (weights[n] == 0) continue;
          Qwin[node + offset[n]] += weights[n] / sumslope * Qwin[node];
          // sumtransfer += weights[n] / sumslope;
        }
      }

      // DEBUG STATEMENT: to keep so far even if commented. Will remove when ok.
      // if(fabs(sumtransfer - 1.) > 0.01 && can_out(node,BCs) == false){
      // 	printf("HAPPENS - %f\n",abs(sumtransfer - 1.));
      // }

      // Calculating the Volumetric discharge based on Manning's friction
      // equation
      GF_FLOAT tQwout =
          (GF_FLOAT)(dxmaxdir / manning[node] *
                     pow(Zw[node] - Z[node], 5. / 3.) * sqrt(maxslope));

      // if(Qwin[node] > 0){
      // 	printf("%f\n", Qwin[node]);
      // }

      // Applying the divergence
      // printf("%f", dt*(Qwin[node] - tQwout)/cell_area);
      Zw[node] =
          max_float(Z[node], Zw[node] + dt * (Qwin[node] - tQwout) / cell_area);
      // printf("vs %f", tQwout);
    }
  }

  // back translate Zw into hw
  for (GF_UINT i = 0; i < nxy(dim); ++i) hw[i] = (GF_FLOAT)(-Z[i] + Zw[i]);

  free(Zw);
  free(Qwin);
  free(Stack);
}

/*
        See topotoolbox.h for full instructions
*/
TOPOTOOLBOX_API
void graphflood_full(GF_FLOAT* Z, GF_FLOAT* hw, uint8_t* BCs,
                     GF_FLOAT* Precipitations, GF_FLOAT* manning, GF_UINT* dim,
                     GF_FLOAT dt, GF_FLOAT dx, bool SFD, bool D8,
                     GF_UINT N_iterations, GF_FLOAT step) {
  // Runs the single flow version of the algorithm
  if (SFD)
    _graphflood_full_sfd(Z, hw, BCs, Precipitations, manning, dim, dt, dx, SFD,
                         D8, N_iterations, step);
  else
    _graphflood_full_mfd(Z, hw, BCs, Precipitations, manning, dim, dt, dx, SFD,
                         D8, N_iterations, step);
}
