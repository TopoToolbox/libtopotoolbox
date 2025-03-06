/*
This file contains routine to build a single flow graph.
A Single flow graph is a data structure describing a Directed Acyclic Graph
(DAG) where nodes only have one receiver max, usually the steepest gradient one.

*/

#define TOPOTOOLBOX_BUILD

#include <float.h>
#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>

#include "gf_utils.h"
#include "graphflood/define_types.h"
#include "pq_priority_flood.h"
#include "queue_pit.h"
#include "topotoolbox.h"

static void recursive_stack(GF_UINT node, GF_UINT* Sdonors, GF_UINT* Stack,
                            uint8_t* NSdonors, GF_UINT* istack, bool D8);

/*
        Computes a single flow graph with minimal characteristics:
        - List of single flow receivers
        - Number of single flow donors
        - the topologically ordered stack (sensu Braun and Willett, 2013)
*/
TOPOTOOLBOX_API
void compute_sfgraph(float* topo, GF_UINT* Sreceivers,
                     GF_FLOAT* distToReceivers, GF_UINT* Sdonors,
                     uint8_t* NSdonors, GF_UINT* Stack, uint8_t* BCs,
                     GF_UINT* dim, float dx, bool D8) {
  // Initialising the offset for neighbouring operations
  GF_INT offset[8];
  (D8 == false) ? generate_offset_D4_flat(offset, dim)
                : generate_offset_D8_flat(offset, dim);
  // // Initialising the offset distance for each neighbour
  float offdx[8];
  (D8 == false) ? generate_offsetdx_D4(offdx, dx)
                : generate_offsetdx_D8(offdx, dx);

  // For all the nodes
  // in row major d0 is row and d1 is col
  // in col major d0 is col and d1 is row
  for (GF_UINT d0 = 0; d0 < dim[0]; ++d0) {
    for (GF_UINT d1 = 0; d1 < dim[1]; ++d1) {
      // Getting flat index of the node
      GF_UINT node = dim2flat(d0, d1, dim);

      // By convention (see fastscape, LSDTT, ...) a no steepest receiver =
      // itself
      Sreceivers[node] = node;
      distToReceivers[node] = 0.;
      NSdonors[node] = 0;

      // Boundary condition checks: the node needs to being able to give to have
      // receivers Note that nodata cannot give so it filter them too
      if (can_give(node, BCs) == false) continue;

      // Targetting the steepest receiver
      // -> Initialising the node to itself (no receivers)
      GF_UINT this_receiver = node;
      GF_FLOAT this_receiverdx = 0.;
      // -> Initialising the slope to 0
      float SD = 0.;

      // for all the neighbours ...
      for (uint8_t n = 0; n < N_neighbour(D8); ++n) {
        // Checking if the neighbour belongs to the grid
        if (check_bound_neighbour(node, n, dim, BCs, D8) == false) {
          continue;
        }
        // flat indices
        GF_UINT nnode = node + offset[n];
        // who can receive
        if (can_receive(nnode, BCs) == false) continue;

        // I check wether their slope is the steepest
        float tS = (topo[node] - topo[nnode]) / offdx[n];

        // if it is
        if (tS > SD) {
          // I save it
          this_receiver = nnode;
          this_receiverdx = offdx[n];
          SD = tS;
        }
      }

      // and the final choice is saved
      Sreceivers[node] = this_receiver;
      distToReceivers[node] = this_receiverdx;
    }
  }

  // Back calculating the number of steepest receivers and inverting the
  // receivers
  for (GF_UINT node = 0; node < dim[0] * dim[1]; ++node) {
    if (node != (GF_UINT)Sreceivers[node]) {
      Sdonors[Sreceivers[node] * N_neighbour(D8) + NSdonors[Sreceivers[node]]] =
          (int32_t)node;
      ++NSdonors[Sreceivers[node]];
    }
  }

  // Finally calculating Braun and Willett 2013
  GF_UINT istack = 0;
  for (GF_UINT node = 0; node < dim[0] * dim[1]; ++node) {
    if (node == (GF_UINT)Sreceivers[node]) {
      recursive_stack(node, Sdonors, Stack, NSdonors, &istack, D8);
    }
  }
}

/*
Recursive function to build Braun and Willett's Stack (single flow topological
ordering) for each donors of a node it successively include them to the stack
and call itself on the donor
*/
static void recursive_stack(GF_UINT node, GF_UINT* Sdonors, GF_UINT* Stack,
                            uint8_t* NSdonors, GF_UINT* istack, bool D8) {
  Stack[*istack] = node;
  ++(*istack);
  for (GF_UINT nd = 0; nd < NSdonors[node]; ++nd) {
    recursive_stack(Sdonors[node * N_neighbour(D8) + nd], Sdonors, Stack,
                    NSdonors, istack, D8);
  }
}

/*
See topotoolbox.h for more details.
Computes both the single flow graph data structure, the topological ordering and
fills the local minimas with Priority Flood + epsilon using Barnes (2014)
receivers sfgraph structure built as the surface field fills up with PF, then
Note that it still runs the Stack ordering from Brun and Willett 2013 at the end
to ensure a stack segmented by watersheds (the PitQueue data structure speeds up
the filling but breaks the node processing in ascending order within
depressions)
*/
TOPOTOOLBOX_API
void compute_sfgraph_priority_flood(float* topo, GF_UINT* Sreceivers,
                                    GF_FLOAT* distToReceivers, GF_UINT* Sdonors,
                                    uint8_t* NSdonors, GF_UINT* Stack,
                                    uint8_t* BCs, GF_UINT* dim, float dx,
                                    bool D8, GF_FLOAT step) {
  // Initialising the offset for neighbouring operations
  GF_INT offset[8];
  (D8 == false) ? generate_offset_D4_flat(offset, dim)
                : generate_offset_D8_flat(offset, dim);
  // // Initialising the offset distance for each neighbour
  float offdx[8];
  (D8 == false) ? generate_offsetdx_D4(offdx, dx)
                : generate_offsetdx_D8(offdx, dx);

  // Reinitialising all the Ndonors to 0
  // initialising the nodes to not closed ( = to be processed)
  uint8_t* closed = (uint8_t*)malloc(nxy(dim) * sizeof(uint8_t));
  for (GF_UINT i = 0; i < nxy(dim); ++i) {
    NSdonors[i] = 0;
    closed[i] = false;
  }

  // PitQueue is a FIFO data structure to fill pits without having to use the
  // more expensive priority queue
  PitQueue pit;
  pitqueue_init(&pit, (int)nxy(dim));

  // The priority queue data structure (keeps stuff sorted)
  PFPQueue open;
  pfpq_init(&open, nxy(dim));

  // temp variable to help with PitQueue
  float PitTop = FLT_MIN;

  // Initialisation phase: initialise the queue with nodethat can drain out of
  // the model Also initialise the sfg data structure
  for (GF_UINT i = 0; i < nxy(dim); ++i) {
    // By convention (see fastscape, LSDTT, ...) a no steepest receiver = itself
    Sreceivers[i] = i;
    distToReceivers[i] = 0.;

    // If flow can leave, I push
    if (can_out(i, BCs)) {
      pfpq_push(&open, i, topo[i]);
      closed[i] = true;
    }

    // Note that no data node are immediately closed as processed
    if (is_nodata(i, BCs)) {
      closed[i] = true;
    }
  }

  // Here we go: Starting the main process
  // Processing stops once all the nodes - nodata have been visited once (i.e.
  // pit fifo and PQ empty)
  GF_UINT node;
  while (pfpq_empty(&open) == false || pit.size > 0) {
    // Selecting the next node
    if (pit.size > 0 && pfpq_empty(&open) == false &&
        pfpq_top_priority(&open) == topo[pit.front]) {
      node = pfpq_pop_and_get_key(&open);
      PitTop = FLT_MIN;

    } else if (pit.size > 0) {
      node = pitqueue_pop_and_get(&pit);
      if (PitTop == FLT_MIN) PitTop = topo[node];
    } else {
      node = pfpq_pop_and_get_key(&open);
      PitTop = FLT_MIN;
    }

    // Check is the Sreceiver has not already been imposed to a processed node
    // if need_update is false, then I don;t chose a new receiver as it means I
    // am in a pit
    bool need_update = Sreceivers[node] == node;

    // Targetting the steepest receiver
    // -> Initialising the node to itself (no receivers)
    GF_UINT this_receiver = node;
    GF_FLOAT this_receiverdx = 0.;
    // -> Initialising the slope to 0
    float SD = 0.;

    // for all the neighbours ...
    for (uint8_t n = 0; n < N_neighbour(D8); ++n) {
      // Checking if the neighbour belongs to the grid
      if (check_bound_neighbour(node, n, dim, BCs, D8) == false) {
        continue;
      }

      // flat indices
      GF_UINT nnode = node + offset[n];

      // if nodata I skip
      if (is_nodata(nnode, BCs)) continue;

      // This section process the graph structure (if needed)
      // Note that it can only be a Sreceiver if already closed
      // otherwise it means it is either a donor or in a pit
      if (can_receive(nnode, BCs) && can_give(node, BCs) && need_update &&
          closed[nnode]) {
        // I check wether their slope is the steepest
        float tS = (topo[node] - topo[nnode]) / offdx[n];

        // if it is
        if (tS > SD) {
          // I save it
          this_receiver = nnode;
          this_receiverdx = offdx[n];
          SD = tS;
        }
      }

      // If the node is closed (i.e. already in a pit or processed) I skip
      if (closed[nnode] == false) {
        // other wise I close it
        closed[nnode] = true;

        // I raise its elevation if is in pit
        // nextafter maskes sure I pick the next floating point data
        // corresponding to the current precision
        if (topo[nnode] <=
            (GF_FLOAT)nextafter((GF_FLOAT)topo[node], (GF_FLOAT)FLT_MAX) +
                step) {
          // raise
          topo[nnode] =
              (GF_FLOAT)nextafter((GF_FLOAT)topo[node], (GF_FLOAT)FLT_MAX) +
              step;
          // put in pit queue
          pitqueue_enqueue(&pit, nnode);
          // Affect current node as neighbours Sreceiver
          Sreceivers[nnode] = node;
          distToReceivers[nnode] = offdx[n];
        } else {
          // ... Not in a pit? then in PQ for next proc
          pfpq_push(&open, nnode, topo[nnode]);
        }
      }
    }

    // Updating SFG data structure if needed
    if (need_update) {
      // and the final choice is saved
      Sreceivers[node] = this_receiver;
      distToReceivers[node] = this_receiverdx;
    }
  }

  // Done with the queues and close, free memory
  pfpq_free(&open);
  pitqueue_free(&pit);
  free(closed);

  // Back calculating the number of steepest receivers and inverting the
  // receivers
  for (GF_UINT tnode = 0; tnode < dim[0] * dim[1]; ++tnode) {
    if (tnode != (GF_UINT)Sreceivers[tnode]) {
      Sdonors[Sreceivers[tnode] * N_neighbour(D8) +
              NSdonors[Sreceivers[tnode]]] = tnode;
      ++NSdonors[Sreceivers[tnode]];
    }
  }

  // Finally calculating Braun and Willett 2013
  GF_UINT istack = 0;
  for (GF_UINT tnode = 0; tnode < dim[0] * dim[1]; ++tnode) {
    if (tnode == (GF_UINT)Sreceivers[tnode]) {
      recursive_stack(tnode, Sdonors, Stack, NSdonors, &istack, D8);
    }
  }
}
