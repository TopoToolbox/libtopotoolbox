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

#include "../morphology/reconstruct.h"
#include "gf_utils.h"
#include "graphflood/define_types.h"
#include "pq_priority_flood.h"
#include "queue_pit.h"
#include "topotoolbox.h"

/*
See topotoolbox.h for details
LSS: runs priority floods (Barnes 2014) + epsilon on a topogrpahy
Fills it in place, and take into D8/D4 topology and custom boundary conditions
*/
TOPOTOOLBOX_API
void compute_priority_flood(float* topo, uint8_t* BCs, GF_UINT* dim, bool D8) {
  // Initialising the offset for neighbouring operations
  GF_INT offset[8];
  (D8 == false) ? generate_offset_D4_flat(offset, dim)
                : generate_offset_D8_flat(offset, dim);

  // initialising the nodes to not closed ( = to be processed)
  uint8_t* closed = (uint8_t*)malloc(nxy(dim) * sizeof(uint8_t));
  for (GF_UINT i = 0; i < nxy(dim); ++i) closed[i] = false;

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

      // If the node is closed (i.e. already in a pit or processed) I skip
      if (closed[nnode] == false) {
        // other wise I close it
        closed[nnode] = true;

        // I raise its elevation if is in pit
        // nextafter maskes sure I pick the next floating point data
        // corresponding to the current precision
        if (topo[nnode] <=
            (GF_FLOAT)nextafter((GF_FLOAT)topo[node], (GF_FLOAT)FLT_MAX)) {
          // raise
          topo[nnode] =
              (GF_FLOAT)nextafter((GF_FLOAT)topo[node], (GF_FLOAT)FLT_MAX);
          // put in pit queue
          pitqueue_enqueue(&pit, nnode);
          // Affect current node as neighbours Sreceiver
        } else {
          // ... Not in a pit? then in PQ for next proc
          pfpq_push(&open, nnode, topo[nnode]);
        }
      }
    }
  }

  // Done with the queues and close, free memory
  pfpq_free(&open);
  pitqueue_free(&pit);
  free(closed);
}

/*
See topotoolbox.h for details
Variant of priority flood + epsilon that also computes a multiple flow
compatible topological ordering (Stack) on the go Note: this version is a bit
slower as it needs to push all nodes in the priority queue including the ones in
pit (no fifo queue) By experience the slowing down factor is small for most
cases ( low amount of depressions, ~10% max slow down)
*/
TOPOTOOLBOX_API
void compute_priority_flood_plus_topological_ordering(float* topo,
                                                      GF_UINT* stack,
                                                      uint8_t* BCs,
                                                      GF_UINT* dim, bool D8) {
  // Initialising the offset for neighbouring operations
  GF_INT offset[8];
  (D8 == false) ? generate_offset_D4_flat(offset, dim)
                : generate_offset_D8_flat(offset, dim);

  // initialising the nodes to not closed ( = to be processed)
  uint8_t* closed = (uint8_t*)malloc(nxy(dim) * sizeof(uint8_t));
  for (GF_UINT i = 0; i < nxy(dim); ++i) closed[i] = false;

  // The priority queue data structure (keeps stuff sorted)
  PFPQueue open;
  pfpq_init(&open, nxy(dim));

  GF_UINT istack = 0;

  // Initialisation phase: initialise the queue with nodethat can drain out of
  // the model Also initialise the sfg data structure
  for (GF_UINT i = 0; i < nxy(dim); ++i) {
    // If flow can leave, I push
    if (can_out(i, BCs)) {
      pfpq_push(&open, i, topo[i]);
      closed[i] = true;
    }

    // Note that no data node are immediately closed as processed
    if (is_nodata(i, BCs)) {
      closed[i] = true;
      stack[istack] = i;
      ++istack;
    }
  }

  // Here we go: Starting the main process
  // Processing stops once all the nodes - nodata have been visited once (i.e.
  // pit fifo and PQ empty)
  GF_UINT node;
  while (pfpq_empty(&open) == false) {
    // printf("DEBUG::A3\n");
    node = pfpq_pop_and_get_key(&open);

    // printf("%u vs %u\n", istack, nxy(dim));
    if (istack < nxy(dim)) {
      stack[istack] = node;
    }

    ++istack;

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

      // If the node is closed (i.e. already in a pit or processed) I skip
      if (closed[nnode] == false) {
        // other wise I close it
        closed[nnode] = true;

        // I raise its elevation if is in pit
        // nextafter maskes sure I pick the next floating point data
        // corresponding to the current precision
        if (topo[nnode] <=
            (GF_FLOAT)(nextafter((GF_FLOAT)topo[node], (GF_FLOAT)FLT_MAX) +
                       1e-4)) {
          // raise
          topo[nnode] =
              (GF_FLOAT)(nextafter((GF_FLOAT)topo[node], (GF_FLOAT)FLT_MAX) +
                         1e-4);
          // put in pqueue
          pfpq_push(&open, nnode, topo[nnode]);
          // Affect current node as neighbours Sreceiver
        } else {
          // ... Not in a pit? then wimply in PQ for next proc
          pfpq_push(&open, nnode, topo[nnode]);
        }
      }
    }
  }

  // Done with the queues and close, free memory
  pfpq_free(&open);
  free(closed);
}
