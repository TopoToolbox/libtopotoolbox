/*
This file contains routine to build a single flow graph.
A Single flow graph is a data structure describing a Directed Acyclic Graph (DAG)
where nodes only have one receiver max, usually the steepest gradient one.

*/

#define TOPOTOOLBOX_BUILD

#include <math.h>
#include <stddef.h>
#include <stdint.h>

#include "../morphology/reconstruct.h"
#include "gf_utils.h"
#include "topotoolbox.h"



static void recursive_stack(int32_t node, int32_t* Sdonors, int32_t* Stack, uint8_t* NSdonors, size_t istack, bool D8);

/*
  Computes a single flow graph with minimal characteristics:
  - List of single flow receivers
  - Number of single flow donors
  - the topologically ordered stack (sensu Braun and Willett, 2013) 
*/
TOPOTOOLBOX_API
void compute_sfgraph(float* topo, int32_t* Sreceivers, int32_t* Sdonors, uint8_t* NSdonors, uint32_t* Stack, uint8_t* BCs, uint32_t* dim, float dx, bool D8) {
  
  printf("DEBUG::1, entering function");


  // // Initialising the offset for neighbouring operations
  // int32_t* offset = NULL;
  // (D8 == false) ? generate_offset_D4_flat(&offset,dim) : generate_offset_D8_flat(&offset,dim);
  // // Initialising the offset distance for each neighbour
  // float* offdx = NULL;
  // (D8 == false) ? generate_offsetdx_D4(&offdx,dx) : generate_offsetdx_D8(&offdx,dx);

  // printf("DEBUG::OFFSET::%s\n", offset[0]);

  // For all the nodes
  // in row major d0 is row and d1 is col
  // in col major d0 is col and d1 is row
  // for(int32_t d0 = 0; d0<dim[0]; ++d0){
  //   for(int32_t d1 = 0; d1<dim[0]; ++d1){

  //     // Getting flat index of the node
  //     int32_t node = dim2flat(d0,d1,dim);

  //     // By convention (see fastscape, LSDTT, ...) a no steepest receiver = itself
  //     Sreceivers[node] = node;

  //     // Boundary condition checks: the node needs to being able to give to have receivers
  //     // Note that nodata cannot give so it filter them too
  //     if(can_give(node,BCs) == false)
  //       continue;

  //     // Targetting the steepest receiver
  //     // -> Initialising the node to itself (no receivers)
  //     int32_t this_receiver = node;
  //     // -> Initialising the slope to 0
  //     float SD = 0.;

  //     // for all the neighbours ...
  //     for(size_t n = 0; n<N_neighbour(D8); ++n){
  //       // flat indices
  //       int32_t nnode = node + offset[n];
  //       // who can receive
  //       if(can_receive(nnode, BCs) == false)
  //         continue;

  //       // I check wether their slope is the steepest
  //       float tS = (topo[node] - topo[nnode])/offdx[n];

  //       // if it is
  //       if(tS > SD){
  //         // I save it
  //         this_receiver = nnode;
  //         SD = 0.;
  //       }
  //     }

  //     // and the final choice is saved
  //     Sreceivers[node] = this_receiver;
  //   }
  // }
  

  // // Back calculating the number of steepest receivers and inverting the receivers
  // for(size_t node = 0; node<dim[0]*dim[1]; ++node){
  //   if(node != Sreceivers[node]){
  //     Sdonors[Sreceivers[node] * N_neighbour(D8) + NSdonors[Sreceivers[node]]] = node;
  //     ++NSdonors[Sreceivers[node]];
  //   }
  // }

  // // Finally calculating Braun and Willett 2013
  // size_t istack = 0;
  // for(size_t node = 0; node<dim[0]*dim[1]; ++node){
  //   if(node == Sreceivers[node]){
  //     recursive_stack(node, Sdonors, Stack, NSdonors, istack, D8);
  //   }

  // }

}



/*
Recursive function to build Braun and Willett's Stack (single flow topological ordering)
for each donors of a node it successively include them to the stack and call itself on the donor
*/
static void recursive_stack(int32_t node, int32_t* Sdonors, int32_t* Stack, uint8_t* NSdonors, size_t istack, bool D8){
  Stack[istack] = node;
  ++istack;
  for(uint32_t nd=0; nd<NSdonors[node]; ++nd){
    recursive_stack(Sdonors[node * N_neighbour(D8) + nd], Sdonors, Stack, NSdonors, istack,D8);
  }
}
