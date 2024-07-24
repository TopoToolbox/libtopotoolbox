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
#include "topotoolbox.h"
#include "gf_utils.h"





/*
  Computes a single flow graph with minimal characteristics:
  - List of single flow receivers
  - Number of single flow donors
  - the topologically ordered stack (sensu Braun and Willett, 2013) 
*/
TOPOTOOLBOX_API
void compute_sfgraph(size_t* Sreceivers, uint8_t* NSdonors, size_t* Stack, uint8_t* BCs, size_t* dim, bool D8) {
  printf("YOLO\n");
  
}
