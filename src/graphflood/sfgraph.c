/*
This file contains routine to build a single flow graph.
A Single flow graph is a data structure describing a Directed Acyclic Graph (DAG)
where nodes only have one receiver max, usually the steepest gradient one.

*/

#define TOPOTOOLBOX_BUILD

#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <float.h>

#include "../morphology/reconstruct.h"
#include "gf_utils.h"
#include "pq_priority_flood.h"
#include "queue_pit.h"
#include "topotoolbox.h"



static void recursive_stack(int32_t node, int32_t* Sdonors, uint32_t* Stack, uint8_t* NSdonors, uint32_t* istack, bool D8);

/*
	Computes a single flow graph with minimal characteristics:
	- List of single flow receivers
	- Number of single flow donors
	- the topologically ordered stack (sensu Braun and Willett, 2013) 
*/
TOPOTOOLBOX_API
void compute_sfgraph(float* topo, int32_t* Sreceivers, int32_t* Sdonors, uint8_t* NSdonors, uint32_t* Stack, uint8_t* BCs, uint32_t* dim, float dx, bool D8) {
	// Initialising the offset for neighbouring operations
	int32_t offset[8];
	(D8 == false) ? generate_offset_D4_flat(offset,dim) : generate_offset_D8_flat(offset, dim);
	// // Initialising the offset distance for each neighbour
	float offdx[8];
	(D8 == false) ? generate_offsetdx_D4(offdx,dx) : generate_offsetdx_D8(offdx,dx);


	// For all the nodes
	// in row major d0 is row and d1 is col
	// in col major d0 is col and d1 is row
	for(uint32_t d0 = 0; d0<dim[0]; ++d0){
		for(uint32_t d1 = 0; d1<dim[1]; ++d1){

			// Getting flat index of the node
			int32_t node = dim2flat(d0,d1,dim);

			// By convention (see fastscape, LSDTT, ...) a no steepest receiver = itself
			Sreceivers[node] = node;

			// Boundary condition checks: the node needs to being able to give to have receivers
			// Note that nodata cannot give so it filter them too
			if(can_give(node,BCs) == false)
				continue;

			// Targetting the steepest receiver
			// -> Initialising the node to itself (no receivers)
			int32_t this_receiver = node;
			// -> Initialising the slope to 0
			float SD = 0.;

			// for all the neighbours ...
			for(size_t n = 0; n<N_neighbour(D8); ++n){
				// Checking if the neighbour belongs to the grid
				if(check_bound_neighbour(node, n, dim, BCs, D8) == false){
					continue;
				}
				// flat indices
				int32_t nnode = node + offset[n];
				// who can receive
				if(can_receive(nnode, BCs) == false)
					continue;

				// I check wether their slope is the steepest
				float tS = (topo[node] - topo[nnode])/offdx[n];

				// if it is
				if(tS > SD){
					// I save it
					this_receiver = nnode;
					SD = tS;
				}
			}

			// and the final choice is saved
			Sreceivers[node] = this_receiver;
		}
	}

	// Back calculating the number of steepest receivers and inverting the receivers
	for(uint32_t node = 0; node<dim[0]*dim[1]; ++node){
		if(node != (uint32_t)Sreceivers[node]){
			Sdonors[Sreceivers[node] * N_neighbour(D8) + NSdonors[Sreceivers[node]]] = (int32_t)node;
			++NSdonors[Sreceivers[node]];
		}
	}

	// Finally calculating Braun and Willett 2013
	uint32_t istack = 0;
	for(uint32_t node = 0; node<dim[0]*dim[1]; ++node){
		if(node == (uint32_t)Sreceivers[node]){
			recursive_stack(node, Sdonors, Stack, NSdonors, &istack, D8);
		}

	}

}



/*
Recursive function to build Braun and Willett's Stack (single flow topological ordering)
for each donors of a node it successively include them to the stack and call itself on the donor
*/
static void recursive_stack(int32_t node, int32_t* Sdonors, uint32_t* Stack, uint8_t* NSdonors, uint32_t* istack, bool D8){
	Stack[*istack] = node;
	++(*istack);
	for(uint32_t nd=0; nd<NSdonors[node]; ++nd){
		recursive_stack(Sdonors[node * N_neighbour(D8) + nd], Sdonors, Stack, NSdonors, istack, D8);
	}
}




TOPOTOOLBOX_API
void compute_sfgraph_priority_flood(float* topo, int32_t* Sreceivers, int32_t* Sdonors, uint8_t* NSdonors, uint32_t* Stack, uint8_t* BCs, uint32_t* dim, float dx, bool D8) {
	
	// Initialising the offset for neighbouring operations
	int32_t offset[8];
	(D8 == false) ? generate_offset_D4_flat(offset,dim) : generate_offset_D8_flat(offset, dim);
	// // Initialising the offset distance for each neighbour
	float offdx[8];
	(D8 == false) ? generate_offsetdx_D4(offdx,dx) : generate_offsetdx_D8(offdx,dx);

	uint8_t* closed = (uint8_t*) malloc( nxy(dim) * sizeof(uint8_t) );
	for(uint32_t i=0; i<nxy(dim); ++i)
		closed[i]=false;

	PitQueue pit;
	pitqueue_init(&pit,nxy(dim));

	PFPQueue open;
	pfpq_init(&open, nxy(dim));

	float PitTop = FLT_MIN;

	uint32_t istack = 0;

	for(uint32_t i=0; i<nxy(dim); ++i){
		if(can_out(i,BCs))
			pfpq_push(&open, i, topo[i]);
		if(is_nodata(i,BCs)){
			closed[i] = true;
			Stack[istack] = i;
			++istack;
		}

	}

	uint32_t node;

	while(pfpq_empty(&open) == false || pit.size>0){

		if(pit.size>0 && pfpq_empty(&open)==false && pfpq_top_priority(&open) == topo[pit.front]){

			node=pfpq_pop_and_get_key(&open);
			PitTop=FLT_MIN;
		
		} else if(pit.size>0){

			node=pitqueue_pop_and_get(&pit);
			if(PitTop==FLT_MIN)
				PitTop=topo[node];
		} else {
			node=pfpq_pop_and_get_key(&open);
			PitTop=FLT_MIN;
		}

		Stack[istack] = node;
		++istack;

		// By convention (see fastscape, LSDTT, ...) a no steepest receiver = itself
		Sreceivers[node] = node;
		// Targetting the steepest receiver
		// -> Initialising the node to itself (no receivers)
		int32_t this_receiver = node;
		// -> Initialising the slope to 0
		float SD = 0.;

		// for all the neighbours ...
		for(size_t n = 0; n<N_neighbour(D8); ++n){

			
			// flat indices
			int32_t nnode = node + (int32_t)offset[n];

			// Checking if the neighbour belongs to the grid
			if(check_bound_neighbour(nnode, n, dim, BCs, D8) == false){
				continue;
			}

			if(is_nodata(nnode,BCs)) continue;
			

			// who can receive 
			if(can_receive(nnode, BCs) && can_give(node,BCs) ){
				

				// I check wether their slope is the steepest
				float tS = (topo[node] - topo[nnode])/offdx[n];

				// if it is
				if(tS > SD){
					// I save it
					this_receiver = nnode;
					SD = tS;
				}
			}

			

			if(closed[node] == false){

				closed[node] = true;

				if(topo[nnode]<=nextafter(topo[node],FLT_MAX)){
					topo[nnode]=nextafter(topo[node],FLT_MAX);
					pitqueue_enqueue(&pit,nnode);
				} else
					pfpq_push(&open,nnode,topo[nnode]);
			}

		}
		// and the final choice is saved
		Sreceivers[node] = this_receiver;
	}

	pfpq_free(&open);
	pitqueue_free(&pit);
	free(closed);

	// Back calculating the number of steepest receivers and inverting the receivers
	for(uint32_t node = 0; node<dim[0]*dim[1]; ++node){
		if(node != (uint32_t)Sreceivers[node]){
			Sdonors[Sreceivers[node] * N_neighbour(D8) + NSdonors[Sreceivers[node]]] = node;
			++NSdonors[Sreceivers[node]];
		}
	}




}


