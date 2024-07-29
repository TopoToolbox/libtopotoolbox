/*
This file contains routine to accumulate flow downstream, a way or another

*/

#define TOPOTOOLBOX_BUILD

#include <math.h>
#include <stddef.h>
#include <stdint.h>

#include "../morphology/reconstruct.h"
#include "gf_utils.h"
#include "topotoolbox.h"

/*
Calculate the drainage area from a topological order
*/
TOPOTOOLBOX_API
void compute_drainage_area_single_flow(float* output, int32_t* Sreceivers, uint32_t* Stack, uint32_t* dim, float dx){

	uint32_t nxy = nxy(dim);
	const float cell_area = dx*dx;
	
	for(size_t i=0; i<nxy;++i){
		size_t ri = nxy - 1 - i;
		int32_t node = Stack[ri];
		if(node == Sreceivers[node])
			continue;
		output[node] += dx*dx;
		output[Sreceivers[node]] += output[node];
	}

}



