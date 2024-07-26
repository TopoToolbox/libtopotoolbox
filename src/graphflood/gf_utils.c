#define TOPOTOOLBOX_BUILD
#include "gf_utils.h"

/*
Generate the offsets for neighbouring operations
D4 version.

Row major representation
| |0| |
|1|x|2|
| |3| |
*/
void generate_offset_D4(int8_t** off0, int8_t** off1){
	*off0 = malloc( 4 * sizeof (int8_t));
	*off0[0] = -1;
	*off0[1] = 0;
	*off0[2] = 0;
	*off0[3] = -1;
	// *off0 = (*int8_t[4]){-1,0,0,1};
	*off1 = malloc( 4 * sizeof (int8_t));
	*off1[0] = 0;
	*off1[1] = -1;
	*off1[2] = 1;
	*off1[3] = 0;
	// *off1 = (*int8_t[4]){0,-1,1,0};
}

/*
Generate the offseting distance for neighbouring operations
D4 version.
Row major representation
| |0| |
|1|x|2|
| |3| |
*/
void generate_offsetdx_D4(float** off, float dx){
	*off = malloc( 4 * sizeof (float));
	*off[0] = dx;
	*off[1] = dx;
	*off[2] = dx;
	*off[3] = dx;
}

/*
Generate the offsets for neighbouring operations 
D8 version:

Row major representation
|0|1|2|
|3|x|4|
|5|6|7|
*/
void generate_offset_D8(int8_t** off0, int8_t** off1){
	*off0 = malloc( 8 * sizeof (int8_t));
	*off0[0] = -1;
	*off0[1] = -1;
	*off0[2] = -1;
	*off0[3] = 0 ;
	*off0[4] = 0 ;
	*off0[5] = 1 ;
	*off0[6] = 1 ;
	*off0[7] = 1 ;
	// *off0 = (*int8_t[8]){-1,-1,-1,0,0,1,1,1};
	*off1 = malloc( 8 * sizeof (int8_t));
	*off1[0] = -1;
	*off1[1] = 0 ;
	*off1[2] = 1 ;
	*off1[3] = -1;
	*off1[4] = 1 ;
	*off1[5] = -1;
	*off1[6] = 0 ;
	*off1[7] = 1 ;
	// *off1 = (*int8_t[8]){-1,0,1,-1,1,-1,0,1};
}

/*
Generate the offseting distance for neighbouring operations
D8 version.
Row major representation
| |0| |
|1|x|2|
| |3| |
*/
void generate_offsetdx_D8(float** off, float dx){
	float diag = sqrt(2) * dx;
	*off = malloc( 4 * sizeof (float));
	*off[0] = diag;
	*off[1] = dx;
	*off[2] = diag;
	*off[3] = dx;
	*off[4] = dx;
	*off[5] = diag;
	*off[6] = dx;
	*off[7] = diag;
}

/*
Generate the flat offsets for neighbouring operations
D4 version.
*/
void generate_offset_D4_flat(int32_t** off, uint32_t* dim){
	*off = malloc( 4 * sizeof (int32_t));
	*off[0] = -dim[0];
	*off[1] = -1;
	*off[2] = 1;
	*off[3] = dim[0];
}

/*
Generate the flat offsets for neighbouring operations
D8 version.
*/
void generate_offset_D8_flat(int32_t** off, uint32_t* dim){
	*off = malloc( 8 * sizeof (int32_t));
	*off[0] = -dim[0] - 1;
	*off[1] = -dim[0] + 0;
	*off[2] = -dim[0] + 1;
	*off[3] = -1;
	*off[4] = 1 ;
	*off[5] = dim[0] - 1;
	*off[6] = dim[0] + 0 ;
	*off[7] = dim[0] + 1 ;
}

/*
Generate the offsets for neighbouring operations 
D8 version:

Row major representation
|0|1|2|
|3|x|4|
|5|6|7|
*/
void generate_offset_D8(int8_t** off0, int8_t** off1){
	*off0 = malloc( 8 * sizeof (int8_t));
	*off0[0] = -1;
	*off0[1] = -1;
	*off0[2] = -1;
	*off0[3] = 0 ;
	*off0[4] = 0 ;
	*off0[5] = 1 ;
	*off0[6] = 1 ;
	*off0[7] = 1 ;
	// *off0 = (*int8_t[8]){-1,-1,-1,0,0,1,1,1};
	*off1 = malloc( 8 * sizeof (int8_t));
	*off1[0] = -1;
	*off1[1] = 0 ;
	*off1[2] = 1 ;
	*off1[3] = -1;
	*off1[4] = 1 ;
	*off1[5] = -1;
	*off1[6] = 0 ;
	*off1[7] = 1 ;
	// *off1 = (*int8_t[8]){-1,0,1,-1,1,-1,0,1};
}


/*
returns the flat index of a node from its dimensions

Example in row major:
this_dim0 is the current row, dim[1] is the number of columns and this_dim1 is the current column

*/
int32_t dim2flat(int32_t this_dim0, int32_t this_dim1, uint32_t* dim){
	return this_dim0 * dim[1] + this_dim1;
}
int32_t dim2flat(uint32_t this_dim0, uint32_t this_dim1, uint32_t* dim){
	return this_dim0 * dim[1] + this_dim1;
}


void flat2dim(int32_t node, uint32_t* this_dim0, uint32_t* this_dim1, uint32_t* dim){
	*this_dim0 = node % dim[1];
	*this_dim1 = node - (*this_dim0) * dim[1];
}

/*
Return the number of neighbours
depending on the DX topology
*/
uint8_t N_neighbour(bool D8){
	if(D8)
		return 8;
	else
		return 4;
}



/*
The following functions helps determining how a given node handles flux.
It uses the standard used in DAGGER/scabbard tools
Each node has a uint8_t code describing its status.
Note that capital letters correspond to the c++ enum class in DAGGER.

The tool is WIP so some of these might not be implemented in every tools yet.

Cannot flow at all (nodata):
NO_FLOW = 0,

Internal Node (can flow in and out of the cell):
FLOW = 1,

2 is legacy and does not do anything

Flow can out there but can also flow to downstream neighbours.
For example, a node at the edge - but that have a neighbour with lower elevation:
CAN_OUT = 3,

Flow can only out when entering the cell:
OUT = 4,

PROBABLY LEGACY
FORCE_OUT = 5,

For cells located at the edge of a DEM, flow can pass through but not leave (local minima if no downstream neighbour)
CANNOT_OUT = 6,

Flow can only flow to potential receivers
exemple, you have a DEM of a river segment and want to input water from a side:
IN = 7,

Same than 7, but cannot even receive from upstream neighbours
FORCE_IN = 8,

periodic border (not implemented yet but will have to be):
PERIODIC_BORDER = 9


uint8_t has 256 possibilities so there is plenty of space for more options
*/





bool can_receive(int32_t node, uint8_t* BCs){
	if(
		BCs[node] == 1 ||
		BCs[node] == 3 ||
		BCs[node] == 4 ||
		BCs[node] == 5 ||
		BCs[node] == 6 ||
		BCs[node] == 7 ||
		BCs[node] == 9

	){
		return true;
	}else{
		return false;
	}
}

bool can_give(int32_t node, uint8_t* BCs){
	if(
		BCs[node] == 1 ||
		BCs[node] == 3 ||
		BCs[node] == 6 ||
		BCs[node] == 7 ||
		BCs[node] == 8 ||
		BCs[node] == 9

	){
		return true;
	}else{
		return false;
	}
}


bool is_nodata(int32_t node, uint8_t* BCs){
	if(
		BCs[node] == 0
	){
		return true;
	}else{
		return false;
	}
}
