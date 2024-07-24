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
void generate_offset_D4(uint8_t** off0, uint8_t** off1){
	*off0 = uint8_t[4]{-1,0,0,1};
	*off1 = uint8_t[4]{0,-1,1,0};
}

/*
Generate the offsets for neighbouring operations 
D8 version:

Row major representation
|0|1|2|
|3|x|4|
|5|6|7|
*/
void generate_offset_D8(uint8_t** off0, uint8_t** off1){
	*off0 = uint8_t[8]{-1,-1,-1,0,0,1,1,1};
	*off1 = uint8_t[8]{-1,0,1,-1,1,-1,0,1};
}


/*
returns the flat index of a node from its dimensions

Example in row major:
this_dim0 is the current row, dim[1] is the number of columns and this_dim1 is the current column

*/
int32_t dim2flat(int32_t this_dim0, int32_t this_dim1, int32_t* dim){
	return this_dim0 * dim[1] + this_dim1;
}


void flat2dim(int32_t* this_dim0, int32_t* this_dim1, int32_t* dim){

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


