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
void generate_offset_D4(int8_t* off0, int8_t* off1) {
  // creating a generic array to assign to the pointer
  off0[0] = -1;
  off0[1] = 0;
  off0[2] = 0;
  off0[3] = -1;
  // creating a generic array to assign to the pointer
  off1[0] = 0;
  off1[1] = -1;
  off1[2] = 1;
  off1[3] = 0;
}

/*
Generate the offseting distance for neighbouring operations
D4 version.
Row major representation
| |0| |
|1|x|2|
| |3| |
*/
void generate_offsetdx_D4(GF_FLOAT* off, GF_FLOAT dx) {
  // creating a generic array to assign to the pointer
  off[0] = dx;
  off[1] = dx;
  off[2] = dx;
  off[3] = dx;
}

/*
Generate the offsets for neighbouring operations
D8 version:

Row major representation
|0|1|2|
|3|x|4|
|5|6|7|
*/
void generate_offset_D8(int8_t* off0, int8_t* off1) {
  // creating a generic array to assign to the pointer
  off0[0] = -1;
  off0[1] = -1;
  off0[2] = -1;
  off0[3] = 0;
  off0[4] = 0;
  off0[5] = 1;
  off0[6] = 1;
  off0[7] = 1;
  // creating a generic array to assign to the pointer
  off1[0] = -1;
  off1[1] = 0;
  off1[2] = 1;
  off1[3] = -1;
  off1[4] = 1;
  off1[5] = -1;
  off1[6] = 0;
  off1[7] = 1;
}

/*
Generate the offseting distance for neighbouring operations
D8 version.
Row major representation
| |0| |
|1|x|2|
| |3| |
*/
void generate_offsetdx_D8(GF_FLOAT* off, GF_FLOAT dx) {
  GF_FLOAT diag = (GF_FLOAT)sqrt(2) * dx;
  // creating a generic array to assign to the pointer
  off[0] = diag;
  off[1] = dx;
  off[2] = diag;
  off[3] = dx;
  off[4] = dx;
  off[5] = diag;
  off[6] = dx;
  off[7] = diag;
}

/*
Generate the flat offsets for neighbouring operations
D4 version.
*/
void generate_offset_D4_flat(GF_INT* off, GF_UINT* dim) {
  // creating a generic array to assign to the pointer
  off[0] = -(GF_INT)dim[1];
  off[1] = -1;
  off[2] = 1;
  off[3] = (GF_INT)dim[1];
}

/*
Generate the flat offsets for neighbouring operations
D8 version.
*/
void generate_offset_D8_flat(GF_INT* off, GF_UINT* dim) {
  // creating a generic array to assign to the pointer
  off[0] = -(GF_INT)dim[1] - 1;
  off[1] = -(GF_INT)dim[1] + 0;
  off[2] = -(GF_INT)dim[1] + 1;
  off[3] = -1;
  off[4] = 1;
  off[5] = dim[1] - 1;
  off[6] = dim[1] + 0;
  off[7] = dim[1] + 1;
}

/*
returns the flat index of a node from its dimensions

Example in row major:
this_dim0 is the current row, dim[1] is the number of columns and this_dim1 is
the current column

*/
// GF_INT dim2flat(GF_INT this_dim0, GF_INT this_dim1, GF_UINT* dim){
// 	return this_dim0 * dim[1] + this_dim1;
// }
GF_UINT dim2flat(GF_UINT this_dim0, GF_UINT this_dim1, GF_UINT* dim) {
  return (this_dim0 * dim[1] + this_dim1);
}

void flat2dim(GF_UINT node, GF_UINT* this_dim0, GF_UINT* this_dim1,
              GF_UINT* dim) {
  *this_dim0 = node % dim[1];
  *this_dim1 = node - (*this_dim0) * dim[1];
}

/*
Return the number of neighbours
depending on the DX topology
*/
uint8_t N_neighbour(bool D8) {
  if (D8)
    return 8;
  else
    return 4;
}

GF_UINT nxy(GF_UINT* dim) { return dim[0] * dim[1]; }

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
For example, a node at the edge - but that have a neighbour with lower
elevation: CAN_OUT = 3,

Flow can only out when entering the cell:
OUT = 4,

PROBABLY LEGACY
FORCE_OUT = 5,

For cells located at the edge of a DEM, flow can pass through but not leave
(local minima if no downstream neighbour) CANNOT_OUT = 6,

Flow can only flow to potential receivers
exemple, you have a DEM of a river segment and want to input water from a side:
IN = 7,

Same than 7, but cannot even receive from upstream neighbours
FORCE_IN = 8,

periodic border (not implemented yet but will have to be):
PERIODIC_BORDER = 9


uint8_t has 256 possibilities so there is plenty of space for more options
*/

/*
I know BCs is not used so far here but that's future proof for periodic
conditions or other edge cases
*/
void check_bound(GF_UINT node, GF_UINT* dim, uint8_t* BCs, bool* valid) {
  if (node >= dim[0] * dim[1]) *valid = false;
}

void check_top_customs(GF_UINT node, uint8_t n, GF_UINT* dim, uint8_t* BCs,
                       bool* valid, bool D8) {
  // Checking if the node is on the first row and trying to target above
  if (node < dim[1]) {
    if ((D8 && n < 3) || (D8 == false && n == 0)) *valid = false;
  }
  // # Done
}

void check_bottom_customs(GF_UINT node, uint8_t n, GF_UINT* dim, uint8_t* BCs,
                          bool* valid, bool D8) {
  // Checking if the node is on the first row and trying to target above
  if (node >= dim[0] * dim[1] - dim[1]) {
    if ((D8 && n >= 5) || (D8 == false && n == 3)) *valid = false;
  }
  // # Done
}

void check_left_customs(GF_UINT node, uint8_t n, GF_UINT* dim, uint8_t* BCs,
                        bool* valid, bool D8) {
  // Checking if the node is on the first row and trying to target above
  if (node % dim[1] == 0) {
    if ((D8 && (n == 0 || n == 3 || n == 5)) || (D8 == false && n == 1))
      *valid = false;
  }
  // # Done
}

void check_right_customs(GF_UINT node, uint8_t n, GF_UINT* dim, uint8_t* BCs,
                         bool* valid, bool D8) {
  // Checking if the node is on the first row and trying to target above
  if (node % dim[1] == dim[1] - 1) {
    if ((D8 && (n == 2 || n == 4 || n == 7)) || (D8 == false && n == 2))
      *valid = false;
  }
  // # Done
}

bool check_bound_neighbour(GF_UINT node, uint8_t n, GF_UINT* dim, uint8_t* BCs,
                           bool D8) {
  bool valid = true;
  check_bound(node, dim, BCs, &valid);
  check_top_customs(node, n, dim, BCs, &valid, D8);
  check_bottom_customs(node, n, dim, BCs, &valid, D8);
  check_left_customs(node, n, dim, BCs, &valid, D8);
  check_right_customs(node, n, dim, BCs, &valid, D8);
  return valid;
}

bool can_receive(GF_UINT node, uint8_t* BCs) {
  if (BCs[node] == 1 || BCs[node] == 3 || BCs[node] == 4 || BCs[node] == 5 ||
      BCs[node] == 6 || BCs[node] == 7 || BCs[node] == 9

  ) {
    return true;
  } else {
    return false;
  }
}

bool can_give(GF_UINT node, uint8_t* BCs) {
  if (BCs[node] == 1 || BCs[node] == 3 || BCs[node] == 6 || BCs[node] == 7 ||
      BCs[node] == 8 || BCs[node] == 9

  ) {
    return true;
  } else {
    return false;
  }
}

bool can_out(GF_UINT node, uint8_t* BCs) {
  if (BCs[node] == 3 || BCs[node] == 4 || BCs[node] == 5

  ) {
    return true;
  } else {
    return false;
  }
}

bool is_nodata(GF_UINT node, uint8_t* BCs) {
  if (BCs[node] == 0) {
    return true;
  } else {
    return false;
  }
}
