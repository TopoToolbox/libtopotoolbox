#pragma once
#define TOPOTOOLBOX_BUILD

#include <math.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "graphflood/define_types.h"
// #include "topotoolbox.h"

/*
Offset helpers:
This bunch of functions generate offset arrays used in neighbouring operations.
It creates a generic framework to iterate in D4 or D8 generically and
future-proofed.

Typically, these are called at the top of a function (cheap to generate if done
once per function)

then they allow a generic neighbouring operation:
for(uint8_t n=0; n<N_neighbour(D8); ++n){
        GF_UINT neighbour = node + offset[n];
        GF_FLOAT dx = offsetdx[n];
}
*/

/*
@Brief Helper function to create the offset arrays for 2D neighbouring in D4
(cardinal only) To be called with int8_t off[N_neighbour(D8)] to feed them;
*/
void generate_offset_D4(int8_t* off0, int8_t* off1);

/*
@Brief Helper function to create the offset arrays for 2D neighbouring in D8
(cardinal + diagonals) To be called with int8_t off[N_neighbour(D8)] to feed
them;
*/
void generate_offset_D8(int8_t* off0, int8_t* off1);

/*
@Brief Helper function to create the offset arrays for 1D neighbouring in D4
(cardinal only) 1D neighbouring uses flat indices: ID = row * ncolumns + column;
for ROW MAJOR or ID = column * nrows + row; for COLUMN MAJOR

@param[out] off: a int32_t[4/8] array to be fed in place;
@param[in] dim: the dimension array (nx,ny or ny,nx depending on row/col major)
*/
void generate_offset_D4_flat(GF_INT* off, GF_UINT* dim);

/*
@Brief Helper function to create the offset arrays for 1D neighbouring in D8
(cardinal only) 1D neighbouring uses flat indices: ID = row * ncolumns + column;
for ROW MAJOR or ID = column * nrows + row; for COLUMN MAJOR

@param[out] off: a GF_UINT[4/8] array to be fed in place;
@param[in] dim: the dimension array (nx,ny or ny,nx depending on row/col major)
*/
void generate_offset_D8_flat(GF_INT* off, GF_UINT* dim);

/*
@brief offset distance for neighbouring operations in D4
Note that all the values are dx, but this is for consistency with D8
*/
void generate_offsetdx_D4(GF_FLOAT* off, GF_FLOAT dx);

/*
@brief offset distance for neighbouring operations in D8
*/
void generate_offsetdx_D8(GF_FLOAT* off, GF_FLOAT dx);

/*
Helper functions for geometry and indices
Note: flat index is a vectorised 1D index to represent row, col or col, row

*/

/*
@Brief Takes the indices ofdim[0] and dim[1] [row,col (row major) or col row
(column major)] and returns the flat index
@param[in] this_dim0: the index in dim[0]
@param[in] this_dim1: the index in dim[1]
@param[in] dim: the dimensions

@return the flat 1D index
*/
GF_UINT dim2flat(GF_UINT this_dim0, GF_UINT this_dim1, GF_UINT* dim);

/*
@Brief Takes the indices ofdim[0] and dim[1] [row,col (row major) or col row
(column major)] and returns the flat index

@param[in] node: the flat index
@param[out] this_dim0: the index in dim[0]
@param[out] this_dim1: the index in dim[1]
@param[in] dim: the dimensions

*/
void flat2dim(GF_UINT node, GF_UINT* this_dim0, GF_UINT* this_dim1,
              GF_UINT* dim);
bool can_receive(GF_UINT node, uint8_t* BCs);
bool can_give(GF_UINT node, uint8_t* BCs);
bool can_out(GF_UINT node, uint8_t* BCs);
bool is_nodata(GF_UINT node, uint8_t* BCs);
uint8_t N_neighbour(bool D8);
GF_UINT nxy(GF_UINT* dim);
bool check_bound_neighbour(GF_UINT node, uint8_t n, GF_UINT* dim, uint8_t* BCs,
                           bool D8);
