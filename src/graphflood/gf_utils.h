#pragma once
#define TOPOTOOLBOX_BUILD

#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include "../morphology/reconstruct.h"
// #include "topotoolbox.h"


void generate_offset_D4(int8_t** off0, int8_t** off1);
void generate_offset_D8(int8_t** off0, int8_t** off1);
void generate_offset_D4_flat(int32_t** off, uint32_t* dim);
void generate_offset_D8_flat(int32_t** off, uint32_t* dim);
void generate_offsetdx_D4(float** off, float dx);
void generate_offsetdx_D8(float** off, float dx);

// int32_t dim2flat(int32_t this_dim1_dim0, int32_t this_dim1, uint32_t* dim);
int32_t dim2flat(uint32_t this_dim0, uint32_t this_dim1, uint32_t* dim);
void flat2dim(int32_t node, uint32_t* this_dim0, uint32_t* this_dim1, uint32_t* dim);
bool can_receive(int32_t node, uint8_t* BCs);
bool can_give(int32_t node, uint8_t* BCs);
bool is_nodata(int32_t node, uint8_t* BCs);
uint8_t N_neighbour(bool D8);
