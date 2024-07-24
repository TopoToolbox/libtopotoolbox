#pragma once
#define TOPOTOOLBOX_BUILD

#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <stdbool.h>

#include "../morphology/reconstruct.h"
#include "topotoolbox.h"


void generate_offset_D4(uint8_t** off0, uint8_t** off1);
void generate_offset_D8(uint8_t** off0, uint8_t** off1);
int32_t dim2flat(int32_t this_dim0, int32_t this_dim1, int32_t* dim);
void flat2dim(int32_t* this_dim0, int32_t* this_dim1, int32_t* dim);

