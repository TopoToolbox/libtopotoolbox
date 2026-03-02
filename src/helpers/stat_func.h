#pragma once

#include <stddef.h>

// ============================================================================
// Online statistics accumulator — running mean, variance, min, max (O(1))
// ============================================================================

typedef struct {
  ptrdiff_t count;
  double sum;
  double sum_sq;
  float min_val;
  float max_val;
} stats_accumulator;

void accumulator_init(stats_accumulator *acc);
void accumulator_add(stats_accumulator *acc, float value);
float accumulator_mean(const stats_accumulator *acc);
float accumulator_stddev(const stats_accumulator *acc);

// ============================================================================
// Percentile accumulator — collects values, sorts on demand, queries by rank
// ============================================================================

typedef struct {
  ptrdiff_t capacity;
  ptrdiff_t count;
  float *values;
} percentile_accumulator;

void percentile_accumulator_init(percentile_accumulator *acc,
                                 ptrdiff_t initial_capacity);
void percentile_accumulator_add(percentile_accumulator *acc, float value);
void percentile_accumulator_sort(percentile_accumulator *acc);
float percentile_accumulator_get(const percentile_accumulator *acc, float p);
void percentile_accumulator_free(percentile_accumulator *acc);
