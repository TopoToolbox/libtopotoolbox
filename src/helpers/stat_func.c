#include "stat_func.h"

#include <float.h>
#include <math.h>
#include <stdlib.h>

// ============================================================================
// stats_accumulator
// ============================================================================

void accumulator_init(stats_accumulator *acc) {
  acc->count = 0;
  acc->sum = 0.0;
  acc->sum_sq = 0.0;
  acc->min_val = FLT_MAX;
  acc->max_val = -FLT_MAX;
}

void accumulator_add(stats_accumulator *acc, float value) {
  if (isnan(value)) return;
  acc->count++;
  acc->sum += value;
  acc->sum_sq += value * value;
  if (value < acc->min_val) acc->min_val = value;
  if (value > acc->max_val) acc->max_val = value;
}

float accumulator_mean(const stats_accumulator *acc) {
  return acc->count > 0 ? (float)(acc->sum / acc->count) : NAN;
}

float accumulator_stddev(const stats_accumulator *acc) {
  if (acc->count <= 1) return 0.0f;
  double variance =
      (acc->sum_sq - acc->sum * acc->sum / acc->count) / (acc->count - 1);
  return variance > 0.0 ? (float)sqrt(variance) : 0.0f;
}

// ============================================================================
// percentile_accumulator
// ============================================================================

static int compare_floats(const void *a, const void *b) {
  float fa = *(const float *)a;
  float fb = *(const float *)b;
  return (fa > fb) - (fa < fb);
}

void percentile_accumulator_init(percentile_accumulator *acc,
                                 ptrdiff_t initial_capacity) {
  acc->capacity = initial_capacity > 0 ? initial_capacity : 16;
  acc->count = 0;
  acc->values = (float *)malloc(acc->capacity * sizeof(float));
}

void percentile_accumulator_add(percentile_accumulator *acc, float value) {
  if (isnan(value)) return;
  if (acc->count >= acc->capacity) {
    acc->capacity *= 2;
    acc->values = (float *)realloc(acc->values, acc->capacity * sizeof(float));
  }
  acc->values[acc->count++] = value;
}

void percentile_accumulator_sort(percentile_accumulator *acc) {
  if (acc->count > 0) {
    qsort(acc->values, acc->count, sizeof(float), compare_floats);
  }
}

float percentile_accumulator_get(const percentile_accumulator *acc, float p) {
  if (acc->count == 0) return NAN;
  if (p < 0.0f) p = 0.0f;
  if (p > 100.0f) p = 100.0f;
  float index = p / 100.0f * (float)(acc->count - 1);
  ptrdiff_t lower = (ptrdiff_t)index;
  ptrdiff_t upper = lower + 1;
  if (upper >= acc->count) return acc->values[acc->count - 1];
  float fraction = index - (float)lower;
  return acc->values[lower] * (1.0f - fraction) + acc->values[upper] * fraction;
}

void percentile_accumulator_free(percentile_accumulator *acc) {
  free(acc->values);
  acc->values = NULL;
  acc->count = 0;
  acc->capacity = 0;
}
