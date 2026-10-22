// Tests of traverse_down_f32_add and traverse_down_f64_add.
//
// A network drawn by hand with known accumulations, and random
// single-flow-direction networks on which the add-only traversal must
// equal traverse_down_*_add_mul with a fraction of ones, bit for bit.
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "topotoolbox.h"

// Nodes 0 and 1 drain to 2, nodes 2 and 3 drain to 4, node 4 drains to
// 5. The edges are listed with their sources in topological order.
static int test_hand_network(void) {
  ptrdiff_t source[5] = {0, 1, 2, 3, 4};
  ptrdiff_t target[5] = {2, 2, 4, 4, 5};
  float expected[6] = {1.0f, 1.0f, 3.0f, 1.0f, 5.0f, 6.0f};

  float output32[6] = {1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f};
  traverse_down_f32_add(output32, source, target, 5);

  double output64[6] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
  traverse_down_f64_add(output64, source, target, 5);

  for (ptrdiff_t i = 0; i < 6; i++) {
    if (output32[i] != expected[i]) {
      printf("traverse_down_f32_add: node %lld is %g, expected %g\n",
             (long long)i, (double)output32[i], (double)expected[i]);
      return 1;
    }
    if (output64[i] != (double)expected[i]) {
      printf("traverse_down_f64_add: node %lld is %g, expected %g\n",
             (long long)i, output64[i], (double)expected[i]);
      return 1;
    }
  }
  return 0;
}

// xorshift64, so that every platform draws the same networks
static uint64_t next_random(uint64_t *state) {
  uint64_t x = *state;
  x ^= x << 13;
  x ^= x >> 7;
  x ^= x << 17;
  *state = x;
  return x;
}

// Every node except the last drains to a random node with a larger
// index, so that the edges, listed by source, are in topological
// order, every node has exactly one outgoing edge, and the last node
// is the only outlet.
static int test_random_network(uint64_t seed, ptrdiff_t node_count) {
  ptrdiff_t edge_count = node_count - 1;
  ptrdiff_t *source = malloc(edge_count * sizeof(ptrdiff_t));
  ptrdiff_t *target = malloc(edge_count * sizeof(ptrdiff_t));
  float *fraction32 = malloc(edge_count * sizeof(float));
  double *fraction64 = malloc(edge_count * sizeof(double));
  float *weights = malloc(node_count * sizeof(float));
  float *add32 = malloc(node_count * sizeof(float));
  float *mul32 = malloc(node_count * sizeof(float));
  double *add64 = malloc(node_count * sizeof(double));
  double *mul64 = malloc(node_count * sizeof(double));

  int result = 0;
  if (!source || !target || !fraction32 || !fraction64 || !weights || !add32 ||
      !mul32 || !add64 || !mul64) {
    printf("allocation failed\n");
    result = 1;
  } else {
    uint64_t state = seed * 2654435761u + 1;
    for (ptrdiff_t u = 0; u < edge_count; u++) {
      source[u] = u;
      target[u] =
          u + 1 +
          (ptrdiff_t)(next_random(&state) % (uint64_t)(node_count - 1 - u));
      fraction32[u] = 1.0f;
      fraction64[u] = 1.0;
    }
    // node weights in [0, 1), so that the test is not only about ones
    double sum_of_weights = 0.0;
    for (ptrdiff_t v = 0; v < node_count; v++) {
      weights[v] = (float)(next_random(&state) % 1000000) / 1000000.0f;
      add32[v] = weights[v];
      mul32[v] = weights[v];
      add64[v] = (double)weights[v];
      mul64[v] = (double)weights[v];
      sum_of_weights += (double)weights[v];
    }

    traverse_down_f32_add(add32, source, target, edge_count);
    traverse_down_f32_add_mul(mul32, fraction32, source, target, edge_count);
    traverse_down_f64_add(add64, source, target, edge_count);
    traverse_down_f64_add_mul(mul64, fraction64, source, target, edge_count);

    if (memcmp(add32, mul32, node_count * sizeof(float)) != 0) {
      printf(
          "seed %llu: traverse_down_f32_add differs from "
          "traverse_down_f32_add_mul with a fraction of ones\n",
          (unsigned long long)seed);
      result = 1;
    }
    if (memcmp(add64, mul64, node_count * sizeof(double)) != 0) {
      printf(
          "seed %llu: traverse_down_f64_add differs from "
          "traverse_down_f64_add_mul with a fraction of ones\n",
          (unsigned long long)seed);
      result = 1;
    }
    // the float accumulation stays close to the double one
    for (ptrdiff_t v = 0; v < node_count; v++) {
      if (fabs((double)add32[v] - add64[v]) > 1e-4 * add64[v]) {
        printf("seed %llu: node %lld accumulates %g in float, %g in double\n",
               (unsigned long long)seed, (long long)v, (double)add32[v],
               add64[v]);
        result = 1;
        break;
      }
    }
    // the outlet drains every weight
    double outlet = add64[node_count - 1];
    if (fabs(outlet - sum_of_weights) > 1e-9 * sum_of_weights) {
      printf(
          "seed %llu: the outlet accumulates %.12g, the weights sum to "
          "%.12g\n",
          (unsigned long long)seed, outlet, sum_of_weights);
      result = 1;
    }
  }

  free(source);
  free(target);
  free(fraction32);
  free(fraction64);
  free(weights);
  free(add32);
  free(mul32);
  free(add64);
  free(mul64);
  return result;
}

int main(void) {
  int result = test_hand_network();
  for (uint64_t seed = 1; seed <= 20 && result == 0; seed++) {
    result = test_random_network(seed, 1000 * (ptrdiff_t)seed);
  }
  return result;
}
