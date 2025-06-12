#define TOPOTOOLBOX_BUILD

#include <math.h>
#include <stddef.h>
#include <stdint.h>

#include "topotoolbox.h"

TOPOTOOLBOX_API
void streamquad_trapz_f32(float *integral, float *integrand, ptrdiff_t *source,
                          ptrdiff_t *target, float *weight,
                          ptrdiff_t edge_count) {
  // Iterate over the edges in reverse topological order (upstream)
  for (ptrdiff_t e = edge_count - 1; e >= 0; e--) {
    // NOTE: source and target are 1-based indices into the node attribute
    // list, so we must subtract 1 here
    ptrdiff_t u = source[e];
    ptrdiff_t v = target[e];
    integral[u] = integral[v] + weight[e] * (integrand[u] + integrand[v]) / 2;
  }
}

TOPOTOOLBOX_API
void streamquad_trapz_f64(double *integral, double *integrand,
                          ptrdiff_t *source, ptrdiff_t *target, float *weight,
                          ptrdiff_t edge_count) {
  // Iterate over the edges in reverse topological order (upstream)
  for (ptrdiff_t e = edge_count - 1; e >= 0; e--) {
    // NOTE: source and target are 1-based indices into the node attribute
    // list, so we must subtract 1 here
    ptrdiff_t u = source[e];
    ptrdiff_t v = target[e];
    integral[u] = integral[v] + weight[e] * (integrand[u] + integrand[v]) / 2;
  }
}

TOPOTOOLBOX_API
void traverse_up_u32_or_and(uint32_t *output, uint32_t *input,
                            ptrdiff_t *source, ptrdiff_t *target,
                            ptrdiff_t edge_count) {
  for (ptrdiff_t e = edge_count - 1; e >= 0; e--) {
    ptrdiff_t u = source[e];
    ptrdiff_t v = target[e];

    output[u] = output[u] | (output[v] & input[e]);
  }
}

TOPOTOOLBOX_API
void traverse_down_u32_or_and(uint32_t *output, uint32_t *input,
                              ptrdiff_t *source, ptrdiff_t *target,
                              ptrdiff_t edge_count) {
  for (ptrdiff_t e = 0; e < edge_count; e++) {
    ptrdiff_t u = source[e];
    ptrdiff_t v = target[e];

    output[v] = output[v] | (output[u] & input[e]);
  }
}

TOPOTOOLBOX_API
void traverse_down_f32_max_add(float *output, float *input, ptrdiff_t *source,
                               ptrdiff_t *target, ptrdiff_t edge_count) {
  for (ptrdiff_t e = 0; e < edge_count; e++) {
    ptrdiff_t u = source[e];
    ptrdiff_t v = target[e];

    output[v] = fmaxf(output[v], output[u] + input[e]);
  }
}

TOPOTOOLBOX_API
void traverse_up_f32_max_add(float *output, float *input, ptrdiff_t *source,
                             ptrdiff_t *target, ptrdiff_t edge_count) {
  for (ptrdiff_t e = edge_count - 1; e >= 0; e--) {
    ptrdiff_t u = source[e];
    ptrdiff_t v = target[e];

    output[u] = fmaxf(output[u], output[v] + input[e]);
  }
}

TOPOTOOLBOX_API
void traverse_down_f32_min_add(float *output, float *input, ptrdiff_t *source,
                               ptrdiff_t *target, ptrdiff_t edge_count) {
  for (ptrdiff_t e = 0; e < edge_count; e++) {
    ptrdiff_t u = source[e];
    ptrdiff_t v = target[e];

    output[v] = fminf(output[v], output[u] + input[e]);
  }
}

TOPOTOOLBOX_API
void traverse_down_f32_add_mul(float *output, float *input, ptrdiff_t *source,
                               ptrdiff_t *target, ptrdiff_t edge_count) {
  for (ptrdiff_t e = 0; e < edge_count; e++) {
    ptrdiff_t u = source[e];
    ptrdiff_t v = target[e];

    output[v] = output[v] + output[u] * input[e];
  }
}

TOPOTOOLBOX_API
void edgelist_degree(uint8_t *indegree, uint8_t *outdegree, ptrdiff_t *source,
                     ptrdiff_t *target, ptrdiff_t node_count,
                     ptrdiff_t edge_count) {
  for (ptrdiff_t v = 0; v < node_count; v++) {
    indegree[v] = 0;
    outdegree[v] = 0;
  }

  for (ptrdiff_t e = 0; e < edge_count; e++) {
    ptrdiff_t u = source[e];
    ptrdiff_t v = target[e];

    indegree[v]++;
    outdegree[u]++;
  }
}

/*
  propagatevaluesupstream is implemented identically for several
  different input data types. While it is possible to reduce this
  duplication using a macro, the various interfaces still need to be
  exposed to downstream users of the library.

  The following input data types are implemented

  float:    propagatevaluesupstream_f32
  double:   propagatevaluesupstream_f64
  uint8_t:  propagatevaluesupstream_u8
  uint32_t: propagatevaluesupstream_u32
  uint64_t: propagatevaluesupstream_u64
  int8_t:   propagatevaluesupstream_u8
  int32_t:  propagatevaluesupstream_u32
  int64_t:  propagatevaluesupstream_u64

  Others can be added as needed.
 */

TOPOTOOLBOX_API void propagatevaluesupstream_f32(float *data, ptrdiff_t *source,
                                                 ptrdiff_t *target,
                                                 ptrdiff_t edge_count) {
  for (ptrdiff_t e = edge_count - 1; e >= 0; e--) {
    ptrdiff_t u = source[e];
    ptrdiff_t v = target[e];
    data[u] = data[v];
  }
}

TOPOTOOLBOX_API void propagatevaluesupstream_f64(double *data,
                                                 ptrdiff_t *source,
                                                 ptrdiff_t *target,
                                                 ptrdiff_t edge_count) {
  for (ptrdiff_t e = edge_count - 1; e >= 0; e--) {
    ptrdiff_t u = source[e];
    ptrdiff_t v = target[e];
    data[u] = data[v];
  }
}

TOPOTOOLBOX_API void propagatevaluesupstream_u8(uint8_t *data,
                                                ptrdiff_t *source,
                                                ptrdiff_t *target,
                                                ptrdiff_t edge_count) {
  for (ptrdiff_t e = edge_count - 1; e >= 0; e--) {
    ptrdiff_t u = source[e];
    ptrdiff_t v = target[e];
    data[u] = data[v];
  }
}

TOPOTOOLBOX_API void propagatevaluesupstream_u32(uint32_t *data,
                                                 ptrdiff_t *source,
                                                 ptrdiff_t *target,
                                                 ptrdiff_t edge_count) {
  for (ptrdiff_t e = edge_count - 1; e >= 0; e--) {
    ptrdiff_t u = source[e];
    ptrdiff_t v = target[e];
    data[u] = data[v];
  }
}

TOPOTOOLBOX_API void propagatevaluesupstream_u64(uint64_t *data,
                                                 ptrdiff_t *source,
                                                 ptrdiff_t *target,
                                                 ptrdiff_t edge_count) {
  for (ptrdiff_t e = edge_count - 1; e >= 0; e--) {
    ptrdiff_t u = source[e];
    ptrdiff_t v = target[e];
    data[u] = data[v];
  }
}

TOPOTOOLBOX_API void propagatevaluesupstream_i8(int8_t *data, ptrdiff_t *source,
                                                ptrdiff_t *target,
                                                ptrdiff_t edge_count) {
  for (ptrdiff_t e = edge_count - 1; e >= 0; e--) {
    ptrdiff_t u = source[e];
    ptrdiff_t v = target[e];
    data[u] = data[v];
  }
}

TOPOTOOLBOX_API void propagatevaluesupstream_i32(int32_t *data,
                                                 ptrdiff_t *source,
                                                 ptrdiff_t *target,
                                                 ptrdiff_t edge_count) {
  for (ptrdiff_t e = edge_count - 1; e >= 0; e--) {
    ptrdiff_t u = source[e];
    ptrdiff_t v = target[e];
    data[u] = data[v];
  }
}

TOPOTOOLBOX_API void propagatevaluesupstream_i64(int64_t *data,
                                                 ptrdiff_t *source,
                                                 ptrdiff_t *target,
                                                 ptrdiff_t edge_count) {
  for (ptrdiff_t e = edge_count - 1; e >= 0; e--) {
    ptrdiff_t u = source[e];
    ptrdiff_t v = target[e];
    data[u] = data[v];
  }
}
