#define TOPOTOOLBOX_BUILD

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
    ptrdiff_t u = source[e] - 1;
    ptrdiff_t v = target[e] - 1;
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
    ptrdiff_t u = source[e] - 1;
    ptrdiff_t v = target[e] - 1;
    integral[u] = integral[v] + weight[e] * (integrand[u] + integrand[v]) / 2;
  }
}
