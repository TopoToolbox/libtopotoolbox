#define TOPOTOOLBOX_BUILD

#include <math.h>
#include <stddef.h>
#include <stdint.h>

#include "priority_queue.h"
#include "topotoolbox.h"

#define SQRT2f 1.41421356237309504880f

// Gray-weighted distance transforms for auxiliary topography

/*
  # Cost computation

  The costs for the gray-weighted distance transform are given by

  cost[i,j] = (maxdepth - (filled_dem[i,j] - original_dem[i,j]))^tweight +
  CarveMinVal;

  where maxdepth is the maximum difference between the filled and
  original DEMs over the connected component of flat pixels to which
  the (i,j) pixel belongs.

  Computing the maximum depth is done with a connected components
  labeling algorithm by Wu et al. (2009) and described by He et
  al. (2017) and related to the Hoshen-Kopelman algorithm (1976). This
  algorithm takes two passes over the data. In the first pass, each
  flat pixel is assigned a provisional label. The labels of
  neighboring flat pixels are considered /equivalent/ because they
  belong to the same connected component. An /equivalence class/ is a
  set of equivalent labels that are all equivalent to each
  other. During the first pass, the algorithm builds a union-find data
  structure that identifies which equivalence class each pixel belongs
  to. This data structure is described below. The second pass resolves
  equivalent labels to a unique label per connected component, with
  that label determined by the pixel with the highest difference
  between filled and original DEMs. A third pass over the data
  computes the costs using the maximum depth for each connected
  component.

  # Union-find implementation

  Equivalence classes of labels are recorded in a union-find data
  structure. To avoid allocating additional memory for this data
  structure, it is built into the return values of the `compute_costs`
  function, namely a ptrdiff_t array of size (dims[0],dims[1]), which
  holds labels for the connected components, and a float array of the
  same size. The float array ultimately holds the costs returned from
  the function, but during the algorithm, it is used to keep track
  of the maximum depth over each connected component.

  The union-find data structure uses linear indices of pixels as
  labels of the connected components. A pixel in the labels array
  points to another pixel in its equivalence class. If it points to
  itself, it is called the /root/ of that equivalence class. While
  scanning through the DEM the first time, the algorithm inserts
  pixels into the union-find data structure by making them the root of
  their own equivalence class. It then /unifies/ the labels of
  neighboring flat pixels by making them point to the same root. The
  new root of two neighboring pixels is chosen so that it has the
  greatest DEM difference of the flat pixels in its equivalence
  class. The first pass determines the correct equivalence classes of
  each flat pixel. Every chain of labels eventually leads to a root
  pixel with the greatest DEM difference in its connected
  components. However, different pixels in the same connected
  component may still be labeled with different labels. A second pass
  is required to label each pixel with the root of its equivalence
  class and the corresponding maximum DEM difference.

  # References

  He, L., Ren, X., Gao, Q., Zhao, X., Yao, B., & Chao, Y. (2017). The
  connected-component labeling problem: A review of state-of-the-art
  algorithms. Pattern Recognition, 70, 25-43.

  Hoshen, J., & Kopelman, R. (1976). Percolation and cluster
  distribution. I. Cluster multiple labeling technique and critical
  concentration algorithm. Physical Review B, 14(8), 3438.

  Wu, K., Otoo, E., & Suzuki, K. (2009). Optimizing two-pass
  connected-component labeling algorithms. Pattern Analysis and
  Applications, 12, 117-135.
*/

// Insert pixel into the union-find data structure with a given
// weight. The pixel is the root of a new disjoint set, so it points
// to itself in the labels array.
static void new_tree(ptrdiff_t *labels, float *weights, ptrdiff_t pixel,
                     float weight) {
  labels[pixel] = pixel;
  weights[pixel] = weight;
}

// This struct lets us return a connected component label and the
// corresponding depth.
typedef struct {
  ptrdiff_t label;
  float weight;
} labeldepth;

// To find the root of the disjoint set containing `pixel`, follow the
// pointers in the `labels` array until reaching a root, where the
// pointer points to itself. Return the label of the root as well as
// its corresponding weight.
static labeldepth find_root(ptrdiff_t *labels, float *weights,
                            ptrdiff_t pixel) {
  labeldepth r = {pixel, 0.0};
  while (labels[r.label] != r.label) {
    r.label = labels[r.label];
  }
  r.weight = weights[r.label];
  return r;
}

// The trees in the union-find structure can get very tall, requiring
// many indirections to find the root of a given pixel. One way of
// managing the tree height is /path compression/, which traverses the
// tree starting from `pixel` and sets all of the labels and weights
// and the weights along the way to `root` and `weight`. The next time
// any of the pixels on the path is accessed, only one indirection is
// required to find the root.
static void set_root(ptrdiff_t *labels, float *weights, ptrdiff_t pixel,
                     ptrdiff_t root, float weight) {
  while (labels[pixel] != pixel) {
    ptrdiff_t v = labels[pixel];
    labels[pixel] = root;
    weights[pixel] = weight;
    pixel = v;
  }
  labels[pixel] = root;
  weights[pixel] = weight;
}

// The heart of the union-find data structure /unifies/ two pixels by
// considering them part of the same /equivalence class/. This is done
// by setting the roots of each pixel to point to the same pixel. The
// pixel that is chosen to be the new root is that with the greatest
// weight.
static labeldepth unify(ptrdiff_t *labels, float *weights, ptrdiff_t pixel1,
                        ptrdiff_t pixel2) {
  labeldepth r = find_root(labels, weights, pixel1);
  if (r.label != pixel2) {
    labeldepth s = find_root(labels, weights, pixel2);
    if (r.weight < s.weight) {
      r = s;
    }
    set_root(labels, weights, pixel2, r.label, r.weight);
  }
  set_root(labels, weights, pixel1, r.label, r.weight);
  return r;
}

/*
  Cost computation

  `costs`, `original_dem` and `filled_dem` should all represent
  two-dimensional single-precision floating point arrays of size
  (dims[0], dims[1]). `flats` represents a two-dimensional 32-bit integer
  array of the same size that is given by the output of
  `identifyflats`. `conncomps` represents a two-dimensional ptrdiff_t
  array of the same size and it is used to store the connected component
  labels.
 */
TOPOTOOLBOX_API
void gwdt_computecosts(float *costs, ptrdiff_t *conncomps, int32_t *flats,
                       float *original_dem, float *filled_dem,
                       ptrdiff_t dims[2]) {
  // flats, original_dem and filled_dem are passed as input
  //
  // conncomps and costs are outputs

  // Forward scan through the image. Compute local difference between
  // original and filled_dem. If the current pixel is flat, record its
  // connected component equivalence class and the maximum elevation
  // difference over its equivalence class in an implicit union-find
  // data structure

  // Offsets to the four neighbors that have already been visited
  // during the scan
  ptrdiff_t forward_j_offset[4] = {-1, -1, -1, 0};
  ptrdiff_t forward_i_offset[4] = {-1, 0, 1, -1};

  // All flats are on the interior of the image, so we don't need to
  // iterate over the borders
  for (ptrdiff_t j = 0; j < dims[1]; j++) {
    for (ptrdiff_t i = 0; i < dims[0]; i++) {
      ptrdiff_t current_pixel = j * dims[0] + i;
      if (!(flats[current_pixel] & 1)) {
        // Current pixel is not a flat
        costs[current_pixel] = 0.0;    // Set cost to zero
        conncomps[current_pixel] = 0;  // Set connected component label to zero
        continue;                      // Skip pixel
      }
      // Current pixel is a flat
      float current_depth =
          filled_dem[current_pixel] - original_dem[current_pixel];

      // Insert the current pixel into the union-find
      new_tree(conncomps, costs, current_pixel, current_depth);

      for (int32_t neighbor = 0; neighbor < 4; neighbor++) {
        ptrdiff_t neighbor_j = j + forward_j_offset[neighbor];
        ptrdiff_t neighbor_i = i + forward_i_offset[neighbor];
        ptrdiff_t neighbor_pixel = neighbor_j * dims[0] + neighbor_i;

        if (neighbor_j < 0 || neighbor_j >= dims[1] || neighbor_i < 0 ||
            neighbor_i >= dims[0]) {
          // This should be unreachable, because no border pixel is
          // also a flat pixel.
          continue;
        }

        if (flats[neighbor_pixel] & 1) {
          // Neighbor is a flat, unify it with the current pixel
          //
          // NOTE(wsk): It is possible to unify pixels somewhat more
          // efficiently by checking whether the neighbors are already
          // unified. This requires each neighboring pixel to be
          // checked slightly differently and is more confusing to
          // read, so it is left for future implementation if
          // performance improvements are needed.
          unify(conncomps, costs, neighbor_pixel, current_pixel);
        }
      }
    }
  }

  // Second scan: find the maximum difference for the connected
  // component of each flat
  // This time, we don't have to loop over the border pixels
  for (ptrdiff_t j = 1; j < dims[1] - 1; j++) {
    for (ptrdiff_t i = 1; i < dims[0] - 1; i++) {
      ptrdiff_t current_pixel = j * dims[0] + i;
      if (!(flats[current_pixel] & 1)) {
        continue;
      }
      labeldepth r = find_root(conncomps, costs, current_pixel);
      conncomps[current_pixel] = r.label;
      costs[current_pixel] = r.weight;
    }
  }

  // Third scan: compute the costs using the max depth for each
  // connected component
  for (ptrdiff_t j = 1; j < dims[1] - 1; j++) {
    for (ptrdiff_t i = 1; i < dims[0] - 1; i++) {
      ptrdiff_t current_pixel = j * dims[0] + i;
      if (!(flats[current_pixel] & 1)) {
        continue;
      }

      float current_depth =
          filled_dem[current_pixel] - original_dem[current_pixel];

      float tweight = 2.0f;
      float CarveMinVal = 0.1f;

      costs[current_pixel] =
          powf(costs[current_pixel] - current_depth, tweight) + CarveMinVal;
    }
  }
}

/*
  Gray-weighted distance transform

  Uses Dijkstra's algorithm with the priority queue implemented in
  priority_queue.c.

  Computes distances using the geodesic time algorithm of Soille
  (1994) and chamfer weights based on a quasi-Euclidean metric.
 */
TOPOTOOLBOX_API
void gwdt(float *dist, ptrdiff_t *prev, float *costs, int32_t *flats,
          ptrdiff_t *heap, ptrdiff_t *back, ptrdiff_t dims[2]) {
  // Initialize the priority queue
  PriorityQueue q = {0};
  q.back = back;
  q.heap = heap;
  q.priorities = dist;
  q.max_size = dims[0] * dims[1];
  for (ptrdiff_t j = 0; j < dims[1]; j++) {
    for (ptrdiff_t i = 0; i < dims[0]; i++) {
      ptrdiff_t idx = j * dims[0] + i;
      // prev points to self for any pixel that doesn't have a lower neighbor
      if (prev != NULL) {
        prev[idx] = idx;
      }

      dist[idx] = 0.0f;
      // Only put flat pixels into the queue
      if (flats[idx] & 1) {
        // Presill pixels are the sources ("seed locations"). They are
        // initialized with distance 1.0. Other flats are initialized at
        // infinity.
        float d = (flats[idx] & 4) ? 1.0f : INFINITY;
        pq_insert(&q, idx, d);
      }
    }
  }

  ptrdiff_t i_offset[8] = {-1, -1, -1, 0, 0, 1, 1, 1};
  ptrdiff_t j_offset[8] = {-1, 0, 1, -1, 1, -1, 0, 1};
  float chamfer[8] = {SQRT2f, 1.0, SQRT2f, 1.0, 1.0, SQRT2f, 1.0, SQRT2f};

  while (!pq_isempty(&q)) {
    ptrdiff_t trial = pq_deletemin(&q);
    float trial_distance = pq_get_priority(&q, trial);

    ptrdiff_t j = trial / dims[0];
    ptrdiff_t i = trial % dims[0];

    for (ptrdiff_t neighbor = 0; neighbor < 8; neighbor++) {
      ptrdiff_t neighbor_i = i + i_offset[neighbor];
      ptrdiff_t neighbor_j = j + j_offset[neighbor];
      ptrdiff_t neighbor_idx = neighbor_j * dims[0] + neighbor_i;

      // Skip pixels outside the boundary or non-flat pixels
      if (neighbor_i < 0 || neighbor_i >= dims[0] || neighbor_j < 0 ||
          neighbor_j >= dims[1] || !(flats[neighbor_idx] & 1)) {
        continue;
      }

      float proposal =
          trial_distance +
          chamfer[neighbor] * (costs[neighbor_idx] + costs[trial]) / 2;
      if (proposal < pq_get_priority(&q, neighbor_idx)) {
        pq_decrease_key(&q, neighbor_idx, proposal);
        if (prev != NULL) {
          prev[neighbor_idx] = trial;
        }
      }
    }
  }
}
