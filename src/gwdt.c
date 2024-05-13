#define TOPOTOOLBOX_BUILD

#include <math.h>
#include <stddef.h>
#include <stdint.h>

#include "topotoolbox.h"

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
  function, namely a ptrdiff_t array of size (nrows,ncols), which
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
  (nrows, ncols). `flats` represents a two-dimensional 32-bit integer
  array of the same size that is given by the output of
  `identifyflats`. `conncomps` represents a two-dimensional ptrdiff_t
  array of the same size and it is used to store the connected component
  labels.
 */
TOPOTOOLBOX_API
void gwdt_computecosts(float *costs, ptrdiff_t *conncomps, int32_t *flats,
                       float *original_dem, float *filled_dem, ptrdiff_t nrows,
                       ptrdiff_t ncols) {
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
  ptrdiff_t forward_col_offset[4] = {-1, -1, -1, 0};
  ptrdiff_t forward_row_offset[4] = {-1, 0, 1, -1};

  // All flats are on the interior of the image, so we don't need to
  // iterate over the borders
  for (ptrdiff_t col = 0; col < ncols; col++) {
    for (ptrdiff_t row = 0; row < nrows; row++) {
      ptrdiff_t current_pixel = col * nrows + row;
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
        ptrdiff_t neighbor_col = col + forward_col_offset[neighbor];
        ptrdiff_t neighbor_row = row + forward_row_offset[neighbor];
        ptrdiff_t neighbor_pixel = neighbor_col * nrows + neighbor_row;

        if (neighbor_col < 0 || neighbor_col >= ncols || row < 0 ||
            row >= ncols) {
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
  for (ptrdiff_t col = 1; col < ncols - 1; col++) {
    for (ptrdiff_t row = 1; row < nrows - 1; row++) {
      ptrdiff_t current_pixel = col * nrows + row;
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
  for (ptrdiff_t col = 1; col < ncols - 1; col++) {
    for (ptrdiff_t row = 1; row < nrows - 1; row++) {
      ptrdiff_t current_pixel = col * nrows + row;
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
