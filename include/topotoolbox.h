/**
   @file topotoolbox.h
   @version 3.0.0

   @brief Public API for libtopotoolbox
 */
#ifndef TOPOTOOLBOX_H
#define TOPOTOOLBOX_H

#ifdef TOPOTOOLBOX_STATICLIB
#define TOPOTOOLBOX_API
#elif defined(TOPOTOOLBOX_BUILD)
#if defined(_WIN32)
#define TOPOTOOLBOX_API __declspec(dllexport)
#else
#define TOPOTOOLBOX_API
#endif
#else
#if defined(_WIN32)
#define TOPOTOOLBOX_API __declspec(dllimport)
#else
#define TOPOTOOLBOX_API
#endif
#endif

#define TOPOTOOLBOX_VERSION_MAJOR 3
#define TOPOTOOLBOX_VERSION_MINOR 0
#define TOPOTOOLBOX_VERSION_PATCH 0

#ifndef TOPOTOOLBOX_OPENMP_VERSION
#define TOPOTOOLBOX_OPENMP_VERSION 0
#endif

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

/**
   @brief Test if topotoolbox is present.

   @details
   Used to ensure that topotoolbox is compiled and linked correctly.

   @return
   Always returns 1.
 */
TOPOTOOLBOX_API
int has_topotoolbox(void);

/**
   @brief Fills sinks in a digital elevation model

   @details
   Uses an algorithm based on grayscale morphological reconstruction.

   @param[out] output The filled DEM
   @parblock
   A pointer to a `float` array of size `dims[0]` x `dims[1]`
   @endparblock

   @param[in] dem The input DEM
   @parblock
   A pointer to a `float` array of size `dims[0]` x `dims[1]`
   @endparblock

   @param[in] dims The dimensions of the arrays
   @parblock
   A pointer to a `ptrdiff_t` array of size 2

   The fastest changing dimension should be provided first. For column-major
   arrays, `dims = {nrows,ncols}`. For row-major arrays, `dims = {ncols,nrows}`.
   @endparblock
 */
TOPOTOOLBOX_API
void fillsinks(float *output, float *dem, ptrdiff_t dims[2]);

/**
   @brief Fills sinks in a digital elevation model

   @details
   Uses an algorithm based on grayscale morphological
   reconstruction. Uses the hybrid algorithm of Vincent (1993) for
   higher performance than fillsinks(), but requires additional memory
   allocation for a FIFO queue.

   # References

   Vincent, Luc. (1993). Morphological grayscale reconstruction in
   image analysis: applications and efficient algorithms. IEEE
   Transactions on Image Processing, Vol. 2, No. 2.
   https://doi.org/10.1109/83.217222

   @param[out] output The filled DEM
   @parblock
   A pointer to a `float` array of size `dims[0]` x `dims[1]`
   @endparblock

   @param queue A pixel queue
   @parblock
   A pointer to a `ptrdiff_t` array of size `dims[0]` x `dims[1]`

   This array is used internally as the backing store for the necessary FIFO
   queue. It does not need to be initialized and can be freed once
   fillsinks_hybrid() returns.
   @endparblock

   @param[in] dem The input DEM
   @parblock
   A pointer to a `float` array of size `dims[0]` x `dims[1]`
   @endparblock

   @param[in] dims The dimensions of the arrays
   @parblock
   A pointer to a `ptrdiff_t` array of size 2

   The fastest changing dimension should be provided first. For column-major
   arrays, `dims = {nrows,ncols}`. For row-major arrays, `dims = {ncols,nrows}`.
   @endparblock
 */
TOPOTOOLBOX_API
void fillsinks_hybrid(float *output, ptrdiff_t *queue, float *dem,
                      ptrdiff_t dims[2]);

/**
   @brief Labels flat, sill and presill pixels in the provided DEM

   @details
   A flat pixel is one surrounded by pixels with the same or higher
   elevations. A sill pixel has the same elevation as a neighboring
   flat pixel but borders a pixel with a lower elevation. A presill
   pixel is a flat pixel that borders a sill pixel.

   The pixels are labeled with a bit field:

   - Bit 0: Set if pixel is a flat
   - Bit 1: Set if pixel is a sill
   - Bit 2: Set if pixel is a presill

   Since all presill pixels are also flats, bits 0 and 2 are both set
   for presill pixels. In other words presill pixels have a value of 5
   in the output array. This allows one to test the identity of pixels
   with bitwise operations:

   ```
   if (output[pixel] & 1) {
     // Pixel is a flat
   }
   if (output[pixel] & 2) {
     // Pixel is a sill
   }
   if (output[pixel] & 4) {
     // Pixel is a presill
   }
   ```

   @param[out] output The pixel labels
   @parblock
   A pointer to an `int32_t` array of size `dims[0]` x `dims[1]`
   @endparblock

   @param[in] dem The input DEM
   @parblock
   A pointer to a `float` array of size `dims[0]` x `dims[1]`
   @endparblock

   @param[in] dims The dimensions of the arrays
   @parblock
   A pointer to a `ptrdiff_t` array of size 2

   The fastest changing dimension should be provided first. For column-major
   arrays, `dims = {nrows,ncols}`. For row-major arrays, `dims = {ncols,nrows}`.
   @endparblock
 */
TOPOTOOLBOX_API
ptrdiff_t identifyflats(int32_t *output, float *dem, ptrdiff_t dims[2]);

/**
   @brief Compute costs for the gray-weighted distance transform

   @details
   The costs used to route flow over flat regions with the
   gray-weighted distance transform are based on the difference between
   original and filled DEMs at each pixel. This difference is
   subtracted from the maximum difference over the flat region to which
   the pixel belongs. It is then squared and a small constant
   (currently set to 0.1) is added.

   This function identifies the connected components of flat pixels,
   provided in the `flats` array, computes the maximum difference
   between the `original_dem` and the `filled_dem` over each connected
   component and then computes the cost for each flat pixel. The
   `costs` are returned, as are the connected components labels
   (`conncomps`). These labels are the linear index of the pixel in the
   connected component with the maximum difference between the filled
   and the output DEMs.

   @param[out] costs The gray-weighted distance transform costs
   @parblock
   A pointer to a `float` array of size `dims[0]` x `dims[1]`
   @endparblock

   @param[out] conncomps Labeled connected components for each flat pixel
   @parblock
   A pointer to a `ptrdiff_t` array of size `dims[0]` x `dims[1]`
   @endparblock

   @param[in] flats Array identifying the flat pixels
   @parblock
   A pointer to an `int32_t` array of size `dims[0]` x `dims[1]`

   The flat pixels must be identified as they are by identifyflats() such that
   `flats[pixels] & 1` is nonzero for any flat pixel.
   @endparblock

   @param[in] original_dem The DEM prior to sink filling
   @parblock
   A pointer to a `float` array of size `dims[0]` x `dims[1]`
   @endparblock

   @param[in] filled_dem The DEM after sink filling
   @parblock
   A pointer to a `float` array of size `dims[0]` x `dims[1]`
   @endparblock

   @param[in] dims The dimensions of the arrays
   @parblock
   A pointer to a `ptrdiff_t` array of size 2

   The fastest changing dimension should be provided first. For column-major
   arrays, `dims = {nrows,ncols}`. For row-major arrays, `dims = {ncols,nrows}`.
   @endparblock
 */
TOPOTOOLBOX_API
void gwdt_computecosts(float *costs, ptrdiff_t *conncomps, int32_t *flats,
                       float *original_dem, float *filled_dem,
                       ptrdiff_t dims[2]);

/**
   @brief Compute the gray-weighted distance transform

   @details
   This gray-weighted distance transform uses Dijkstra's algorithm to
   compute the geodesic time of Soille (1994) using the provided
   costs raster. The `flats` array, which could be generated by
   identifyflats(), controls which pixels are considered in the
   distance transform. Any pixel such that `(flats[pixel] & 1) != 0`
   is considered by the algorithm, and the algorithm starts at source
   pixels where `(flats[pixel] & 4) != 0`. All other pixels are
   considered barriers through which paths cannot pass.

   Chamfer weights are multiplied by the cost for each edge based on a
   Euclidean metric: corner pixels are multiplied by sqrt(2).

   # References

   Soille, Pierre (1994). Generalized geodesy via geodesic
   time. Pattern Recognition Letters 15, 1235-1240.

   @param[out] dist The computed gray-weighted distance transform
   @parblock
   A pointer to a `float` array of size `dims[0]` x `dims[1]`
   @endparblock

   @param[out] prev Backlinks along the geodesic path
   @parblock
   A pointer to a `ptrdiff_t` array of size `dims[0]` x `dims[1]`

   If backlinks are not required, a null pointer can be passed here: it is
   checked for NULL before being accessed.
   @endparblock

   @param[in] costs The input costs computed by gwdt_computecosts()
   @parblock
   A pointer to a `float` array of size `dims[0]` x `dims[1]`
   @endparblock

   @param[in] flats Array identifying the flat pixels
   @parblock
   A pointer to an `int32_t` array of size `dims[0]` x `dims[1]`

   The flat pixels must be identified as they are by identifyflats().
   @endparblock

   @param heap Storage for the priority queue
   @parblock
   A pointer to a `ptrdiff_t` array of size `dims[0]` x `dims[1]`
   @endparblock

   @param back Storage for the priority queue
   @parblock
   A pointer to a `ptrdiff_t` array of indices `dims[0]` x `dims[1]`
   @endparblock

   @param[in] dims The dimensions of the arrays
   @parblock
   A pointer to a `ptrdiff_t` array of size 2

   The fastest changing dimension should be provided first. For column-major
   arrays, `dims = {nrows,ncols}`. For row-major arrays, `dims = {ncols,nrows}`.
   @endparblock
 */
TOPOTOOLBOX_API
void gwdt(float *dist, ptrdiff_t *prev, float *costs, int32_t *flats,
          ptrdiff_t *heap, ptrdiff_t *back, ptrdiff_t dims[2]);

/**
   @brief Compute excess topography with 2D varying threshold slopes
   using the fast sweeping method

   @details
   The excess topography (Blöthe et al. 2015) is computed by solving an
   eikonal equation (Anand et al. 2023) constrained to lie below the
   original DEM. Where the slope of the DEM is greater than the
   threshold slope, the eikonal solver limits the output topography to
   that slope, but where the slope of the DEM is lower that the
   threshold slope, the output follows the DEM.

   The eikonal equation is solved using the fast sweeping method (Zhao
   2004), which iterates over the DEM in alternating directions and
   updates the topography according to an upwind discretization of the
   gradient. To constrain the solution by the original DEM, the
   output topography is initiated with the DEM and only updates lower than the
   DEM are accepted.

   The fast sweeping method is simpler than the fast marching method
   (excesstopography_fmm2d()), requires less memory, and can be faster,
   particularly when the threshold slopes are constant or change
   infrequently across the domain.

   # References

   Anand, Shashank Kumar, Matteo B. Bertagni, Arvind Singh and Amilcare
   Porporato (2023). Eikonal equation reproduces natural landscapes with
   threshold hillslopes. Geophysical Research Letters, 50, 21.

   Blöthe, Jan Henrik, Oliver Korup and Wolfgang Schwanghart
   (2015). Large landslides lie low: Excess topography in the
   Himalaya-Karakoram ranges. Geology, 43, 6, 523-526.

   Zhao, Hongkai (2004). A fast sweeping method for eikonal
   equations. Mathematics of Computation, 74, 250, 603-627.

   @param[out] excess The solution of the constrained eikonal equation
   @parblock
   A pointer to a `float` array of size `dims[0]` x `dims[1]`

   To compute the excess topography, subtract this array elementwise from the
   DEM.
   @endparblock

   @param[in] dem The input digital elevation model.
   @parblock
   A pointer to a `float` array of size `dims[0]` x `dims[1]`
   @endparblock

   @param[in] threshold_slopes The threshold slopes at each grid cell.
   @parblock
   A pointer to a `float` array of size `dims[0]` x `dims[1]`
   @endparblock

   @param[in] cellsize The spacing between grid cells
   @parblock
   A `float`

   The spacing is assumed to be constant and identical in the x- and y-
   directions.
   @endparblock

   @param[in] dims The dimensions of the arrays
   @parblock
   A pointer to a `ptrdiff_t` array of size 2

   The fastest changing dimension should be provided first. For column-major
   arrays, `dims = {nrows,ncols}`. For row-major arrays, `dims = {ncols,nrows}`.
   @endparblock
 */
TOPOTOOLBOX_API
void excesstopography_fsm2d(float *excess, float *dem, float *threshold_slopes,
                            float cellsize, ptrdiff_t dims[2]);

/**
   @brief Compute excess topography with 2D varying threshold slopes
   using the fast marching method

   @details
   The excess topography (Blöthe et al. 2015) is computed by solving an
   eikonal equation (Anand et al. 2023) constrained to lie below the
   original DEM. Where the slope of the DEM is greater than the
   threshold slope, the eikonal solver limits the output topography to
   that slope, but where the slope of the DEM is lower that the
   threshold slope, the output follows the DEM.

   The eikonal equation is solved using the fast marching method
   (Sethian 1996), which uses a priority queue to propagate slopes from
   the lowest elevation pixels to the highest according to an upwind
   discretization of the gradient. To constrain the solution by the
   original DEM, the output topography is initiated with the DEM and
   only updates lower than the DEM are accepted.

   The fast marching method is more complicated than the fast sweeping
   method (excesstopography_fmm2d()) and requires pre-allocated memory
   for the priority queue. It is faster than the fast sweeping method
   when the threshold slopes change frequently.

   # References

   Anand, Shashank Kumar, Matteo B. Bertagni, Arvind Singh and Amilcare
   Porporato (2023). Eikonal equation reproduces natural landscapes with
   threshold hillslopes. Geophysical Research Letters, 50, 21.

   Blöthe, Jan Henrik, Oliver Korup and Wolfgang Schwanghart
   (2015). Large landslides lie low: Excess topography in the
   Himalaya-Karakoram ranges. Geology, 43, 6, 523-526.

   Sethian, James (1996). A fast marching level set method for
   monotonically advancing fronts. Proceedings of the National Academy
   of Sciences, 93, 4, 1591-1595.

   @param[out] excess The solution of the constrained eikonal equation
   @parblock
   A pointer to a `float` array of size `dims[0]` x `dims[1]`

   To compute the excess topography, subtract this array elementwise from the
   DEM.
   @endparblock

   @param heap Storage for the priority queue
   @parblock
   A pointer to a `ptrdiff_t` array of size `dims[0]` x `dims[1]`
   @endparblock

   @param back Storage for the priority queue
   @parblock
   A pointer to a `ptrdiff_t` array of indices `dims[0]` x `dims[1]`
   @endparblock

   @param[in] dem The input digital elevation model.
   @parblock
   A pointer to a `float` array of size `dims[0]` x `dims[1]`
   @endparblock

   @param[in] threshold_slopes The threshold slopes at each grid cell.
   @parblock
   A pointer to a `float` array of size `dims[0]` x `dims[1]`
   @endparblock

   @param[in] cellsize The spacing between grid cells
   @parblock
   A `float`

   The spacing is assumed to be constant and identical in the x- and y-
   directions.
   @endparblock

   @param[in] dims The dimensions of the arrays
   @parblock
   A pointer to a `ptrdiff_t` array of size 2

   The fastest changing dimension should be provided first. For column-major
   arrays, `dims = {nrows,ncols}`. For row-major arrays, `dims = {ncols,nrows}`.
   @endparblock

 */
TOPOTOOLBOX_API
void excesstopography_fmm2d(float *excess, ptrdiff_t *heap, ptrdiff_t *back,
                            float *dem, float *threshold_slopes, float cellsize,
                            ptrdiff_t dims[2]);

/**
   @brief Compute excess topography with three-dimensionally variable
   lithology using the fast marching method

   @details
   The excess topography is computed by solving an eikonal equation with
   the fast marching method (see excesstopography_fmm2d() for more
   information). The threshold slope at a grid cell is computed from
   that cell's position within the three-dimensional lithology.

   The lithology consists of a set of discrete layers, each of which has
   its own threshold slope, which is specified by the caller in the
   `threshold_slopes` array from the bottom layer to the top layer. The
   layer geometry is provided by the caller in the `lithstack` array,
   which holds the elevation of the top surface of each layer at each
   grid cell.

   The algorithm proceeds similarly to the regular fast marching method,
   using a priority queue to update grid cells from bottom to
   top. Whenever a cell is updated, the eikonal solver is used to
   propose a new elevation using the threshold slopes from each layer
   from bottom to top. The first elevation that is proposed that lies
   below the top surface of the layer whose slope is used in the
   proposal is accepted as the provisional height for that grid cell.

   @param[out] excess The solution of the constrained eikonal equation
   @parblock
   A pointer to a `float` array of size `dims[0]` x `dims[1]`

   To compute the excess topography, subtract this array elementwise from the
   DEM.
   @endparblock

   @param heap Storage for the priority queue
   @parblock
   A pointer to a `ptrdiff_t` array of size `dims[0]` x `dims[1]`
   @endparblock

   @param back Storage for the priority queue
   @parblock
   A pointer to a `ptrdiff_t` array of indices `dims[0]` x `dims[1]`
   @endparblock

   @param[in] dem The input digital elevation model
   @parblock
   A pointer to a `float` array of size `dims[0]` x `dims[1]`
   @endparblock

   @param[in] lithstack The input lithology.
   @parblock
   A pointer to a `float` array of size `nlayers` x `dims[0]` x `dims[1]`

   The value of `lithstack[layer,row,col]` is the elevation of the top
   surface of the given layer.  Note that the first dimension is the
   layer, so that the layers of each cell are stored contiguously.
   @endparblock

   @param[in] threshold_slopes The threshold slopes for each layer
   @parblock
   A pointer to a `float` array of size `nlayers`
   @endparblock

   @param[in] cellsize The spacing between grid cells
   @parblock
   A `float`

   The spacing is assumed to be constant and identical in the x- and y-
   directions.
   @endparblock

   @param[in] dims The horizontal dimensions of the arrays
   @parblock
   A pointer to a `ptrdiff_t` array of size 2

   The fastest changing dimension should be provided first. For column-major
   arrays, `dims = {nrows,ncols}`. For row-major arrays, `dims = {ncols,nrows}`.
   @endparblock
   @param[in] nlayers The number of layers in lithstack and threshold_slopes

 */
TOPOTOOLBOX_API
void excesstopography_fmm3d(float *excess, ptrdiff_t *heap, ptrdiff_t *back,
                            float *dem, float *lithstack,
                            float *threshold_slopes, float cellsize,
                            ptrdiff_t dims[2], ptrdiff_t nlayers);

/**
   @brief Route flow over the DEM using the D8 method

   @details
   The flow routing is solved using a D8 steepest descent
   algorithm. Flat regions, which are identified in the `flats` array
   using the indexing scheme of identifyflats(), are routed by carving
   using the auxiliary topography provided in the `dist` array, which
   can be generated by gwdt().

   @param[out] source The source pixel for each edge
   @parblock
   A pointer to a `ptrdiff_t` array of size `dims[0]` x `dims[1]`

   The source pixels are sorted topologically so that pixel 0
   indicates the first source to the flow network.
   @endparblock

   @param[out] direction The flow directions as a bit field.
   @parblock
   A pointer to a `uint8_t` array of size `dims[0]` x `dims[1]`

   The 8 bits (0-7) identify the downstream neighbor of pixel `(i,j)`
   as follows:

   ```
   --- | j-1 | j | j+1|
   ----+-----+---+----|
   i-1 | 5   | 6 | 7  |
   i   | 4   |   | 0  |
   i+1 | 3   | 2 | 1  |
   ```

   For example, a pixel with its downstream neighbor at `(i+1,j-1)`
   has a value in the direction array of `0b00001000 = 8`. A value of
   0 indicates that the pixel has no downstream neighbors and is
   either a sink or an outlet. A value of 255 (all bits set) is used
   internally as a sentinel value and its presence in output data is a
   sign of an error.
   @endparblock

   @param[in] dem The input digital elevation model
   @parblock
   A pointer to a `float` array of size `dims[0]` x `dims[1]`
   @endparblock

   @param[in] dist The auxiliary topography for routing over floats
   @parblock
   A pointer to a `float` array of size `dims[0]` x `dims[1]`

   This will typically be generated by gwdt() as the output `dist`.
   @endparblock

   @param[in] flats Array identifying the flat pixels
   @parblock
   A pointer to an `int32_t` array of size `dims[0]` x `dims[1]`

   The flat pixels must be identified as they are by identifyflats().
   @endparblock

   @param[in] dims The dimensions of the arrays
   @parblock
   A pointer to a `ptrdiff_t` array of size 2

   The fastest changing dimension should be provided first. For column-major
   arrays, `dims = {nrows,ncols}`. For row-major arrays, `dims = {ncols,nrows}`.
   @endparblock
 */
TOPOTOOLBOX_API
void flow_routing_d8_carve(ptrdiff_t *source, uint8_t *direction, float *dem,
                           float *dist, int32_t *flats, ptrdiff_t dims[2]);

/**
   @brief Compute downstream pixel indices from flow directions

   The `source` and `direction` outputs from flow_routing_d8_carve()
   implicitly define the downstream targets of each edge in the flow
   network. This function computes the linear indices of those
   downstream targets and stores them in the `target` array.

   @param[out] target The target pixel for each edge
   @parblock
   A pointer to a `ptrdiff_t` array of size `dims[0]` x `dims[1]`
   @endparblock

   @param[in] source The source pixel for each edge
   @parblock
   A pointer to a `ptrdiff_t` array of size `dims[0]` x `dims[1]`
   @endparblock

   @param[in] direction The flow directions as a bit field
   @parblock
   A pointer to a `uint8_t` array of size `dims[0]` x `dims[1]`

   The flow directions should be encoded as they are in flow_routing_d8_carve().
   @endparblock

   @param[in] dims The dimensions of the arrays
   @parblock
   A pointer to a `ptrdiff_t` array of size 2

   The fastest changing dimension should be provided first. For column-major
   arrays, `dims = {nrows,ncols}`. For row-major arrays, `dims = {ncols,nrows}`.
   @endparblock
 */
TOPOTOOLBOX_API
void flow_routing_targets(ptrdiff_t *target, ptrdiff_t *source,
                          uint8_t *direction, ptrdiff_t dims[2]);

/**
   @brief Compute flow accumulation

   Accumulates flow by summing contributing areas along flow paths. Uses the
   `source` and `direction` outputs of flow_routing_d8_carve().

   @param[out] acc The computed flow accumulation
   @parblock
   A pointer to a `float` array of size `dims[0]` x `dims[1]`
   @endparblock

   @param[in] source The source pixel for each edge
   @parblock
   A pointer to a `ptrdiff_t` array of size `dims[0]` x `dims[1]`

   The source pixels must be in a topological order.
   @endparblock

   @param[in] direction The flow directions as a bit field
   @parblock
   A pointer to a `uint8_t` array of size `dims[0]` x `dims[1]`

   The flow directions should be encoded as they are in flow_routing_d8_carve().
   @endparblock

   @param[in] weights Initial water depths
   @parblock
   A pointer to a `float` array of size `dims[0]` x `dims[1]`

   The initial weights can be used to represent spatially variable
   precipitation.

   If a null pointer is passed, a default weight of 1.0 for every
   pixel is used. In this case the resulting flow accumulation is the
   upstream area in number of pixels.
   @endparblock

   @param[in] dims The dimensions of the arrays
   @parblock
   A pointer to a `ptrdiff_t` array of size 2

   The fastest changing dimension should be provided first. For column-major
   arrays, `dims = {nrows,ncols}`. For row-major arrays, `dims = {ncols,nrows}`.
   @endparblock
 */
TOPOTOOLBOX_API
void flow_accumulation(float *acc, ptrdiff_t *source, uint8_t *direction,
                       float *weights, ptrdiff_t dims[2]);

#include "graphflood/define_types.h"

/**
   @brief Computes a single flow graph:
   Receivers/Donors using the steepest descent method and topological
   ordering following a modified Braun and Willett (2013)

   @param[in]  topo: the topographic surface
   @param[out] Sreceivers: array of steepest receiver vectorised index
   @param[out] distToReceivers: array of distance to steepest receiver
   vectorised index
   @param[out] Sdonors: array of donors to steepest receiver vectorised
   index (index * (8 or 4) + 0:NSdonors[index] to get them)
   @param[out] NSdonors: array of number of steepest donors (nodes having
   this one as steepest receivers)
   @param[out] Stack: topologically ordered list of nodes, from the
   baselevel to the sources
   @param[in]  BCs: codes for boundary conditions and no data management,
   see gf_utils.h or examples for the meaning
   @param[in]  dims: [rows,columns] if row major and [columns, rows] if
   column major
   @param[in]  dx: spatial step
   @param[in]  D8: true for topology including cardinals + diagonals, false
   for cardinals only
*/
TOPOTOOLBOX_API
void compute_sfgraph(GF_FLOAT *topo, GF_UINT *Sreceivers,
                     GF_FLOAT *distToReceivers, GF_UINT *Sdonors,
                     uint8_t *NSdonors, GF_UINT *Stack, uint8_t *BCs,
                     GF_UINT *dim, GF_FLOAT dx, bool D8);

/**
   @brief Compute the graphflood single flow graph and fills local minima using
Priority Floods - Barnes 2014 (see compute_sfgraph for details)
*/
TOPOTOOLBOX_API
void compute_sfgraph_priority_flood(GF_FLOAT *topo, GF_UINT *Sreceivers,
                                    GF_FLOAT *distToReceivers, GF_UINT *Sdonors,
                                    uint8_t *NSdonors, GF_UINT *Stack,
                                    uint8_t *BCs, GF_UINT *dim, GF_FLOAT dx,
                                    bool D8);

/**
   @brief Fills the depressions in place in the topography using Priority
   Floods Barnes (2014, modified to impose a minimal slope)

   @param[inout]  topo: array of surface elevation
   @param[in]     BCs: codes for boundary conditions and no data
   management, see gf_utils.h or examples for the meaning
   @param[in]     dims: [rows,columns] if row major and [columns, rows] if
   column major
   @param[in]     D8: true for topology including cardinals + diagonals,
   false for cardinals only
*/
TOPOTOOLBOX_API
void compute_priority_flood(float *topo, uint8_t *BCs, GF_UINT *dim, bool D8);

/**
   @brief Fills the depressions in place in the topography using Priority
   Floods Barnes (2014, modified to impose a minimal slope) This variant
   computes the topological order on the go (slightly slower as it uses a
   priority queue for all the nodes including in depressions)

   @param[inout]  topo: array of surface elevation
   @param[in]     BCs: codes for boundary conditions and no data
   management, see gf_utils.h or examples for the meaning
   @param[in]     dims: [rows,columns] if row major and [columns, rows] if
   column major
   @param[in]     D8: true for topology including cardinals + diagonals,
   false for cardinals only
*/
TOPOTOOLBOX_API
void compute_priority_flood_plus_topological_ordering(float *topo,
                                                      GF_UINT *stack,
                                                      uint8_t *BCs,
                                                      GF_UINT *dim, bool D8);

/**
   @brief Accumulate single flow drainage area downstream from a calculated
   graphflood single flow graph
   @param[out] output: the field of drainage area
   @param[in]  Sreceivers: array of steepest receiver vectorised index
   @param[in]  Stack: topologically ordered list of nodes, from the
   baselevel to the sources
   @param[in]  dims: [rows,columns] if row major and [columns, rows] if
   column major
   @param[in]  dx: spatial step
*/
TOPOTOOLBOX_API
void compute_drainage_area_single_flow(GF_FLOAT *output, GF_UINT *Sreceivers,
                                       GF_UINT *Stack, GF_UINT *dim,
                                       GF_FLOAT dx);

/**
   @brief Accumulate single flow drainage area downstream from a calculated
   graphflood single flow graph weighted by an arbitrary input (e.g.
   Precipitation rates to get effective discharge)
   @param[out] output: the field of drainage area
   @param[in]  weights: node-wise weights
   @param[in]  Sreceivers: array of steepest receiver vectorised index
   @param[in]  Stack: topologically ordered list of nodes, from the
   baselevel to the sources
   @param[in]  dims: [rows,columns] if row major and [columns, rows] if
   column major
   @param[in]  dx: spatial step
*/
TOPOTOOLBOX_API
void gradient8(float *output, float *dem, float cellsize, int use_mp,
               ptrdiff_t dims[2]);
/**
   @brief Compute  the gradient for each cell in the provided DEM array.
   The gradient is calculated as the maximum slope between the cell and its
   8 neighboring cells. The result can be output in different units based on the
   `unit` parameter, and the computation can be parallelized using OpenMP.

   @param[out] output: Array to store the computed gradient values for each
   cell. It should have the same dimensions as the DEM.
   @param[in]  dem: Input digital elevation model as a 2D array flattened into a
   1D array. This array represents the elevation values of each cell.
   @param[in]  cellsize: The spatial resolution of the DEM (i.e., the size of
   each cell).
   @param[in]  use_mp: If set to 1, enables parallel processing using OpenMP.
                       If set to 0, the function runs on a single thread.
   @param[in]  dims: An array specifying the dimensions of the DEM.
                     It should contain two values: [rows, columns].
*/
TOPOTOOLBOX_API
void compute_weighted_drainage_area_single_flow(GF_FLOAT *output,
                                                GF_FLOAT *weights,
                                                GF_UINT *Sreceivers,
                                                GF_UINT *Stack, GF_UINT *dim,
                                                GF_FLOAT dx);

/**
   @brief Run N iteration of graphflood as described in Gailleton et al.,
   2024. From an input field of topography, optional original flow depth,
   mannings friction coefficient and precipitation rates, calculates the field
   of flow depth following a steady flow assumption.

   @param[in]     Z: surface topography
   @param[inout]  hw: field of flow depth
   @param[in]     BCs: codes for boundary conditions and no data
   management, see gf_utils.h or examples for the meaning
   @param[in]     Precipitations: Precipitation rates
   @param[in]     manning: friction coefficient
   @param[in]     dims: [rows,columns] if row major and [columns, rows] if
   column major
   @param[in]     dt: time step
   @param[in]     dx: spatial step
   @param[in]     SFD: single flow direction if True, multiple flow if
   false
*/
TOPOTOOLBOX_API
void graphflood_full(GF_FLOAT *Z, GF_FLOAT *hw, uint8_t *BCs,
                     GF_FLOAT *Precipitations, GF_FLOAT *manning, GF_UINT *dim,
                     GF_FLOAT dt, GF_FLOAT dx, bool SFD, bool D8,
                     GF_UINT N_iterations);

#endif  // TOPOTOOLBOX_H
