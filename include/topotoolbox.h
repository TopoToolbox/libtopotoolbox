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
   @brief Minimum value filter

   @remark
   All arrays passed onto the function must be non-overlapping regions of
   memory.

   @note
   All structuring elements must be of the same size.

   @note
   Pixels at the border of the input image are not set to 0, NAN or +/- INFINITY
   except when those are their input values. If the structuring element extends
   over the input array boundaries, parts of it are disregarded and only a
   subset of the structuring element is used for the computation.

   @param[out] output The filtered DEM
   @parblock
   A pointer to a `float` array of size `io_dims[0]` x `io_dims[1]`
   @endparblock

   @param[in] dem The input DEM
   @parblock
   A pointer to a `float` array of size `io_dims[0]` x `io_dims[1]`
   @endparblock

   @param[in] structuring_element The structuring elements
   @parblock
   A pointer to a `uint8_t` array of size `se_dims[0]` x `se_dims[1]` x
   `se_dims[2]`
   @endparblock

   @param[in] io_dims The dimensions of the input/output DEM arrays
   @parblock
   A pointer to a `ptrdiff_t` array of size 2

   The fastest changing dimension should be provided first. For column-major
   arrays, `io_dims = {nrows,ncols}`. For row-major arrays, `io_dims =
   {ncols,nrows}`.
   @endparblock

   @param[in] se_dims The dimensions of the SE arrays
   @parblock
   A pointer to a `ptrdiff_t` array of size 3

   The fastest changing dimension should be provided first. For column-major
   arrays, `se_dims = {nrows,ncols,SE entity}`.
   For row-major arrays, `se_dims = {ncols,nrows,SE entity}`.
   @endparblock
 */
TOPOTOOLBOX_API
void min_filter(float *output, float *dem, uint8_t *structuring_element,
                ptrdiff_t io_dims[2], ptrdiff_t se_dims[3]);

/**
   @brief Maximum value filter

   @remark
   All arrays passed onto the function must be non-overlapping regions of
   memory.

   @note
   All structuring elements must be of the same size.

   @note
   Pixels at the border of the input image are not set to 0, NAN or +/- INFINITY
   except when those are their input values. If the structuring element extends
   over the input array boundaries, parts of it are disregarded and only a
   subset of the structuring element is used for the computation.

   @param[out] output The filtered DEM
   @parblock
   A pointer to a `float` array of size `io_dims[0]` x `io_dims[1]`
   @endparblock

   @param[in] dem The input DEM
   @parblock
   A pointer to a `float` array of size `io_dims[0]` x `io_dims[1]`
   @endparblock

   @param[in] structuring_element The structuring elements
   @parblock
   A pointer to a `uint8_t` array of size `se_dims[0]` x `se_dims[1]` x
   `se_dims[2]`
   @endparblock

   @param[in] io_dims The dimensions of the input/output DEM arrays
   @parblock
   A pointer to a `ptrdiff_t` array of size 2

   The fastest changing dimension should be provided first. For column-major
   arrays, `io_dims = {nrows,ncols}`. For row-major arrays, `io_dims =
   {ncols,nrows}`.
   @endparblock

   @param[in] se_dims The dimensions of the SE arrays
   @parblock
   A pointer to a `ptrdiff_t` array of size 3

   The fastest changing dimension should be provided first. For column-major
   arrays, `se_dims = {nrows,ncols,SE entity}`.
   For row-major arrays, `se_dims = {ncols,nrows,SE entity}`.
   @endparblock
 */
TOPOTOOLBOX_API
void max_filter(float *output, float *dem, uint8_t *structuring_element,
                ptrdiff_t io_dims[2], ptrdiff_t se_dims[3]);

/**
   @brief Minimum value filter, optimized for square and full structuring
   elements

   @remark
   All arrays passed onto the function must be non-overlapping regions of
   memory.

   @note
   Pixels at the border of the input image are not set to 0, NAN or +/- INFINITY
   except when those are their input values. If the structuring element extends
   over the input array boundaries, parts of it are disregarded and only a
   subset of the structuring element is used for the computation.

   @param[out] output The filtered DEM
   @parblock
   A pointer to a `float` array of size `io_dims[0]` x `io_dims[1]`
   @endparblock

   @param[in] dem The input DEM
   @parblock
   A pointer to a `float` array of size `io_dims[0]` x `io_dims[1]`
   @endparblock

   @param[in] tmp Temporary value storage
   @parblock
   This function trades computation time for memory space and requires
   an additional `float` array of size `io_dims[0]` x `io_dims[1]`
   to accomodate intermediate values.
   @endparblock

   @param[in] width The structuring elements
   @parblock
   The width and height of a squre structuring element which is filled with 1s
   @endparblock

   @param[in] io_dims The dimensions of the input/output DEM arrays
   @parblock
   A pointer to a `ptrdiff_t` array of size 2

   The fastest changing dimension should be provided first. For column-major
   arrays, `io_dims = {nrows,ncols}`. For row-major arrays, `io_dims =
   {ncols,nrows}`.
   @endparblock
 */
TOPOTOOLBOX_API
void min_filter_square(float *output, float *dem, float *tmp, uint8_t width,
                       ptrdiff_t io_dims[2]);

/**
   @brief Maximum value filter, optimized for square and full structuring
   elements

   @remark
   All arrays passed onto the function must be non-overlapping regions of
   memory.

   @note
   Pixels at the border of the input image are not set to 0, NAN or +/- INFINITY
   except when those are their input values. If the structuring element extends
   over the input array boundaries, parts of it are disregarded and only a
   subset of the structuring element is used for the computation.

   @param[out] output The filtered DEM
   @parblock
   A pointer to a `float` array of size `io_dims[0]` x `io_dims[1]`
   @endparblock

   @param[in] dem The input DEM
   @parblock
   A pointer to a `float` array of size `io_dims[0]` x `io_dims[1]`
   @endparblock

   @param[in] tmp Temporary value storage
   @parblock
   This function trades computation time for memory space and requires
   an additional `float` array of size `io_dims[0]` x `io_dims[1]`
   to accomodate intermediate values.
   @endparblock

   @param[in] width The structuring elements
   @parblock
   The width and height of a squre structuring element which is filled with 1s
   @endparblock

   @param[in] io_dims The dimensions of the input/output DEM arrays
   @parblock
   A pointer to a `ptrdiff_t` array of size 2

   The fastest changing dimension should be provided first. For column-major
   arrays, `io_dims = {nrows,ncols}`. For row-major arrays, `io_dims =
   {ncols,nrows}`.
   @endparblock
 */
TOPOTOOLBOX_API
void max_filter_square(float *output, float *dem, float *tmp, uint8_t width,
                       ptrdiff_t io_dims[2]);

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

   @param[in] bc Array used to set boundary conditions
   @parblock
   A pointer to a `uint8_t` array of size `dims[0]` x `dims[1]`

   `bc` is used to control which pixels get filled. Pixels that are
   set equal to 1 are fixed to their value in the input DEM while
   pixels equal to 0 are filled. For the standard fillsinks operation,
   bc equals 1 on the boundaries of the DEM and 0 on the interior. Set
   bc equal to 1 for NaN pixels to ensure that they are treated as
   sinks.
   @endparblock

   @param[in] dims The dimensions of the arrays
   @parblock
   A pointer to a `ptrdiff_t` array of size 2

   The fastest changing dimension should be provided first. For column-major
   arrays, `dims = {nrows,ncols}`. For row-major arrays, `dims = {ncols,nrows}`.
   @endparblock
 */
TOPOTOOLBOX_API
void fillsinks(float *output, float *dem, uint8_t *bc, ptrdiff_t dims[2]);

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

   @param[in] bc Array used to set boundary conditions
   @parblock
   A pointer to a `uint8_t` array of size `dims[0]` x `dims[1]`

   `bc` is used to control which pixels get filled. Pixels that are
   set equal to 1 are fixed to their value in the input DEM while
   pixels equal to 0 are filled. For the standard fillsinks operation,
   bc equals 1 on the boundaries of the DEM and 0 on the interior. Set
   bc equal to 1 for NaN pixels to ensure that they are treated as
   sinks.
   @endparblock

   @param[in] dims The dimensions of the arrays
   @parblock
   A pointer to a `ptrdiff_t` array of size 2

   The fastest changing dimension should be provided first. For column-major
   arrays, `dims = {nrows,ncols}`. For row-major arrays, `dims = {ncols,nrows}`.
   @endparblock
 */
TOPOTOOLBOX_API
void fillsinks_hybrid(float *output, ptrdiff_t *queue, float *dem, uint8_t *bc,
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

   @param[out] node Nodes of the flow graph in topological order
   @parblock
   A pointer to a `ptrdiff_t` array of size `dims[0]` x `dims[1]`

   The pixels of the dem are sorted topologically that every pixel
   that drains to pixel v comes before v in the array.
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

   @param[in] order The memory order of the underlying array
   @parblock
   0 for column-major, 1 for row-major

   @endparblock
 */
TOPOTOOLBOX_API
void flow_routing_d8_carve(ptrdiff_t *node, uint8_t *direction, float *dem,
                           float *dist, int32_t *flats, ptrdiff_t dims[2],
                           unsigned int order);

/**
   @brief Compute downstream pixel indices from flow directions

   The `node` and `direction` outputs from flow_routing_d8_carve()
   form a topologically sorted adjacency list representation of the
   flow graph. This function constructs an edge list stores it in the
   `source` and `target` arrays.

   @param[out] source The source pixel for each edge
   @parblock
   A pointer to a `ptrdiff_t` array of size `dims[0]` x `dims[1]`
   @endparblock

   @param[out] target The target pixel for each edge
   @parblock
   A pointer to a `ptrdiff_t` array of size `dims[0]` x `dims[1]`
   @endparblock

   @param[in] node Nodes of the flow graph in topological order
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

   @param[in] order The memory order of the underlying array
   @parblock
   0 for column-major, 1 for row-major

   @endparblock

   @return The number of valid edges contained in the `source` and `target`
   arrays
 */
TOPOTOOLBOX_API
ptrdiff_t flow_routing_d8_edgelist(ptrdiff_t *source, ptrdiff_t *target,
                                   ptrdiff_t *node, uint8_t *direction,
                                   ptrdiff_t dims[2], unsigned int order);

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

/**
   @brief Compute flow accumulation based on a weighted edge list

   Accumulates flow by summing contributing areas along flow paths.

   @param[out] acc The computed flow accumulation
   @parblock
   A pointer to a `float` array of size `dims[0]` x `dims[1]`
   @endparblock

   @param[in] source The source pixel for each edge
   @parblock
   A pointer to a `ptrdiff_t` array of size `edge_count`

   The source pixels must be in a topological order.
   @endparblock

   @param[in] target The target pixel for each edge
   @parblock
   A pointer to a `ptrdiff_t` array of size `edge_count`
   @endparblock

   @param[in] fraction The fraction of flow transported along each edge
   @parblock
   A pointer to a `float` array of size `edge_count`

   The fraction for each edge should be a value between zero and one,
   and the fractions for every edge with the same source pixel should
   sum to one.
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

   @param[in] edge_count The number of edges in the edge list

   @param[in] dims The dimensions of the arrays
   @parblock
   A pointer to a `ptrdiff_t` array of size 2

   The fastest changing dimension should be provided first. For column-major
   arrays, `dims = {nrows,ncols}`. For row-major arrays, `dims = {ncols,nrows}`.
   @endparblock
 */
TOPOTOOLBOX_API
void flow_accumulation_edgelist(float *acc, ptrdiff_t *source,
                                ptrdiff_t *target, float *fraction,
                                float *weights, ptrdiff_t edge_count,
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
void gradient8(float *output, float *dem, float cellsize, int use_mp,
               ptrdiff_t dims[2]);

/**
   @brief Performs a grayscale reconstruction of the `mask` image by the
   `marker` image using the sequential reconstruction algorithm of Vincent
   (1993).

   @param[out] marker: The marker array is updated with the result in-place.
   @param[in]  mask: The array holding the mask (for example dem).
   @param[in]  dims: An array specifying the dimensions of both used arrays.
                     It should contain two values: [rows, columns].
*/
TOPOTOOLBOX_API
void reconstruct(float *marker, float *mask, ptrdiff_t dims[2]);

/**
   @brief Performs a ggrayscale reconstruction using the hybrid
   algorithm of Vincent(1993).

   @param[out] marker: The marker array is updated with the result in-place.
   @param[in]  queue: A ptrdiff_t array of the same size as the marker and mask
   arrays.
   @param[in]  mask: The array holding the mask (for example dem).
   @param[in]  dims: An array specifying the dimensions of all used arrays.
                     It should contain two values: [rows, columns].
*/
TOPOTOOLBOX_API
void reconstruct_hybrid(float *marker, ptrdiff_t *queue, float *mask,
                        ptrdiff_t dims[2]);

/**
   @brief Integrate a `float` quantity over a stream network using
   trapezoidal integration.

   @param[out] integral The integrated output
   @parblock
   A pointer to a `float` array representing a node attribute list

   If the stream network has N nodes, this array should have a length
   N. This value must not be less than the largest value in either the
   `source` or `target` arrays.
   @endparblock

   @param[in] integrand The quantity to be integrated
   @parblock
   A pointer to a `float` array representing a node attribute list

   If the stream network has N nodes, this array should have a length
   N. This value must not be less than the largest value in either the
   `source` or `target` arrays.
   @endparblock

   @param[in] source The source node of each edge in the stream
                     network
   @parblock
   A pointer to a `ptrdiff_t` array of size `edge_count`

   The source nodes must be in topological order. The labels must
   correspond to the 0-based indices of the nodes in the
   node-attribute lists `integral` and `integrand`.
   @endparblock

   @param[in] target The target nodes of each edge in the stream
                     network
   @parblock
   A pointer to a `ptrdiff_t` array of size `edge_count`

   The labels must correspond to the 0-based indices of the nodes in
   the node-attribute lists `integral` and `integrand`.
   @endparblock

   @param[in] weight The weight assigned to each edge in the stream network
   @parblock
   A pointer to a `float` array of size `edge_count`

   For most applications of integration along the stream network, this
   will be the geometric distance between the source and target pixels
   in the desired units.
   @endparblock

   @param[in] edge_count The number of edges in the stream network
 */
TOPOTOOLBOX_API
void streamquad_trapz_f32(float *integral, float *integrand, ptrdiff_t *source,
                          ptrdiff_t *target, float *weight,
                          ptrdiff_t edge_count);

TOPOTOOLBOX_API
void streamquad_trapz_f64(double *integral, double *integrand,
                          ptrdiff_t *source, ptrdiff_t *target, float *weight,
                          ptrdiff_t edge_count);

/**
   @brief Upstream traversal in the (or, and) semiring

   Traverses the stream network in the upstream direction,
   accumulating the values in `output` with edge weights given by the
   `weights` array. The weights are applied with a bitwise and and
   accumulated with a bitwise or.

   ```
   for (e = (u,v) in edges:
     output[u] = output[u] | (output[v] & input[e]);
   ```

   @param[out] output The accumulated output
   @parblock
   A pointer to a `uint32_t` array representing a node attribute list

   This array should have a length equal to the number of nodes in the
   stream network.  This value must not be less than the largest value
   in either the `source` or `target` arrays.
   @endparblock

   @param[in] weights The edge weights
   @parblock
   A pointer to a `uint32_t` array of size `edge_count`
   @endparblock

   @param[in] source The source node of each edge in the stream
                     network
   @parblock
   A pointer to a `ptrdiff_t` array of size `edge_count`

   This array should have length `edge_count`.

   The source nodes must be in topological order. The labels must
   correspond to the 0-based indices of the node-attribute list
   `output`
   @endparblock

   @param[in] target The target nodes of each edge in the stream
                     network
   @parblock
   A pointer to a `ptrdiff_t` array of size `edge_count`

   The labels must correspond to the 0-based indices of the
   node-attribute lists `output`
   @endparblock

   @param[in] edge_count The number of edges in the stream network
 */
TOPOTOOLBOX_API
void traverse_up_u32_or_and(uint32_t *output, uint32_t *weights,
                            ptrdiff_t *source, ptrdiff_t *target,
                            ptrdiff_t edge_count);

/**
   @brief Downstream traversal in the Boolean ({0,1}, or, and) semiring

   Accumulates the input node attribute list downstream using bitwise or
   with multiplication given by bitwise and:

   for e = (u=>v) in edges:
     output[v] = output[v] | output[u] & input[e];

   @param[out] output The accumulated output
   @parblock
   A pointer to a `uint32_t` array representing a node attribute list

   This array should have a length equal to the number of nodes in the
   stream network.  This value must not be less than the largest value
   in either the `source` or `target` arrays.
   @endparblock

   @param[in] input The edge weights
   @parblock
   A pointer to a `uint32_t` array of size `edge_count`
   @endparblock

   @param[in] source The source node of each edge in the stream
                     network
   @parblock
   A pointer to a `ptrdiff_t` array of size `edge_count`

   The source nodes must be in topological order. The labels must
   correspond to the 0-based indices of the node-attribute list
   `output`.
   @endparblock

   @param[in] target The target nodes of each edge in the stream
                     network
   @parblock
   A pointer to a `ptrdiff_t` array of size `edge_count`

   The labels must correspond to the 0-based indices of the
   node-attribute list `output`.
   @endparblock

   @param[in] edge_count The number of edges in the stream network
 */
TOPOTOOLBOX_API
void traverse_down_u32_or_and(uint32_t *output, uint32_t *input,
                              ptrdiff_t *source, ptrdiff_t *target,
                              ptrdiff_t edge_count);
/**
   @brief Downstream traversal with max-plus

   Accumulates the edge weights in the `input` edge attribute list
   using a max-plus update:

   for (e = (u,v)) in edges:
     output[v] = max(output[v], output[u] + input[e]);

   With input giving the length of each edge and output initialized to
   zero, this will compute the longest path from the outlet to a
   source node.

   @param[out] output The accumulated output
   @parblock
   A pointer to a `float` array representing a node attribute list

   This array should have a length equal to the number of nodes in the
   stream network.  This value must not be less than the largest value
   in either the `source` or `target` arrays.
   @endparblock

   @param[in] input The edge weights
   @parblock
   A pointer to a `float` array of size `edge_count`
   @endparblock

   @param[in] source The source node of each edge in the stream
                     network
   @parblock
   A pointer to a `ptrdiff_t` array of size `edge_count`

   The source nodes must be in topological order. The labels must
   correspond to the 0-based indices of the node-attribute list
   `output`.
   @endparblock

   @param[in] target The target nodes of each edge in the stream
                     network
   @parblock
   A pointer to a `ptrdiff_t` array of size `edge_count`

   The labels must correspond to the 0-based indices of the
   node-attribute list `output`.
   @endparblock

   @param[in] edge_count The number of edges in the stream network
 */
TOPOTOOLBOX_API
void traverse_down_f32_max_add(float *output, float *input, ptrdiff_t *source,
                               ptrdiff_t *target, ptrdiff_t edge_count);

/**
   @brief Upstream traversal with max-plus

   Accumulates the edge weights in the `input` edge attribute list
   using a max-plus update.

   With input giving the length of each edge and output initialized to
   zero, this will compute the longest path from the outlet to a
   source node.

   @param[out] output The accumulated output
   @parblock
   A pointer to a `float` array representing a node attribute list

   This array should have a length equal to the number of nodes in the
   stream network.  This value must not be less than the largest value
   in either the `source` or `target` arrays.
   @endparblock

   @param[in] input The edge weights
   @parblock
   A pointer to a `float` array of size `edge_count`
   @endparblock

   @param[in] source The source node of each edge in the stream
                     network
   @parblock
   A pointer to a `ptrdiff_t` array of size `edge_count`

   The source nodes must be in topological order. The labels must
   correspond to the 0-based indices of the node-attribute list
   `output`.
   @endparblock

   @param[in] target The target nodes of each edge in the stream
                     network
   @parblock
   A pointer to a `ptrdiff_t` array of size `edge_count`

   The labels must correspond to the 0-based indices of the
   node-attribute list `output`.
   @endparblock

   @param[in] edge_count The number of edges in the stream network
 */
TOPOTOOLBOX_API
void traverse_up_f32_max_add(float *output, float *input, ptrdiff_t *source,
                             ptrdiff_t *target, ptrdiff_t edge_count);

/**
   @brief Downstream traversal with max-mul

   Accumulates the edge weights in the `input` edge attribute list
   using a max-mul update.

   @param[out] output The accumulated output
   @parblock
   A pointer to a `float` array representing a node attribute list

   This array should have a length equal to the number of nodes in the
   stream network.  This value must not be less than the largest value
   in either the `source` or `target` arrays.
   @endparblock

   @param[in] input The edge weights
   @parblock
   A pointer to a `float` array of size `edge_count`
   @endparblock

   @param[in] source The source node of each edge in the stream
                     network
   @parblock
   A pointer to a `ptrdiff_t` array of size `edge_count`

   The source nodes must be in topological order. The labels must
   correspond to the 0-based indices of the node-attribute list
   `output`.
   @endparblock

   @param[in] target The target nodes of each edge in the stream
                     network
   @parblock
   A pointer to a `ptrdiff_t` array of size `edge_count`

   The labels must correspond to the 0-based indices of the
   node-attribute list `output`.
   @endparblock

   @param[in] edge_count The number of edges in the stream network
 */
TOPOTOOLBOX_API
void traverse_down_f32_max_mul(float *output, float *input, ptrdiff_t *source,
                               ptrdiff_t *target, ptrdiff_t edge_count);

/**
   @brief Downstream traversal with max-mul and maximum tracking

   Accumulates the edge weights in the `input` edge attribute list
   using a max-mul update.

   The idx array can be used to track the index of the maximum
   value. It should be initialized with the index of each node in the
   array. If the upstream value multiplied by the input is greater
   than the downstream value, the upstream index is copied into the
   downstream index.

   @param[out] output The accumulated output
   @parblock
   A pointer to a `float` array representing a node attribute list

   This array should have a length equal to the number of nodes in the
   stream network.  This value must not be less than the largest value
   in either the `source` or `target` arrays.
   @endparblock

   @param[out] idx The index of the maximum value
   @parblock
   A pointer to an `int64_t` array representing a node attribute list

   This array should have a length equal to the number of nodes in the
   stream network.  This value must not be less than the largest value
   in either the `source` or `target` arrays.
   @endparblock

   @param[in] input The edge weights
   @parblock
   A pointer to a `float` array of size `edge_count`
   @endparblock

   @param[in] source The source node of each edge in the stream
                     network
   @parblock
   A pointer to a `ptrdiff_t` array of size `edge_count`

   The source nodes must be in topological order. The labels must
   correspond to the 0-based indices of the node-attribute list
   `output`.
   @endparblock

   @param[in] target The target nodes of each edge in the stream
                     network
   @parblock
   A pointer to a `ptrdiff_t` array of size `edge_count`

   The labels must correspond to the 0-based indices of the
   node-attribute list `output`.
   @endparblock

   @param[in] edge_count The number of edges in the stream network
 */
TOPOTOOLBOX_API
void traverse_down_f32_max_mul_arg(float *output, int64_t *idx, float *input,
                                   ptrdiff_t *source, ptrdiff_t *target,
                                   ptrdiff_t edge_count);

/**
   @brief Upstream traversal with max-mul

   Accumulates the edge weights in the `input` edge attribute list
   using a max-mul update.

   @param[out] output The accumulated output
   @parblock
   A pointer to a `float` array representing a node attribute list

   This array should have a length equal to the number of nodes in the
   stream network.  This value must not be less than the largest value
   in either the `source` or `target` arrays.
   @endparblock

   @param[in] input The edge weights
   @parblock
   A pointer to a `float` array of size `edge_count`
   @endparblock

   @param[in] source The source node of each edge in the stream
                     network
   @parblock
   A pointer to a `ptrdiff_t` array of size `edge_count`

   The source nodes must be in topological order. The labels must
   correspond to the 0-based indices of the node-attribute list
   `output`.
   @endparblock

   @param[in] target The target nodes of each edge in the stream
                     network
   @parblock
   A pointer to a `ptrdiff_t` array of size `edge_count`

   The labels must correspond to the 0-based indices of the
   node-attribute list `output`.
   @endparblock

   @param[in] edge_count The number of edges in the stream network
 */
TOPOTOOLBOX_API
void traverse_up_f32_max_mul(float *output, float *input, ptrdiff_t *source,
                             ptrdiff_t *target, ptrdiff_t edge_count);

/**
   @brief Upstream traversal with max-mul and maximum tracking

   Accumulates the edge weights in the `input` edge attribute list
   using a max-mul update.

   The idx array can be used to track the index of the maximum
   value. It should be initialized with the index of each node in the
   array. If the upstream value multiplied by the input is greater
   than the downstream value, the upstream index is copied into the
   downstream index.

   @param[out] output The accumulated output
   @parblock
   A pointer to a `float` array representing a node attribute list

   This array should have a length equal to the number of nodes in the
   stream network.  This value must not be less than the largest value
   in either the `source` or `target` arrays.
   @endparblock

   @param[out] idx The index of the maximum value
   @parblock
   A pointer to an `int64_t` array representing a node attribute list

   This array should have a length equal to the number of nodes in the
   stream network.  This value must not be less than the largest value
   in either the `source` or `target` arrays.
   @endparblock

   @param[in] input The edge weights
   @parblock
   A pointer to a `float` array of size `edge_count`
   @endparblock

   @param[in] source The source node of each edge in the stream
                     network
   @parblock
   A pointer to a `ptrdiff_t` array of size `edge_count`

   The source nodes must be in topological order. The labels must
   correspond to the 0-based indices of the node-attribute list
   `output`.
   @endparblock

   @param[in] target The target nodes of each edge in the stream
                     network
   @parblock
   A pointer to a `ptrdiff_t` array of size `edge_count`

   The labels must correspond to the 0-based indices of the
   node-attribute list `output`.
   @endparblock

   @param[in] edge_count The number of edges in the stream network
 */
TOPOTOOLBOX_API
void traverse_up_f32_max_mul_arg(float *output, int64_t *idx, float *input,
                                 ptrdiff_t *source, ptrdiff_t *target,
                                 ptrdiff_t edge_count);

/**
   @brief Downstream traversal with min-plus

   Accumulates the edge weights in the `input` edge attribute list
   using a min-plus update:

   for (e = (u,v)) in edges:
     output[v] = min(output[v], output[u] + input[e]);

   With input giving the length of each edge and output initialized to
   zero, this will compute the shortest path from the outlet to a
   source node.

   @param[out] output The accumulated output
   @parblock
   A pointer to a `float` array representing a node attribute list

   This array should have a length equal to the number of nodes in the
   stream network.  This value must not be less than the largest value
   in either the `source` or `target` arrays.
   @endparblock

   @param[in] input The edge weights
   @parblock
   A pointer to a `float` array of size `edge_count`
   @endparblock

   @param[in] source The source node of each edge in the stream
                     network
   @parblock
   A pointer to a `ptrdiff_t` array of size `edge_count`

   The source nodes must be in topological order. The labels must
   correspond to the 0-based indices of the node-attribute list
   `output`.
   @endparblock

   @param[in] target The target nodes of each edge in the stream
                     network
   @parblock
   A pointer to a `ptrdiff_t` array of size `edge_count`

   The labels must correspond to the 0-based indices of the
   node-attribute list `output`.
   @endparblock

   @param[in] edge_count The number of edges in the stream network
 */
TOPOTOOLBOX_API
void traverse_down_f32_min_add(float *output, float *input, ptrdiff_t *source,
                               ptrdiff_t *target, ptrdiff_t edge_count);

/**
   @brief Downstream traversal with multiply-add

   Accumulates the edge weights in the `input` edge attribute list
   using a multiply-add update:

   for (e = (u,v)) in edges:
     output[v] = output[v] + output[u] * input[e];

   With input giving the flow fraction of each edge and output
   initialized to ones, this will compute the flow accumulation.

   @param[out] output The accumulated output
   @parblock
   A pointer to a `float` array representing a node attribute list

   This array should have a length equal to the number of nodes in the
   stream network.  This value must not be less than the largest value
   in either the `source` or `target` arrays.
   @endparblock

   @param[in] input The edge weights
   @parblock
   A pointer to a `float` array of size `edge_count`
   @endparblock

   @param[in] source The source node of each edge in the stream
                     network
   @parblock
   A pointer to a `ptrdiff_t` array of size `edge_count`

   The source nodes must be in topological order. The labels must
   correspond to the 0-based indices of the node-attribute list
   `output`.
   @endparblock

   @param[in] target The target nodes of each edge in the stream
                     network
   @parblock
   A pointer to a `ptrdiff_t` array of size `edge_count`

   The labels must correspond to the 0-based indices of the
   node-attribute list `output`.
   @endparblock

   @param[in] edge_count The number of edges in the stream network
 */
TOPOTOOLBOX_API
void traverse_down_f32_add_mul(float *output, float *input, ptrdiff_t *source,
                               ptrdiff_t *target, ptrdiff_t edge_count);

/**
   @brief Compute the in- and outdegrees of each node in a graph

   @param[out] indegree The indegree of each node in the stream network
   @parblock
   A pointer to a `uint8_t` array representing a node attribute list

   This array should have a length equal to the number of nodes in the
   stream network.  This value must not be less than the largest value
   in either the `source` or `target` arrays.
   @endparblock

   @param[out] outdegree The outdegree of each node in the stream network
   @parblock
   A pointer to a `uint8_t` array representing a node attribute list

   This array should have a length equal to the number of nodes in the
   stream network.  This value must not be less than the largest value
   in either the `source` or `target` arrays.
   @endparblock

   @param[in] source The source node of each edge in the stream
                     network
   @parblock
   A pointer to a `ptrdiff_t` array of size `edge_count`

   The source nodes must be in topological order. The labels must
   correspond to the 0-based indices of the node-attribute list
   `indegree` and `outdegree`.
   @endparblock

   @param[in] target The target nodes of each edge in the stream
                     network
   @parblock
   A pointer to a `ptrdiff_t` array of size `edge_count`

   The labels must correspond to the 0-based indices of the
   node-attribute lists `indegree` and `outdegree`.
   @endparblock

   @param[in] node_count The number of nodes in the stream network
   @param[in] edge_count The number of edges in the stream network
 */
TOPOTOOLBOX_API
void edgelist_degree(uint8_t *indegree, uint8_t *outdegree, ptrdiff_t *source,
                     ptrdiff_t *target, ptrdiff_t node_count,
                     ptrdiff_t edge_count);

/**
   @brief Compute stream order of Stream Object using Strahler method

   Accumulates the edge weights in the `input` edge attribute list
   using a multiply-add update:

   for (e = (u,v)) in edges:
     output[v] = output[v] + output[u] * input[e];

   With input giving the flow fraction of each edge and output
   initialized to ones, this will compute the flow accumulation.

   @param[out] output The accumulated output
   @parblock
   A pointer to a `float` array representing a node attribute list

   This array should have a length equal to the number of nodes in the
   stream network.  This value must not be less than the largest value
   in either the `source` or `target` arrays.
   @endparblock

   @param[in] input The edge weights
   @parblock
   A pointer to a `float` array of size `edge_count`
   @endparblock

   @param[in] source The source node of each edge in the stream
                     network
   @parblock
   A pointer to a `ptrdiff_t` array of size `edge_count`

   The source nodes must be in topological order. The labels must
   correspond to the 0-based indices of the node-attribute list
   `output`.
   @endparblock

   @param[in] target The target nodes of each edge in the stream
                     network
   @parblock
   A pointer to a `ptrdiff_t` array of size `edge_count`

   The labels must correspond to the 0-based indices of the
   node-attribute list `output`.
   @endparblock

   @param[in] edge_count The number of edges in the stream network
 */

TOPOTOOLBOX_API
void traverse_down_f32_strahler(float *output, float *input, ptrdiff_t *source,
                                ptrdiff_t *target, ptrdiff_t edge_count);

/**
   @brief Propagate `float` values upstream

   The values in the `data` node attribute list are copied to their
   upstream neighbors according to the edges contained in the `source`
   and `target` arrays.

   Nodes with no downstream neighbors should be initialized in the
   `data` array with the values to be propagated upstream. The
   upstream values will be overwritten.

   @param[in, out] data The data array
   @parblock
   A pointer to an array of the desired type representing a node attribute list.

   This must be at least as large as the largest index in either the
   `source` or `target` arrays.
   @endparblock

   @param[in] source The source node of each edge in the stream
                     network
   @parblock
   A pointer to a `ptrdiff_t` array of size `edge_count`

   The source nodes must be in topological order. The labels must
   correspond to the 0-based indices of the node-attribute list
   `data`.
   @endparblock

   @param[in] target The target nodes of each edge in the stream
                     network
   @parblock
   A pointer to a `ptrdiff_t` array of size `edge_count`

   The labels must correspond to the 0-based indices of the
   node-attribute list `data`.
   @endparblock

   @param[in] edge_count The number of edges in the stream network
 */
TOPOTOOLBOX_API void propagatevaluesupstream_f32(float *data, ptrdiff_t *source,
                                                 ptrdiff_t *target,
                                                 ptrdiff_t edge_count);

/**
   @brief Propagate `double` values upstream

   @copydetails propagatevaluesupstream_f32()
 */
TOPOTOOLBOX_API void propagatevaluesupstream_f64(double *data,
                                                 ptrdiff_t *source,
                                                 ptrdiff_t *target,
                                                 ptrdiff_t edge_count);

/**
   @brief Propagate `uint8_t` values upstream

   @copydetails propagatevaluesupstream_f32()
 */
TOPOTOOLBOX_API void propagatevaluesupstream_u8(uint8_t *data,
                                                ptrdiff_t *source,
                                                ptrdiff_t *target,
                                                ptrdiff_t edge_count);

/**
   @brief Propagate `uint32_t` values upstream

   @copydetails propagatevaluesupstream_f32()
 */
TOPOTOOLBOX_API void propagatevaluesupstream_u32(uint32_t *data,
                                                 ptrdiff_t *source,
                                                 ptrdiff_t *target,
                                                 ptrdiff_t edge_count);

/**
   @brief Propagate `uint64_t` values upstream

   @copydetails propagatevaluesupstream_f32()
 */
TOPOTOOLBOX_API void propagatevaluesupstream_u64(uint64_t *data,
                                                 ptrdiff_t *source,
                                                 ptrdiff_t *target,
                                                 ptrdiff_t edge_count);

/**
   @brief Propagate `int8_t` values upstream

   @copydetails propagatevaluesupstream_f32()
 */
TOPOTOOLBOX_API void propagatevaluesupstream_i8(int8_t *data, ptrdiff_t *source,
                                                ptrdiff_t *target,
                                                ptrdiff_t edge_count);

/**
   @brief Propagate `int32_t` values upstream

   @copydetails propagatevaluesupstream_f32()
 */
TOPOTOOLBOX_API void propagatevaluesupstream_i32(int32_t *data,
                                                 ptrdiff_t *source,
                                                 ptrdiff_t *target,
                                                 ptrdiff_t edge_count);

/**
   @brief Propagate `int64_t` values upstream

   @copydetails propagatevaluesupstream_f32()
 */
TOPOTOOLBOX_API void propagatevaluesupstream_i64(int64_t *data,
                                                 ptrdiff_t *source,
                                                 ptrdiff_t *target,
                                                 ptrdiff_t edge_count);
/*
  Graphflood
*/

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
   @param[in]  dim: [rows,columns] if row major and [columns, rows] if
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
                                    bool D8, GF_FLOAT step);

/**
   @brief Fills the depressions in place in the topography using Priority
   Floods Barnes (2014, modified to impose a minimal slope)

   @param[inout]  topo: array of surface elevation
   @param[in]     BCs: codes for boundary conditions and no data
   management, see gf_utils.h or examples for the meaning
   @param[in]     dim: [rows,columns] if row major and [columns, rows] if
   column major
   @param[in]     D8: true for topology including cardinals + diagonals,
   @param[in]     step: delta_Z to apply minimum elevation increase and avoid
   flats, false for cardinals only
*/
TOPOTOOLBOX_API
void compute_priority_flood(GF_FLOAT *topo, uint8_t *BCs, GF_UINT *dim, bool D8,
                            GF_FLOAT step);

/**
   @brief Fills the depressions in place in the topography using Priority
   Floods Barnes (2014, modified to impose a minimal slope) This variant
   computes the topological order on the go (slightly slower as it uses a
   priority queue for all the nodes including in depressions)

   @param[inout]  topo: array of surface elevation
   @param[in]  Stack: topologically ordered list of nodes, from the
   baselevel to the sources
   @param[in]     BCs: codes for boundary conditions and no data
   management, see gf_utils.h or examples for the meaning
   @param[in]     dim: [rows,columns] if row major and [columns, rows] if
   column major
   @param[in]     D8: true for topology including cardinals + diagonals,
   @param[in]     step: delta_Z to apply minimum elevation increase and avoid
   flats, false for cardinals only
*/
TOPOTOOLBOX_API
void compute_priority_flood_plus_topological_ordering(GF_FLOAT *topo,
                                                      GF_UINT *Stack,
                                                      uint8_t *BCs,
                                                      GF_UINT *dim, bool D8,
                                                      GF_FLOAT step);

/**
   @brief Accumulate single flow drainage area downstream from a calculated
   graphflood single flow graph
   @param[out] output: the field of drainage area
   @param[in]  Sreceivers: array of steepest receiver vectorised index
   @param[in]  Stack: topologically ordered list of nodes, from the
   baselevel to the sources
   @param[in]  dim: [rows,columns] if row major and [columns, rows] if
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
   @param[in]  dim: [rows,columns] if row major and [columns, rows] if
   column major
   @param[in]  dx: spatial step
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
   @param[in]     dim: [rows,columns] if row major and [columns, rows] if
   column major
   @param[in]     dt: time step
   @param[in]     dx: spatial step
   @param[in]     SFD: single flow direction if True, multiple flow if
   false
   @param[in]     D8: true for topology including cardinals + diagonals,
   false for cardinals only
   @param[in]     N_iterations: number of iterations of the flooding algorithm
   @param[in]     step: delta_Z to apply minimum elevation increase and avoid
   flats, false for cardinals only
*/
TOPOTOOLBOX_API
void graphflood_full(GF_FLOAT *Z, GF_FLOAT *hw, uint8_t *BCs,
                     GF_FLOAT *Precipitations, GF_FLOAT *manning, GF_UINT *dim,
                     GF_FLOAT dt, GF_FLOAT dx, bool SFD, bool D8,
                     GF_UINT N_iterations, GF_FLOAT step);

/**
   @brief Calculate steady-state flow metrics from topography and water depths
   using the GraphFlood algorithm as described in Gailleton et al., 2024.
   From input fields of topography, water depth, Manning's friction coefficient
   and precipitation rates, computes various flow metrics assuming steady-state
   conditions with multiple flow direction routing.

   @param[in]     Z: surface topography [m]
   @param[in]     hw: field of water depth [m]
   @param[in]     BCs: codes for boundary conditions and no data management,
                  see gf_utils.h or examples for the meaning
   @param[in]     Precipitations: precipitation rates [m/s]
   @param[in]     manning: Manning's friction coefficient [s/m^(1/3)]
   @param[in]     dim: [rows,columns] if row major and [columns, rows] if
                  column major
   @param[in]     dx: spatial step [m]
   @param[in]     D8: true for topology including cardinals + diagonals,
                  false for cardinals only
   @param[in]     step: delta_Z to apply minimum elevation increase and avoid
                  flats during priority flooding
   @param[out]    Qi: input discharge based on drainage area accumulation [m³/s]
   @param[out]    Qo: output discharge calculated from Manning's equation [m³/s]
   @param[out]    qo: discharge per unit width [m²/s]
   @param[out]    u: flow velocity [m/s]
   @param[out]    Sw: water surface slope [-]

   @note At convergence, Qi should approximately equal Qo (mass conservation).
   Differences may indicate numerical issues, local instabilities, or physical
   processes such as ponding in lakes or depressions.
*/
TOPOTOOLBOX_API
void graphflood_metrics(GF_FLOAT *Z, GF_FLOAT *hw, uint8_t *BCs,
                        GF_FLOAT *Precipitations, GF_FLOAT *manning,
                        GF_UINT *dim, GF_FLOAT dx, bool D8, GF_FLOAT step,
                        GF_FLOAT *Qi, GF_FLOAT *Qo, GF_FLOAT *qo, GF_FLOAT *u,
                        GF_FLOAT *Sw);

/**
   @brief Run dynamic induced graph flood simulation using wavefront propagation
   from specified input discharge locations. Processes cells in descending
   hydraulic elevation order, dynamically building the flow graph as it
   propagates downstream.

   @param[in]     Z: surface topography [m]
   @param[inout]  hw: field of water depth [m]
   @param[in]     BCs: codes for boundary conditions and no data management
   @param[in]     Precipitations: precipitation rates [m/s]
   @param[in]     manning: Manning's friction coefficient [s/m^(1/3)]
   @param[in]     input_Qw: input discharge at specific cells [m³/s]
                  (cells with value > 0 are used as starting points)
   @param[inout]  Qwin: accumulated input discharge array [m³/s]
                  (pre-allocated array of size dim[0] * dim[1], reinitialized
   each iteration)
   @param[in]     dim: [rows,columns] if row major and [columns, rows] if
                  column major
   @param[in]     dt: time step [s]
   @param[in]     dx: spatial step [m]
   @param[in]     D8: true for topology including cardinals + diagonals,
                  false for cardinals only
   @param[in]     N_iterations: number of iterations of the flooding algorithm
*/
TOPOTOOLBOX_API
void graphflood_dynamic_graph(GF_FLOAT *Z, GF_FLOAT *hw, uint8_t *BCs,
                              GF_FLOAT *Precipitations, GF_FLOAT *manning,
                              GF_FLOAT *input_Qw, GF_FLOAT *Qwin, GF_UINT *dim,
                              GF_FLOAT dt, GF_FLOAT dx, bool D8,
                              GF_UINT N_iterations);

/**
   @brief Compute input discharge array for dynamic graph from drainage area
   threshold. Identifies channel heads where drainage area crosses a threshold
   and assigns precipitation-weighted discharge at those entry points.

   @param[out]    input_Qw: input discharge array to fill [m³/s]
                  (pre-allocated array of size dim[0] * dim[1])
   @param[in]     Z: surface topography [m]
   @param[in]     hw: field of water depth [m]
   @param[in]     BCs: codes for boundary conditions and no data management
   @param[in]     Precipitations: precipitation rates [m/s]
   @param[in]     area_threshold: drainage area threshold for entry points [m²]
   @param[in]     dim: [rows,columns] if row major and [columns, rows] if
                  column major
   @param[in]     dx: spatial step [m]
   @param[in]     D8: true for topology including cardinals + diagonals,
                  false for cardinals only
   @param[in]     step: delta_Z to apply minimum elevation increase and avoid
                  flats during priority flooding
*/
TOPOTOOLBOX_API
void compute_input_Qw_from_area_threshold(GF_FLOAT *input_Qw, GF_FLOAT *Z,
                                          GF_FLOAT *hw, uint8_t *BCs,
                                          GF_FLOAT *Precipitations,
                                          GF_FLOAT area_threshold, GF_UINT *dim,
                                          GF_FLOAT dx, bool D8, GF_FLOAT step);

/**
   @brief Label drainage basins based on the flow directions provided
   by a topologically sorted edge list.

   @param[out] basins The drainage basin label
   @parblock
   A pointer to a `ptrdiff_t` array of size `dims[0]` x `dims[1]`
   @endparblock

   @param[out] source The source pixel for each edge
   @parblock
   A pointer to a `ptrdiff_t` array of size `edge_count`
   @endparblock

   @param[in] target The target pixel for each edge
   @parblock
   A pointer to a `ptrdiff_t` array of size `edge_count`
   @endparblock

   @param[in] edge_count The number of edges in the flow network
   @parblock
   A ptrdiff_t representing the length of the `source` and `target` arrays
   @endparblock

   @param[in] dims The dimensions of the arrays
   @parblock
   A pointer to a `ptrdiff_t` array of size 2

   The fastest changing dimension should be provided first. For column-major
   arrays, `dims = {nrows,ncols}`. For row-major arrays, `dims = {ncols,nrows}`.
   @endparblock
 */
TOPOTOOLBOX_API
void drainagebasins(ptrdiff_t *basins, ptrdiff_t *source, ptrdiff_t *target,
                    ptrdiff_t edge_count, ptrdiff_t dims[2]);

/**
   @brief Compute the gradient of a DEM using a second-order finite difference
approximation

   @param[out] p0 The gradient in the first dimension
   @parblock
   A pointer to a float array of size `dims[0] * dims[1]`

   `p0` will contain the gradient in the first dimension of the array
   as specified by `dims`. Whether this corresponds to the positive or
   negative x or y direction depends on the memory layout and
   coordinate system of the array.
   @endparblock

   @param[out] p1 The gradient in the second dimension
   @parblock
   A pointer to a float array of size `dims[0] * dims[1]`

   `p1` will contain the gradient in the second dimension of the array
   as specified by `dims`. Whether this corresponds to the positive or
   negative x or y direction depends on the memory layout and
   coordinate system of the array.
   @endparblock

   @param[in] dem The input digital elevation model
   @parblock
   A pointer to a float array of size `dims[0] * dims[1]`
   @endparblock

   @param[in] cellsize The horizontal resolution of the digital elevation model

   @param[in] dims The dimensions of the arrays
   @parblock
   A pointer to a `ptrdiff_t` array of size 2

   The fastest changing dimension should be provided first. For column-major
   arrays, `dims = {nrows,ncols}`. For row-major arrays, `dims = {ncols,nrows}`.
   @endparblock
 */
TOPOTOOLBOX_API
void gradient_secondorder(float *p0, float *p1, float *dem, float cellsize,
                          ptrdiff_t dims[2]);

/**
   @brief Compute a hillshade of the supplied digital elevation model

   @param[out] output The output hillshade
   @parblock
   A pointer to a float array of size `dims[0] * dims[1]`.
   @endparblock

   @param[out] dx The surface gradient in the first dimension
   @parblock
   A pointer to a float array of size `dims[0] * dims[1]`.
   @endparblock

   @param[out] dy The surface gradient in the second dimension
   @parblock
   A pointer to a float array of size `dims[0] * dims[1]`.
   @endparblock

   @param[in] dem The input digital elevation model
   @parblock
   A pointer to a float array of size `dims[0] * dims[1]`.
   @endparblock

   @param[in] azimuth The azimuth angle of the light source
   @parblock
   `hillshade` expects the azimuth angle to be provided in radians
   from the first dimension of the grid towards the second dimension
   of the grid. This is most likely not identical to the geographic
   azimuth traditionally measured clockwise from north. Callers should
   take care to rotate their desired azimuth angle into the grid
   coordinate system depending on the georeferencing of the grid and
   the memory layout of the array.
   @endparblock

   @param[in] altitude The altitude angle of the light source
   @parblock
   The altitude angle must be provided in radians above the horizon.
   @endparblock

   @param[in] cellsize The spatial resolution of the DEM

   @param[in] dims The dimensions of the arrays
   @parblock
   A pointer to a `ptrdiff_t` array of size 2

   The fastest changing dimension should be provided first. For column-major
   arrays, `dims = {nrows,ncols}`. For row-major arrays, `dims = {ncols,nrows}`.
   @endparblock
 */
TOPOTOOLBOX_API
void hillshade(float *output, float *dx, float *dy, float *dem, float azimuth,
               float altitude, float cellsize, ptrdiff_t dims[2]);

/**
   @brief Compute a hillshade of the supplied digital elevation model

   The output of `hillshade_fused` should be identical to that of
   `hillshade`, but the gradient, surface normal and shading loops
   have been fused to remove the need for intermediate arrays. This
   causes a slight performance degradation, but significantly reduces
   the amount of memory needed.

   @param[out] output The output hillshade
   @parblock
   A pointer to a float array of size `dims[0] * dims[1]`.
   @endparblock

   @param[in] dem The input digital elevation model
   @parblock
   A pointer to a float array of size `dims[0] * dims[1]`.
   @endparblock

   @param[in] azimuth The azimuth angle of the light source
   @parblock
   `hillshade` expects the azimuth angle to be provided in radians
   from the first dimension of the grid towards the second dimension
   of the grid. This is most likely not identical to the geographic
   azimuth traditionally measured clockwise from north. Callers should
   take care to rotate their desired azimuth angle into the grid
   coordinate system depending on the georeferencing of the grid and
   the memory layout of the array.
   @endparblock

   @param[in] altitude The altitude angle of the light source
   @parblock
   The altitude angle must be provided in radians above the horizon.
   @endparblock

   @param[in] cellsize The spatial resolution of the DEM

   @param[in] dims The dimensions of the arrays
   @parblock
   A pointer to a `ptrdiff_t` array of size 2

   The fastest changing dimension should be provided first. For column-major
   arrays, `dims = {nrows,ncols}`. For row-major arrays, `dims = {ncols,nrows}`.
   @endparblock
 */
TOPOTOOLBOX_API
void hillshade_fused(float *output, float *dem, float azimuth, float altitude,
                     float cellsize, ptrdiff_t dims[2]);

/**
   @brief Compute the lower convex envelope of a stream profile

   @param[inout] elevation A node attribute list of elevations
   @parblock
   A pointer to a float array representing a node attribute list.

   `elevation` should be initialized with the elevation of the stream
   profile whose lower convex envelope will be computed. It will be
   modified in place.
   @endparblock

   @param[in] knickpoints A logical node attribute list specifying knickpoints
   @parblock

   A pointer to a uint8_t array representing a node attribute list.

   The resulting elevation profile will be nonconvex only at nodes
   whose value in this array is nonzero.
   @endparblock

   @param[in] distance A node attribute list containing upstream distances
   @parblock

   A pointer to a float array representing a node attribute list. Each
   entry should contain the distance upstream from an outlet of the
   corresponding node.
   @endparblock

   @param[inout] ix A node attribute list used as an intermediate array
   @parblock

   A pointer to a ptrdiff_t array representing a node attribute
   list. It is used as an intermediate array and is initialized as
   needed within `lowerenv`.

   @endparblock

   @param[inout] onenvelope A node attribute list used as an intermediate array
   @parblock

   A pointer to a uint8_t array representing a node attribute list. It
   is used an an intermediate array and is initialized as needed
   within `lowerenv`.

   @endparblock

   @param[in] source The source node of each edge in the stream network
   @parblock

   A pointer to a ptrdiff_t array representing an edge attribute list.

   @endparblock

   @param[in] target The target node of each edge in the stream network
   @parblock

   A pointer to a ptrdiff_t array representing an edge attribute list.

   @endparblock

   @param[in] edge_count The number of edges in the stream network
   @param[in] node_count The number of nodes in the stream network
 */
TOPOTOOLBOX_API
void lowerenv(float *elevation, uint8_t *knickpoints, float *distance,
              ptrdiff_t *ix, uint8_t *onenvelope, ptrdiff_t *source,
              ptrdiff_t *target, ptrdiff_t edge_count, ptrdiff_t node_count);

/**
   @brief Compute number of bins for swath profile

   @details
   Helper function to compute the number of bins needed for a swath
   profile given the half-width and bin resolution.

   @param[in] half_width Half-width of the swath in meters
   @param[in] bin_resolution Spacing between bin centers in meters

   @return Number of bins from -half_width to +half_width
 */
TOPOTOOLBOX_API
ptrdiff_t swath_compute_nbins(float half_width, float bin_resolution);

/**
   @brief Compute distance map from DEM pixels to track

   @details
   For each DEM pixel, computes the signed perpendicular distance to the
   nearest track segment. Negative distances indicate pixels to the left
   of the track direction, positive to the right.

   Uses Euclidean distance in pixel space scaled by cellsize.

   @param[out] distance Distance map
   @parblock
   A pointer to a `float` array of size `dims[0]` x `dims[1]`

   Each element contains the signed perpendicular distance in meters from
   the pixel to the nearest track segment.
   @endparblock

   @param[out] nearest_segment Nearest segment map (optional)
   @parblock
   A pointer to a `ptrdiff_t` array of size `dims[0]` x `dims[1]`, or NULL

   If provided, each element will contain the index of the nearest track
   segment (0 to n_track_points-2). Pass NULL to skip this output.
   @endparblock

   @param[in] track_i Track point coordinates in fast dimension
   @parblock
   A pointer to a `float` array of size `n_track_points`

   Contains pixel coordinates (can be sub-pixel values) in the fast-changing
   dimension. For column-major arrays, these are row indices. For row-major
   arrays, these are column indices.
   @endparblock

   @param[in] track_j Track point coordinates in slow dimension
   @parblock
   A pointer to a `float` array of size `n_track_points`

   Contains pixel coordinates (can be sub-pixel values) in the slow-changing
   dimension. For column-major arrays, these are column indices. For row-major
   arrays, these are row indices.
   @endparblock

   @param[in] n_track_points Number of points in the track
   @parblock
   Must be at least 2 to form segments.
   @endparblock

   @param[in] dims The dimensions of the DEM arrays
   @parblock
   A pointer to a `ptrdiff_t` array of size 2

   The fastest changing dimension should be provided first. For column-major
   arrays, `dims = {nrows,ncols}`. For row-major arrays, `dims = {ncols,nrows}`.
   @endparblock

   @param[in] cellsize Physical size of one pixel in meters
 */
TOPOTOOLBOX_API
void swath_distance_map(float *distance, ptrdiff_t *nearest_segment,
                        const float *track_i, const float *track_j,
                        ptrdiff_t n_track_points, ptrdiff_t dims[2],
                        float cellsize);

/**
   @brief Compute binned swath profile perpendicular to track

   @details
   Aggregates DEM elevations by perpendicular distance to the track,
   creating a single averaged cross-sectional profile. All pixels along
   the entire track length are binned by their perpendicular distance,
   producing statistics for each distance bin.

   This is useful for analyzing typical valley cross-sections or ridge
   profiles averaged over a long track.

   @param[out] bin_distances Distance of each bin center from track
   @parblock
   A pointer to a `float` array of size `n_bins`

   Contains the distance in meters from the track for each bin center.
   Negative values are to the left of the track, positive to the right.
   @endparblock

   @param[out] bin_means Mean elevation in each bin
   @parblock
   A pointer to a `float` array of size `n_bins`

   Contains the mean elevation in meters for all pixels in each distance bin.
   @endparblock

   @param[out] bin_stddevs Standard deviation of elevation in each bin
   @parblock
   A pointer to a `float` array of size `n_bins`

   Contains the standard deviation of elevations in meters for each bin.
   @endparblock

   @param[out] bin_mins Minimum elevation in each bin
   @parblock
   A pointer to a `float` array of size `n_bins`

   Contains the minimum elevation in meters for each bin.
   @endparblock

   @param[out] bin_maxs Maximum elevation in each bin
   @parblock
   A pointer to a `float` array of size `n_bins`

   Contains the maximum elevation in meters for each bin.
   @endparblock

   @param[out] bin_counts Number of pixels in each bin
   @parblock
   A pointer to a `ptrdiff_t` array of size `n_bins`

   Contains the count of valid (non-NaN) pixels contributing to each bin.
   @endparblock

   @param[out] bin_medians Median elevation in each bin
   @parblock
   A pointer to a `float` array of size `n_bins`, or NULL to skip

   Contains the median (50th percentile) elevation in meters for each bin.
   @endparblock

   @param[out] bin_q1 First quartile elevation in each bin
   @parblock
   A pointer to a `float` array of size `n_bins`, or NULL to skip

   Contains the 25th percentile elevation in meters for each bin.
   @endparblock

   @param[out] bin_q3 Third quartile elevation in each bin
   @parblock
   A pointer to a `float` array of size `n_bins`, or NULL to skip

   Contains the 75th percentile elevation in meters for each bin.
   @endparblock

   @param[in] percentile_list List of percentiles to compute
   @parblock
   A pointer to an `int` array containing percentile values (0-100), or NULL

   Each value should be between 0 and 100. Pass NULL if no custom percentiles needed.
   @endparblock

   @param[in] n_percentiles Number of percentiles in percentile_list
   @parblock
   Number of elements in percentile_list. Ignored if percentile_list is NULL.
   @endparblock

   @param[out] bin_percentiles Custom percentiles for each bin
   @parblock
   A pointer to a `float` array of size `n_bins` x `n_percentiles`, or NULL

   Layout: bin_percentiles[bin * n_percentiles + p] contains the p-th percentile
   for the given bin. Both percentile_list and this must be non-NULL to compute.
   @endparblock

   @param[in] dem The input DEM
   @parblock
   A pointer to a `float` array of size `dims[0]` x `dims[1]`

   NaN values in the DEM are ignored in statistics computation.
   @endparblock

   @param[in] track_i Track point coordinates in fast dimension
   @parblock
   A pointer to a `float` array of size `n_track_points`
   @endparblock

   @param[in] track_j Track point coordinates in slow dimension
   @parblock
   A pointer to a `float` array of size `n_track_points`
   @endparblock

   @param[in] n_track_points Number of points in the track (must be >= 2)

   @param[in] dims The dimensions of the DEM arrays
   @parblock
   A pointer to a `ptrdiff_t` array of size 2
   @endparblock

   @param[in] cellsize Physical size of one pixel in meters

   @param[in] half_width Half-width of the swath in meters
   @parblock
   Only pixels within this perpendicular distance from the track are included.
   @endparblock

   @param[in] bin_resolution Spacing between bin centers in meters

   @param[in] n_bins Number of bins (use swath_compute_nbins to calculate)

   @param[in] normalize Normalize elevations to track elevation
   @parblock
   If non-zero, elevations are normalized by subtracting the mean track
   elevation before binning, then adding it back to the output. This
   highlights relative elevation patterns perpendicular to the track.
   @endparblock
 */
TOPOTOOLBOX_API
void swath_transverse(float *bin_distances, float *bin_means,
                      float *bin_stddevs, float *bin_mins, float *bin_maxs,
                      ptrdiff_t *bin_counts, float *bin_medians, float *bin_q1,
                      float *bin_q3, const int *percentile_list,
                      ptrdiff_t n_percentiles, float *bin_percentiles,
                      const float *dem, const float *track_i,
                      const float *track_j, ptrdiff_t n_track_points,
                      ptrdiff_t dims[2], float cellsize, float half_width,
                      float bin_resolution, ptrdiff_t n_bins, int normalize);


/**
   @brief Compute per-point swath profile along track

   @details
   For each track point, computes statistics of DEM elevations within a
   perpendicular swath at that point. This shows how the swath profile
   changes along the track length.

   This is useful for analyzing how valley depth, width, or asymmetry
   varies along a river or ridge line.

   @param[out] point_means Mean elevation for each track point
   @parblock
   A pointer to a `float` array of size `n_track_points`

   Contains the mean elevation in meters of all pixels within the swath
   at each track point.
   @endparblock

   @param[out] point_stddevs Standard deviation for each track point
   @parblock
   A pointer to a `float` array of size `n_track_points`

   Contains the standard deviation of elevations in meters for each point's swath.
   @endparblock

   @param[out] point_mins Minimum elevation for each track point
   @parblock
   A pointer to a `float` array of size `n_track_points`

   Contains the minimum elevation in meters within each point's swath.
   @endparblock

   @param[out] point_maxs Maximum elevation for each track point
   @parblock
   A pointer to a `float` array of size `n_track_points`

   Contains the maximum elevation in meters within each point's swath.
   @endparblock

   @param[out] point_counts Number of pixels for each track point
   @parblock
   A pointer to a `ptrdiff_t` array of size `n_track_points`

   Contains the count of valid pixels within each point's swath.
   @endparblock

   @param[out] point_medians Median elevation for each track point
   @parblock
   A pointer to a `float` array of size `n_track_points`, or NULL to skip

   Contains the median (50th percentile) elevation in meters for each point.
   @endparblock

   @param[out] point_q1 First quartile elevation for each track point
   @parblock
   A pointer to a `float` array of size `n_track_points`, or NULL to skip

   Contains the 25th percentile elevation in meters for each point.
   @endparblock

   @param[out] point_q3 Third quartile elevation for each track point
   @parblock
   A pointer to a `float` array of size `n_track_points`, or NULL to skip

   Contains the 75th percentile elevation in meters for each point.
   @endparblock

   @param[in] percentile_list List of percentiles to compute
   @parblock
   A pointer to an `int` array containing percentile values (0-100), or NULL

   Each value should be between 0 and 100. Pass NULL if no custom percentiles needed.
   @endparblock

   @param[in] n_percentiles Number of percentiles in percentile_list
   @parblock
   Number of elements in percentile_list. Ignored if percentile_list is NULL.
   @endparblock

   @param[out] point_percentiles Custom percentiles for each track point
   @parblock
   A pointer to a `float` array of size `n_track_points` x `n_percentiles`, or NULL

   Layout: point_percentiles[point * n_percentiles + p] contains the p-th percentile
   for the given point. Both percentile_list and this must be non-NULL to compute.
   @endparblock

   @param[in] dem The input DEM
   @parblock
   A pointer to a `float` array of size `dims[0]` x `dims[1]`
   @endparblock

   @param[in] track_i Track point coordinates in fast dimension
   @parblock
   A pointer to a `float` array of size `n_track_points`
   @endparblock

   @param[in] track_j Track point coordinates in slow dimension
   @parblock
   A pointer to a `float` array of size `n_track_points`
   @endparblock

   @param[in] n_track_points Number of points in the track (must be >= 2)

   @param[in] dims The dimensions of the DEM arrays
   @parblock
   A pointer to a `ptrdiff_t` array of size 2
   @endparblock

   @param[in] cellsize Physical size of one pixel in meters

   @param[in] half_width Half-width of the swath in meters
   @parblock
   Defines the perpendicular extent of the swath on each side of the track.
   @endparblock

   @param[in] binning_distance Along-track binning distance in meters
   @parblock
   For each track point, all track points within ±binning_distance (measured
   along the track) form a local sub-track. Pixels are gathered using
   perpendicular distance to this sub-track, same as the transverse method.

   If binning_distance < cellsize: the sub-track reduces to a single
   orthogonal line at the track point (minimal along-track extent).

   If binning_distance >= cellsize: the sub-track covers a longer portion,
   gathering more pixels and producing smoother statistics.
   @endparblock
 */
TOPOTOOLBOX_API
void swath_longitudinal(float *point_means, float *point_stddevs,
                         float *point_mins, float *point_maxs,
                         ptrdiff_t *point_counts, float *point_medians,
                         float *point_q1, float *point_q3,
                         const int *percentile_list, ptrdiff_t n_percentiles,
                         float *point_percentiles, const float *dem,
                         const float *track_i, const float *track_j,
                         ptrdiff_t n_track_points, ptrdiff_t dims[2],
                         float cellsize, float half_width,
                         float binning_distance);

/**
   @brief Get pixel coordinates associated with a single track point.

   @details
   Uses the same sub-track + perpendicular distance logic as swath_longitudinal.
   For the given track point, finds all grid pixels within half_width of the
   local sub-track (defined by ±binning_distance along the track).

   Caller allocates pixels_i and pixels_j large enough (safe upper bound:
   dims[0] * dims[1], though actual count will be much smaller). The function
   returns the number of pixels written.

   @param[out] pixels_i Fast-dimension coordinates of associated pixels (node indices)
   @param[out] pixels_j Slow-dimension coordinates of associated pixels (node indices)
   @param[in] track_i Track point coordinates in fast dimension (pixel indices)
   @param[in] track_j Track point coordinates in slow dimension (pixel indices)
   @param[in] n_track_points Number of points in the track (must be >= 2)
   @param[in] point_index Index into track_i/track_j for the query point
   @param[in] dims Dimensions of the DEM grid (pixel counts)
   @param[in] cellsize Physical size of one pixel (meters/pixel)
   @param[in] half_width Perpendicular half-width of the swath (meters)
   @param[in] binning_distance Along-track binning distance (meters)

   @return Number of pixels written to pixels_i and pixels_j
 */
TOPOTOOLBOX_API
ptrdiff_t swath_get_point_pixels(ptrdiff_t *pixels_i, ptrdiff_t *pixels_j,
                                  const float *track_i, const float *track_j,
                                  ptrdiff_t n_track_points,
                                  ptrdiff_t point_index, ptrdiff_t dims[2],
                                  float cellsize, float half_width,
                                  float binning_distance);
#endif  // TOPOTOOLBOX_H
