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

#include <stddef.h>
#include <stdint.h>

/**
   Used to ensure that topotoolbox is compiled and linked
   correctly. Always returns 1.
 */
TOPOTOOLBOX_API
int has_topotoolbox(void);

/**
   @brief Fills sinks in a digital elevation model

   Uses an algorithm based on grayscale morphological reconstruction.

   @param[out] output The filled DEM
   @param[in]  dem    The input DEM
   @param[in]  nrows  The size of both DEMs in the fastest changing dimension
   @param[in]  ncols  The size of both DEMs in the slowest changing dimension
 */
TOPOTOOLBOX_API
void fillsinks(float *output, float *dem, ptrdiff_t nrows, ptrdiff_t ncols);

/**
   @brief Labels flat, sill and presill pixels in the provided DEM

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

   @param[out] output The integer-valued output array with pixel labels
   @param[in]  dem    The input DEM
   @param[in]  nrows  The size of both DEMs in the fastest changing dimension
   @param[in]  ncols  The size of both DEMs in the slowest changing dimension
 */
TOPOTOOLBOX_API
ptrdiff_t identifyflats(int32_t *output, float *dem, ptrdiff_t nrows,
                        ptrdiff_t ncols);

/**
  @brief Compute costs for the gray-weighted distance transform

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

  @param[out] costs        The costs for the gray-weighted distance transform
  @param[out] conncomps    Labeled connected components
  @param[in]  flats        Array identifying the flat pixels
  @param[in]  original_dem The DEM prior to sink filling
  @param[in]  filled_dem   The DEM after sink filling
  @param[in]  nrows  The size of both DEMs in the fastest changing dimension
  @param[in]  ncols  The size of both DEMs in the slowest changing dimension
 */
TOPOTOOLBOX_API
void gwdt_computecosts(float *costs, ptrdiff_t *conncomps, int32_t *flats,
                       float *original_dem, float *filled_dem, ptrdiff_t nrows,
                       ptrdiff_t ncols);

/**
 @brief Compute excess topography with 2D varying threshold slopes
 using the fast sweeping method

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

 @param[out] excess           The solution of the constrained eikonal equation.
                              To compute the excess topography, subtract this
                              array elementwise from the DEM. A float array of
                              size (nrows x ncols).
 @param[in]  dem              The input digital elevation model. A float array
                              of size (nrows x ncols).
 @param[in]  threshold_slopes The threshold slopes (tangent of the critical
                              angle) at each grid cell. A float array of size
                              (nrows x ncols).
 @param[in]  cellsize         The spacing between grid cells, assumed to be
                              constant and identical in the x- and y- directions
 @param[in]  nrows            The size of the input and output DEMs and the
                              threshold_slopes array in the fastest changing
                              dimension
 @param[in]  ncols            The size of the input and output DEMs and the
                              threshold_slopes array in the slowest changing
                              dimension
 */
TOPOTOOLBOX_API
void excesstopography_fsm2d(float *excess, float *dem, float *threshold_slopes,
                            float cellsize, ptrdiff_t nrows, ptrdiff_t ncols);

/**
 @brief Compute excess topography with 2D varying threshold slopes
 using the fast marching method

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

 @param[out] excess           The solution of the constrained eikonal equation.
                              To compute the excess topography, subtract this
                              array elementwise from the DEM. A float array of
                              size (nrows x ncols).
 @param[in]  heap             A ptrdiff_t array of indices (nrows x ncols) used
                              for implementing the priority queue.
 @param[in]  back             A ptrdiff_t array of indices (nrows x ncols) used
                              for implementing the priority queue.
 @param[in]  dem              The input digital elevation model. A float array
                              of size (nrows x ncols).
 @param[in]  threshold_slopes The threshold slopes (tangent of the critical
                              angle) at each grid cell. A float array of size
                              (nrows x ncols).
 @param[in]  cellsize         The spacing between grid cells, assumed to be
                              constant and identical in the x- and y- directions
 @param[in]  nrows            The size of the input and output DEMs and the
                              threshold_slopes array in the fastest changing
                              dimension
 @param[in]  ncols            The size of the input and output DEMs and the
                              threshold_slopes array in the slowest changing
                              dimension
 */
TOPOTOOLBOX_API
void excesstopography_fmm2d(float *excess, ptrdiff_t *heap, ptrdiff_t *back,
                            float *dem, float *threshold_slopes, float cellsize,
                            ptrdiff_t nrows, ptrdiff_t ncols);

/**
 @brief Compute excess topography with three-dimensionally variable
 lithology using the fast marching method

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

 @param[out] excess           The solution of the constrained eikonal equation.
                              To compute the excess topography, subtract this
                              array elementwise from the DEM. A float array of
                              size (nrows x ncols).
 @param[in]  heap             A ptrdiff_t array of indices (nrows x ncols) used
                              for implementing the priority queue.
 @param[in]  back             A ptrdiff_t array of indices (nrows x ncols) used
                              for implementing the priority queue.
 @param[in]  dem              The input digital elevation model. A float array
                              of size (nrows x ncols).
 @param[in] lithstack         The input lithology. A three-dimensional float
                              array of size (nlayers x nrows x ncols). The value
                              of `lithstack[layer,row,col]` is the elevation of
                              the top surface of the given layer. Note that the
                              first dimension is the layer, so that the layers
                              of each cell are stored contiguously.
 @param[in]  threshold_slopes The threshold slopes (tangent of the critical
                              angle) for each layer. A float array of size
                              (nlayers).
 @param[in]  cellsize         The spacing between grid cells, assumed to be
                              constant and identical in the x- and y- directions
 @param[in]  nrows            The size of the input and output DEMs in the
                              fastest changing dimension.
 @param[in]  ncols            The size of the input and output DEMs and the
                              in the slowest changing
                              dimension.
 @param[in]  nlayers          The number of layers in the lithstack and
                              threshold_slopes arrays.
 */
TOPOTOOLBOX_API
void excesstopography_fmm3d(float *excess, ptrdiff_t *heap, ptrdiff_t *back,
                            float *dem, float *lithstack,
                            float *threshold_slopes, float cellsize,
                            ptrdiff_t nrows, ptrdiff_t ncols,
                            ptrdiff_t nlayers);

#endif  // TOPOTOOLBOX_H
