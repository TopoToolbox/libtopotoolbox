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

 The output is the solution of the constrained eikonal equation. To
 compute the excess topography, the caller must subtract the output from the
 original DEM.

 The threshold_slopes array should be the same size as the DEM and the
 output array and contain the threshold slope, the tangent of the
 critical angle.

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
 Himalaya-Karakoram ranges.  Geology, 43, 6, 523-526.

 Zhao, Hongkai (2004). A fast sweeping method for eikonal
 equations. Mathematics of Computation, 74, 250, 603-627.

 @param[out] excess           The topography constrained by the threshold slopes
 @param[in]  dem              The input digital elevation model
 @param[in]  threshold_slopes The threshold slopes at each grid cell
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

TOPOTOOLBOX_API
void excesstopography_fmm2d(float *excess, ptrdiff_t *heap, ptrdiff_t *back,
                            float *dem, float *threshold_slopes, float cellsize,
                            ptrdiff_t nrows, ptrdiff_t ncols);

TOPOTOOLBOX_API
void excesstopography_fmm3d(float *excess, ptrdiff_t *heap, ptrdiff_t *back,
                            float *dem, float *lithstack,
                            float *threshold_slopes, float cellsize,
                            ptrdiff_t nrows, ptrdiff_t ncols,
                            ptrdiff_t nlayers);

#endif  // TOPOTOOLBOX_H
