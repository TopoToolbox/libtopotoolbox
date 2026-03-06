#pragma once

#include <stddef.h>
#include <stdint.h>

/*
  Polyline geometry helpers — libtopotoolbox (internal)
  ======================================================

  Generic, swath-agnostic functions for manipulating rasterised and
  floating-point polylines.  Used by swaths.c and exposed as public API.

  Internal helpers (used by swaths.c):
    point_to_segment_distance  — signed perpendicular distance, point→segment
    bresenham_d8 / bresenham_d4 — integer rasterisation, D8 and D4

  Public API (also declared in topotoolbox.h):
    thin_rasterised_line_to_D8
    rasterize_path
    simplify_line
*/

// Signed perpendicular distance from point (px,py) to segment (ax,ay)→(bx,by).
// proj_x/proj_y: nearest point on segment; lambda: clamped parameter in [0,1].
// Returns positive left of a→b, negative right.
float point_to_segment_distance(float px, float py, float ax, float ay,
                                float bx, float by, float *proj_x,
                                float *proj_y, float *lambda);

// Bresenham D8 rasterisation: one pixel per step (diagonal or cardinal).
// skip_first: if nonzero, the start pixel is not written.
// Returns number of pixels written to out_i/out_j.
ptrdiff_t bresenham_d8(ptrdiff_t *out_i, ptrdiff_t *out_j, ptrdiff_t i0,
                       ptrdiff_t j0, ptrdiff_t i1, ptrdiff_t j1,
                       int skip_first);

// Bresenham D4 rasterisation: only cardinal moves (diagonal → two cardinals).
ptrdiff_t bresenham_d4(ptrdiff_t *out_i, ptrdiff_t *out_j, ptrdiff_t i0,
                       ptrdiff_t j0, ptrdiff_t i1, ptrdiff_t j1,
                       int skip_first);
