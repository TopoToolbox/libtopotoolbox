#ifndef TOPOTOOLBOX_H
#define TOPOTOOLBOX_H

#if defined(TOPOTOOLBOX_BUILD)
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

/*
  has_topotoolbox()

  Returns 1. Used to ensure that topotoolbox is compiled and linked
  correctly.
 */
TOPOTOOLBOX_API
int has_topotoolbox(void);

TOPOTOOLBOX_API
void fillsinks(float *output, float *dem, ptrdiff_t nrows, ptrdiff_t ncols);

#endif  // TOPOTOOLBOX_H
