#include <stddef.h>
#include <stdint.h>

/*
  PCG4D hash function


  Jarzynski, Mark and Olano, Marc. (2020). Hash functions for GPU
  rendering. Journal of Computer Graphics Techniques. Vol 9,
  No. 3. 21-38.
 */
float pcg4d(uint32_t a, uint32_t b, uint32_t c, uint32_t d);
