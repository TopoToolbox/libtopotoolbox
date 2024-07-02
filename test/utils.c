#include "utils.h"

#include <stddef.h>
#include <stdint.h>

float pcg4d(uint32_t a, uint32_t b, uint32_t c, uint32_t d) {
  uint32_t x = a * 1664525u + 1013904223u;
  uint32_t y = b * 1664525u + 1013904223u;
  uint32_t z = c * 1664525u + 1013904223u;
  uint32_t w = d * 1664525u + 1013904223u;

  x += y * w;
  y += z * x;
  z += x * y;
  w += y * z;

  x ^= x >> 16;
  y ^= y >> 16;
  z ^= z >> 16;
  w ^= w >> 16;

  x += y * w;
  y += z * x;
  z += x * y;
  w += y * z;

  return (float)(w >> 8) / (1 << 24);
}
