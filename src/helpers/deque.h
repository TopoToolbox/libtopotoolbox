#ifndef DEQUE_H
#define DEQUE_H

#include <stddef.h>

// Monotonic deque for O(1) amortised sliding-window min/max.
//
// Stores (index, value) pairs.  head and tail are plain cursors that only
// ever advance right, so the backing arrays need capacity == window_size.
// The caller provides the backing storage (no internal allocation).
//
// Invariant: values are monotonically increasing (min) or decreasing (max)
// from front to back.  The front is always the extremum of the current window.

typedef struct {
  ptrdiff_t *idx;
  float *val;
  ptrdiff_t head, tail;
} SlidingExtremaDQ;

// Initialise with caller-provided arrays of length >= capacity.
void sedq_init(SlidingExtremaDQ *dq, ptrdiff_t *idx_buf, float *val_buf);

// Push (i, v): evict dominated tail entries, then append.
// Use sedq_push_min for a min-deque (front = minimum of window).
void sedq_push_min(SlidingExtremaDQ *dq, ptrdiff_t i, float v);
// Use sedq_push_max for a max-deque (front = maximum of window).
void sedq_push_max(SlidingExtremaDQ *dq, ptrdiff_t i, float v);

// Evict the front entry if its index equals i (call when block i leaves
// window).
void sedq_evict(SlidingExtremaDQ *dq, ptrdiff_t i);

// Front index / value.  Only call when !sedq_empty().
ptrdiff_t sedq_front_idx(const SlidingExtremaDQ *dq);
float sedq_front_val(const SlidingExtremaDQ *dq);

int sedq_empty(const SlidingExtremaDQ *dq);

#endif  // DEQUE_H
