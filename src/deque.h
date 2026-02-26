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
    float     *val;
    ptrdiff_t  head, tail;
} MonoDeque;

// Initialise with caller-provided arrays of length >= capacity.
void mdeque_init(MonoDeque *dq, ptrdiff_t *idx_buf, float *val_buf);

// Push (i, v): evict dominated tail entries, then append.
// Use push_min for a min-deque (front = minimum of window).
void mdeque_push_min(MonoDeque *dq, ptrdiff_t i, float v);
// Use push_max for a max-deque (front = maximum of window).
void mdeque_push_max(MonoDeque *dq, ptrdiff_t i, float v);

// Evict the front entry if its index equals i (call when block i leaves window).
void mdeque_evict(MonoDeque *dq, ptrdiff_t i);

// Front index / value.  Only call when !mdeque_empty().
ptrdiff_t mdeque_front_idx(const MonoDeque *dq);
float     mdeque_front_val(const MonoDeque *dq);

int mdeque_empty(const MonoDeque *dq);

#endif // DEQUE_H
