#ifndef TOPOTOOLBOX_PRIORITY_QUEUE_H
#define TOPOTOOLBOX_PRIORITY_QUEUE_H

#include <stddef.h>
#include <stdint.h>

/*
  Priority queue implemented as a binary min-heap

  Priorities (i.e. elevations) are stored in their correct positions
  in a float array while indices are stored in a ptrdiff_t array.
 */
typedef struct {
  float *priorities;
  ptrdiff_t *back;  // Back pointers: heap[back[idx]] == idx
  ptrdiff_t *heap;
  ptrdiff_t max_size;  // Not sure if we really need this
  ptrdiff_t count;
} PriorityQueue;

// Create a priority queue
//
// If flag is nonzero, assume that the arrays have been pre-filled and
// run the heap algorithm to ensure that the queue is ready to use.
PriorityQueue pq_create(ptrdiff_t count, ptrdiff_t *heap, ptrdiff_t *back,
                        float *priorities, int32_t flag);

// Returns 1 if queue is empty, 0 if queue has elements
int32_t pq_isempty(PriorityQueue *q);

float pq_get_priority(PriorityQueue *q, ptrdiff_t key);

// Insert key with priority
void pq_insert(PriorityQueue *q, ptrdiff_t key, float priority);

// Delete the minimum key from the queue and return it
//
// The corresponding priority can still be retrieved with get_priority
ptrdiff_t pq_deletemin(PriorityQueue *q);

void pq_decrease_key(PriorityQueue *q, ptrdiff_t idx, float new_priority);
void pq_increase_key(PriorityQueue *q, ptrdiff_t idx, float new_priority);

// For testing
int32_t pq_hasminheap(PriorityQueue *q);

#endif  // TOPOTOOLBOX_PRIORITY_QUEUE_H
