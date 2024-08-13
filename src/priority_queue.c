#include "priority_queue.h"

#include <stddef.h>
#include <stdint.h>

float pq_get_priority(PriorityQueue *q, ptrdiff_t key) {
  return q->priorities[key];
}

int32_t pq_isempty(PriorityQueue *q) { return q->count == 0; }

static ptrdiff_t isleaf(PriorityQueue *q, ptrdiff_t position) {
  return (q->count / 2 <= position) && (position < q->count);
}
static ptrdiff_t left_child(ptrdiff_t position) { return 2 * position + 1; }
static ptrdiff_t right_child(ptrdiff_t position) { return 2 * position + 2; }
static ptrdiff_t parent(ptrdiff_t position) { return (position - 1) / 2; }

static void swap(PriorityQueue *q, ptrdiff_t node1, ptrdiff_t node2) {
  ptrdiff_t a = q->heap[node1];
  ptrdiff_t b = q->heap[node2];

  // Update the back pointers
  q->back[a] = node2;
  q->back[b] = node1;

  // Update the heap
  q->heap[node1] = b;
  q->heap[node2] = a;
}

static void siftup(PriorityQueue *q, ptrdiff_t position) {
  while (position > 0) {
    ptrdiff_t node = parent(position);
    if (q->priorities[q->heap[node]] < q->priorities[q->heap[position]]) {
      return;
    }
    swap(q, position, node);
    position = node;
  }
}

static void siftdown(PriorityQueue *q, ptrdiff_t position) {
  while (!isleaf(q, position)) {
    ptrdiff_t child = left_child(position);

    // Which child has the smallest value?
    if ((child + 1 < q->count) &&
        (q->priorities[q->heap[child + 1]] < q->priorities[q->heap[child]])) {
      child = child + 1;
    }

    if (q->priorities[q->heap[child]] > q->priorities[q->heap[position]]) {
      return;
    }
    swap(q, position, child);
    position = child;
  }
}

static void build_heap(PriorityQueue *q) {
  for (ptrdiff_t i = parent(q->count - 1); i >= 0; i--) {
    siftdown(q, i);
  }
}

PriorityQueue pq_create(ptrdiff_t max_size, ptrdiff_t *heap, ptrdiff_t *back,
                        float *priorities, int32_t flag) {
  PriorityQueue q = {0};
  q.back = back;
  q.heap = heap;
  q.priorities = priorities;
  q.max_size = max_size;

  if (flag) {
    q.count = max_size;
    build_heap(&q);
  }

  return q;
}

void pq_insert(PriorityQueue *q, ptrdiff_t key, float priority) {
  q->back[key] = q->count;
  q->heap[q->count] = key;
  q->priorities[key] = priority;
  siftup(q, q->count++);
}

ptrdiff_t pq_deletemin(PriorityQueue *q) {
  ptrdiff_t root = q->heap[0];
  q->count--;
  if (q->count > 0) {
    // Otherwise we just deleted the root
    swap(q, 0, q->count);
    siftdown(q, 0);
  }
  q->back[root] = -1;  // Sentinel value for deleted values
  return root;
}

void pq_decrease_key(PriorityQueue *q, ptrdiff_t idx, float new_priority) {
  if (new_priority < q->priorities[idx]) {
    ptrdiff_t node = q->back[idx];
    q->priorities[idx] = new_priority;
    siftup(q, node);
  }
}

int32_t pq_hasminheap(PriorityQueue *q) {
  // Ensure that the min-heap property is satisfied by the priority queue
  int32_t flag = 1;
  for (ptrdiff_t i = 0; i < q->count / 2; i++) {
    if (left_child(i) < q->count &&
        q->priorities[q->heap[left_child(i)]] <= q->priorities[q->heap[i]]) {
      flag &= 0;
    }
    if (right_child(i) < q->count &&
        q->priorities[q->heap[right_child(i)]] <= q->priorities[q->heap[i]]) {
      flag &= 0;
    }
  }
  return flag;
}
