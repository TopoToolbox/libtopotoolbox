#pragma once

#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "graphflood/define_types.h"

// Define the element structure in the max-heap priority queue
typedef struct {
  GF_UINT key;        // The key associated with the element
  GF_FLOAT priority;  // The priority of the element (higher values have higher
                      // priority)
  GF_FLOAT Qw;        // Transported discharge quantity
} MaxHeapElement;

// Max-Heap Priority Queue structure
typedef struct {
  MaxHeapElement* data;
  GF_UINT size;
  GF_UINT capacity;
} MaxHeapPQueue;

// Initialize the max-heap priority queue with a given capacity
static inline bool maxheap_init(MaxHeapPQueue* pq, GF_UINT capacity) {
  memset(pq, 0, sizeof(*pq));
  pq->data = (MaxHeapElement*)malloc(capacity * sizeof(MaxHeapElement));
  if (pq->data == NULL) {
    return false;  // Allocation failed
  }
  pq->size = 0;
  pq->capacity = capacity;
  return true;
}

// Free the max-heap priority queue
static inline void maxheap_free(MaxHeapPQueue* pq) {
  free(pq->data);
  pq->data = NULL;
  pq->capacity = 0;
  pq->size = 0;
}

// Swap two elements in the queue
static inline void maxheap_swap(MaxHeapElement* a, MaxHeapElement* b) {
  MaxHeapElement temp = *a;
  *a = *b;
  *b = temp;
}

// Push an element into the max-heap priority queue
static inline bool maxheap_push(MaxHeapPQueue* pq, GF_UINT key,
                                GF_FLOAT priority) {
  if (pq->size >= pq->capacity) {
    return false;  // Priority queue is full
  }

  // Insert the new element at the end
  GF_UINT index = pq->size++;
  pq->data[index].key = key;
  pq->data[index].priority = priority;
  pq->data[index].Qw = 0.0;

  // Heapify up (for max heap, parent should be >= children)
  while (index > 0) {
    GF_UINT parent = (index - 1) / 2;
    if (pq->data[index].priority <= pq->data[parent].priority) {
      break;
    }
    maxheap_swap(&pq->data[index], &pq->data[parent]);
    index = parent;
  }

  return true;
}

// Push an element with discharge into the max-heap priority queue
static inline bool maxheap_push_with_qw(MaxHeapPQueue* pq, GF_UINT key,
                                        GF_FLOAT priority, GF_FLOAT Qw) {
  if (pq->size >= pq->capacity) {
    return false;  // Priority queue is full
  }

  // Insert the new element at the end
  GF_UINT index = pq->size++;
  pq->data[index].key = key;
  pq->data[index].priority = priority;
  pq->data[index].Qw = Qw;

  // Heapify up (for max heap, parent should be >= children)
  while (index > 0) {
    GF_UINT parent = (index - 1) / 2;
    if (pq->data[index].priority <= pq->data[parent].priority) {
      break;
    }
    maxheap_swap(&pq->data[index], &pq->data[parent]);
    index = parent;
  }

  return true;
}

// Check if the max-heap is empty
static inline bool maxheap_empty(MaxHeapPQueue* pq) {
  return (pq->size > 0) ? false : true;
}

// Get the key of the top element without removing it
static inline GF_UINT maxheap_top_key(MaxHeapPQueue* pq) {
  return (pq->size > 0) ? pq->data[0].key : 0;
}

// Get the priority of the top element without removing it
static inline GF_FLOAT maxheap_top_priority(MaxHeapPQueue* pq) {
  return (pq->size > 0) ? pq->data[0].priority : 0.0f;
}

// Pop the top element from the max-heap priority queue and get its key
static inline GF_UINT maxheap_pop_and_get_key(MaxHeapPQueue* pq) {
  if (pq->size == 0) {
    return 0;  // Priority queue is empty
  }

  GF_UINT key = pq->data[0].key;
  // Replace the root with the last element
  pq->data[0] = pq->data[--pq->size];

  // Heapify down (for max heap, parent should be >= children)
  GF_UINT index = 0;
  while (true) {
    GF_UINT left = 2 * index + 1;
    GF_UINT right = 2 * index + 2;
    GF_UINT largest = index;

    if (left < pq->size &&
        pq->data[left].priority > pq->data[largest].priority) {
      largest = left;
    }
    if (right < pq->size &&
        pq->data[right].priority > pq->data[largest].priority) {
      largest = right;
    }
    if (largest == index) {
      break;
    }

    maxheap_swap(&pq->data[index], &pq->data[largest]);
    index = largest;
  }

  return key;
}

// Pop the top element from the max-heap priority queue and get its key and Qw
static inline GF_UINT maxheap_pop_and_get_key_qw(MaxHeapPQueue* pq,
                                                 GF_FLOAT* Qw_out) {
  if (pq->size == 0) {
    *Qw_out = 0.0;
    return 0;  // Priority queue is empty
  }

  GF_UINT key = pq->data[0].key;
  *Qw_out = pq->data[0].Qw;
  // Replace the root with the last element
  pq->data[0] = pq->data[--pq->size];

  // Heapify down (for max heap, parent should be >= children)
  GF_UINT index = 0;
  while (true) {
    GF_UINT left = 2 * index + 1;
    GF_UINT right = 2 * index + 2;
    GF_UINT largest = index;

    if (left < pq->size &&
        pq->data[left].priority > pq->data[largest].priority) {
      largest = left;
    }
    if (right < pq->size &&
        pq->data[right].priority > pq->data[largest].priority) {
      largest = right;
    }
    if (largest == index) {
      break;
    }

    maxheap_swap(&pq->data[index], &pq->data[largest]);
    index = largest;
  }

  return key;
}
