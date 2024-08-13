#pragma once

#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "graphflood/define_types.h"

// Define the element structure in the priority queue
typedef struct {
  GF_UINT key;        // The key associated with the element
  GF_FLOAT priority;  // The priority of the element (lower values have higher
                      // priority)
} PFElement;

// Priority Queue structure
typedef struct {
  PFElement* data;
  GF_UINT size;
  GF_UINT capacity;
} PFPQueue;

// Initialize the priority queue with a given capacity
static inline bool pfpq_init(PFPQueue* pq, GF_UINT capacity) {
  // This line tries to suppress warnings about potentially uninitialised
  // varibles in the struct
  memset(pq, 0, sizeof(*pq));
  pq->data = (PFElement*)malloc(capacity * sizeof(PFElement));
  if (pq->data == NULL) {
    return false;  // Allocation failed
  }
  pq->size = 0;
  pq->capacity = capacity;
  return true;
}

// Free the priority queue
static inline void pfpq_free(PFPQueue* pq) {
  free(pq->data);
  pq->data = NULL;
  pq->capacity = 0;
  pq->size = 0;
}

// Swap two elements in the queue
static inline void pfpq_swap(PFElement* a, PFElement* b) {
  PFElement temp = *a;
  *a = *b;
  *b = temp;
}

// Push an element into the priority queue
static inline bool pfpq_push(PFPQueue* pq, GF_UINT key, GF_FLOAT priority) {
  if (pq->size >= pq->capacity) {
    return false;  // Priority queue is full
  }

  // Insert the new element at the end
  GF_UINT index = pq->size++;
  pq->data[index].key = key;
  pq->data[index].priority = priority;

  // Heapify up
  while (index > 0) {
    GF_UINT parent = (index - 1) / 2;
    if (pq->data[index].priority >= pq->data[parent].priority) {
      break;
    }
    pfpq_swap(&pq->data[index], &pq->data[parent]);
    index = parent;
  }

  return true;
}

// Get the key of the top element without removing it
static inline GF_UINT pfpq_top_key(PFPQueue* pq) {
  return (pq->size > 0) ? pq->data[0].key : 0;
}

// Get the priority of the top element without removing it
static inline GF_FLOAT pfpq_top_priority(PFPQueue* pq) {
  return (pq->size > 0) ? pq->data[0].priority : 0.0f;
}

// Get the priority of the top element without removing it
static inline bool pfpq_empty(PFPQueue* pq) {
  return (pq->size > 0) ? false : true;
}

// Pop the top element from the priority queue and get its key
static inline bool pfpq_pop(PFPQueue* pq, GF_UINT* key, GF_FLOAT* priority) {
  if (pq->size == 0) {
    return false;  // Priority queue is empty
  }

  *key = pq->data[0].key;
  *priority = pq->data[0].priority;

  // Replace the root with the last element
  pq->data[0] = pq->data[--pq->size];

  // Heapify down
  GF_UINT index = 0;
  while (true) {
    GF_UINT left = 2 * index + 1;
    GF_UINT right = 2 * index + 2;
    GF_UINT smallest = index;

    if (left < pq->size &&
        pq->data[left].priority < pq->data[smallest].priority) {
      smallest = left;
    }
    if (right < pq->size &&
        pq->data[right].priority < pq->data[smallest].priority) {
      smallest = right;
    }
    if (smallest == index) {
      break;
    }

    pfpq_swap(&pq->data[index], &pq->data[smallest]);
    index = smallest;
  }

  return true;
}

// Pop the top element from the priority queue and get its key
static inline GF_UINT pfpq_pop_and_get_key(PFPQueue* pq) {
  if (pq->size == 0) {
    return false;  // Priority queue is empty
  }

  GF_UINT key = pq->data[0].key;
  // Replace the root with the last element
  pq->data[0] = pq->data[--pq->size];

  // Heapify down
  GF_UINT index = 0;
  while (true) {
    GF_UINT left = 2 * index + 1;
    GF_UINT right = 2 * index + 2;
    GF_UINT smallest = index;

    if (left < pq->size &&
        pq->data[left].priority < pq->data[smallest].priority) {
      smallest = left;
    }
    if (right < pq->size &&
        pq->data[right].priority < pq->data[smallest].priority) {
      smallest = right;
    }
    if (smallest == index) {
      break;
    }

    pfpq_swap(&pq->data[index], &pq->data[smallest]);
    index = smallest;
  }

  return key;
}
