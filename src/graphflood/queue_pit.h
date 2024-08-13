/*
Fifo data structure I use in priority flood

Simple implementation adapted from ChatGPT's example


B.G.

*/
#pragma once

#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>

#include "graphflood/define_types.h"

// Define the queue structure
typedef struct {
  GF_UINT* buffer;
  int front;
  int rear;
  int size;
  int capacity;
} PitQueue;

// Function to initialize the queue with a given capacity
static inline bool pitqueue_init(PitQueue* q, int capacity) {
  memset(q, 0, sizeof(*q));
  q->buffer = (GF_UINT*)malloc(capacity * sizeof(GF_UINT));
  if (q->buffer == NULL) {
    return false;  // Allocation failed
  }
  q->front = 0;
  q->rear = -1;
  q->size = 0;
  q->capacity = capacity;
  return true;
}

// Function to check if the queue is empty
static inline bool pitqueue_is_empty(PitQueue* q) { return (q->size == 0); }

// Function to check if the queue is full
static inline bool pitqueue_is_full(PitQueue* q) {
  return (q->size == q->capacity);
}

// Function to enqueue an element into the queue
static inline bool pitqueue_enqueue(PitQueue* q, GF_UINT value) {
  if (pitqueue_is_full(q)) {
    return false;  // Queue is full
  }
  q->rear = (q->rear + 1) % q->capacity;
  q->buffer[q->rear] = value;
  q->size++;
  return true;
}

// Function to dequeue an element from the queue
static inline bool pitqueue_dequeue(PitQueue* q, GF_UINT* value) {
  if (pitqueue_is_empty(q)) {
    return false;  // Queue is empty
  }
  *value = q->buffer[q->front];
  q->front = (q->front + 1) % q->capacity;
  q->size--;
  return true;
}

// Function to dequeue an element from the queue
static inline GF_UINT pitqueue_pop_and_get(PitQueue* q) {
  GF_UINT value = q->buffer[q->front];
  q->front = (q->front + 1) % q->capacity;
  q->size--;
  return value;
}

// Function to free the allocated memory for the queue
static inline void pitqueue_free(PitQueue* q) {
  free(q->buffer);
  q->buffer = NULL;
  q->capacity = 0;
  q->size = 0;
  q->front = 0;
  q->rear = -1;
}
