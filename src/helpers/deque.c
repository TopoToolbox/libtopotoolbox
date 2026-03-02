#include "deque.h"

void sedq_init(SlidingExtremaDQ *dq, ptrdiff_t *idx_buf, float *val_buf) {
  dq->idx = idx_buf;
  dq->val = val_buf;
  dq->head = 0;
  dq->tail = -1;
}

void sedq_push_min(SlidingExtremaDQ *dq, ptrdiff_t i, float v) {
  while (dq->head <= dq->tail && dq->val[dq->tail] >= v) dq->tail--;
  dq->idx[++dq->tail] = i;
  dq->val[dq->tail] = v;
}

void sedq_push_max(SlidingExtremaDQ *dq, ptrdiff_t i, float v) {
  while (dq->head <= dq->tail && dq->val[dq->tail] <= v) dq->tail--;
  dq->idx[++dq->tail] = i;
  dq->val[dq->tail] = v;
}

void sedq_evict(SlidingExtremaDQ *dq, ptrdiff_t i) {
  if (dq->head <= dq->tail && dq->idx[dq->head] == i) dq->head++;
}

ptrdiff_t sedq_front_idx(const SlidingExtremaDQ *dq) {
  return dq->idx[dq->head];
}
float sedq_front_val(const SlidingExtremaDQ *dq) { return dq->val[dq->head]; }
int sedq_empty(const SlidingExtremaDQ *dq) { return dq->head > dq->tail; }
