#include "deque.h"

void mdeque_init(MonoDeque *dq, ptrdiff_t *idx_buf, float *val_buf) {
    dq->idx  = idx_buf;
    dq->val  = val_buf;
    dq->head = 0;
    dq->tail = -1;
}

void mdeque_push_min(MonoDeque *dq, ptrdiff_t i, float v) {
    while (dq->head <= dq->tail && dq->val[dq->tail] >= v) dq->tail--;
    dq->idx[++dq->tail] = i;
    dq->val[  dq->tail] = v;
}

void mdeque_push_max(MonoDeque *dq, ptrdiff_t i, float v) {
    while (dq->head <= dq->tail && dq->val[dq->tail] <= v) dq->tail--;
    dq->idx[++dq->tail] = i;
    dq->val[  dq->tail] = v;
}

void mdeque_evict(MonoDeque *dq, ptrdiff_t i) {
    if (dq->head <= dq->tail && dq->idx[dq->head] == i) dq->head++;
}

ptrdiff_t mdeque_front_idx(const MonoDeque *dq) { return dq->idx[dq->head]; }
float     mdeque_front_val(const MonoDeque *dq) { return dq->val[dq->head]; }
int       mdeque_empty(const MonoDeque *dq)     { return dq->head > dq->tail; }
