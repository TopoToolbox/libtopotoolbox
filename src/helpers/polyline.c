#define TOPOTOOLBOX_BUILD

#include "polyline.h"

#include <float.h>
#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "priority_queue.h"
#include "topotoolbox.h"

// 8-connected neighbor offsets, used by thin_rasterised_line_to_D8.
static const int k_di8[8] = {-1, -1, -1, 0, 0, 1, 1, 1};
static const int k_dj8[8] = {-1, 0, 1, -1, 1, -1, 0, 1};

// ============================================================================
// Signed perpendicular distance from a point to a segment (pixel space)
// ============================================================================
// Returns positive if the point is left of a→b, negative if right.
// proj_x/proj_y: closest point on segment; lambda: clamped parameter in [0,1].
float point_to_segment_distance(float px, float py, float ax, float ay,
                                float bx, float by, float *proj_x,
                                float *proj_y, float *lambda) {
  float dx = bx - ax;
  float dy = by - ay;
  float dpx = px - ax;
  float dpy = py - ay;
  float seg_length_sq = dx * dx + dy * dy;

  if (seg_length_sq < 1e-10f) {
    *lambda = 0.0f;
    *proj_x = ax;
    *proj_y = ay;
    return sqrtf(dpx * dpx + dpy * dpy);
  }

  float t = (dpx * dx + dpy * dy) / seg_length_sq;
  if (t < 0.0f)
    t = 0.0f;
  else if (t > 1.0f)
    t = 1.0f;

  *lambda = t;
  *proj_x = ax + t * dx;
  *proj_y = ay + t * dy;

  float to_p_x = px - *proj_x;
  float to_p_y = py - *proj_y;
  float dist = sqrtf(to_p_x * to_p_x + to_p_y * to_p_y);
  float cross = dpx * dy - dpy * dx;

  return cross >= 0.0f ? dist : -dist;
}

// ============================================================================
// Bresenham rasterization — D8 (diagonal allowed) and D4 (cardinal only)
// ============================================================================

ptrdiff_t bresenham_d8(ptrdiff_t *out_i, ptrdiff_t *out_j, ptrdiff_t i0,
                       ptrdiff_t j0, ptrdiff_t i1, ptrdiff_t j1,
                       int skip_first) {
  ptrdiff_t di = i1 - i0;
  ptrdiff_t dj = j1 - j0;
  ptrdiff_t si = di > 0 ? 1 : (di < 0 ? -1 : 0);
  ptrdiff_t sj = dj > 0 ? 1 : (dj < 0 ? -1 : 0);
  ptrdiff_t adi = di < 0 ? -di : di;
  ptrdiff_t adj = dj < 0 ? -dj : dj;

  ptrdiff_t count = 0;
  ptrdiff_t ci = i0, cj = j0;

  if (!skip_first) {
    out_i[count] = ci;
    out_j[count] = cj;
    count++;
  }

  if (adi >= adj) {
    ptrdiff_t err = adi / 2;
    for (ptrdiff_t step = 0; step < adi; step++) {
      err -= adj;
      if (err < 0) {
        cj += sj;
        err += adi;
      }
      ci += si;
      out_i[count] = ci;
      out_j[count] = cj;
      count++;
    }
  } else {
    ptrdiff_t err = adj / 2;
    for (ptrdiff_t step = 0; step < adj; step++) {
      err -= adi;
      if (err < 0) {
        ci += si;
        err += adj;
      }
      cj += sj;
      out_i[count] = ci;
      out_j[count] = cj;
      count++;
    }
  }

  return count;
}

ptrdiff_t bresenham_d4(ptrdiff_t *out_i, ptrdiff_t *out_j, ptrdiff_t i0,
                       ptrdiff_t j0, ptrdiff_t i1, ptrdiff_t j1,
                       int skip_first) {
  ptrdiff_t di = i1 - i0;
  ptrdiff_t dj = j1 - j0;
  ptrdiff_t si = di > 0 ? 1 : (di < 0 ? -1 : 0);
  ptrdiff_t sj = dj > 0 ? 1 : (dj < 0 ? -1 : 0);
  ptrdiff_t adi = di < 0 ? -di : di;
  ptrdiff_t adj = dj < 0 ? -dj : dj;

  ptrdiff_t count = 0;
  ptrdiff_t ci = i0, cj = j0;

  if (!skip_first) {
    out_i[count] = ci;
    out_j[count] = cj;
    count++;
  }

  if (adi >= adj) {
    ptrdiff_t err = adi / 2;
    for (ptrdiff_t step = 0; step < adi; step++) {
      err -= adj;
      if (err < 0) {
        cj += sj;
        out_i[count] = ci;
        out_j[count] = cj;
        count++;
        err += adi;
      }
      ci += si;
      out_i[count] = ci;
      out_j[count] = cj;
      count++;
    }
  } else {
    ptrdiff_t err = adj / 2;
    for (ptrdiff_t step = 0; step < adj; step++) {
      err -= adi;
      if (err < 0) {
        ci += si;
        out_i[count] = ci;
        out_j[count] = cj;
        count++;
        err += adj;
      }
      cj += sj;
      out_i[count] = ci;
      out_j[count] = cj;
      count++;
    }
  }

  return count;
}

// ============================================================================
// Elbow removal — thin a rasterised polyline to 1-pixel D8 width
// ============================================================================
TOPOTOOLBOX_API
ptrdiff_t thin_rasterised_line_to_D8(float *centre_line_i, float *centre_line_j,
                                     float *centre_width, ptrdiff_t n_centre,
                                     ptrdiff_t dims[2]) {
  if (n_centre <= 2) return n_centre;

  static const int c_di[4] = {-1, 0, 0, 1};
  static const int c_dj[4] = {0, -1, 1, 0};
  static const int d_di[4] = {-1, -1, 1, 1};
  static const int d_dj[4] = {-1, 1, -1, 1};

  ptrdiff_t total = dims[0] * dims[1];
  uint8_t *on_cl = (uint8_t *)calloc(total, 1);
  uint8_t *vis = (uint8_t *)calloc(total, 1);
  ptrdiff_t *path = (ptrdiff_t *)malloc(n_centre * sizeof(ptrdiff_t));
  if (!on_cl || !vis || !path) {
    free(on_cl);
    free(vis);
    free(path);
    return n_centre;
  }

  for (ptrdiff_t k = 0; k < n_centre; k++)
    on_cl[(ptrdiff_t)centre_line_j[k] * dims[0] + (ptrdiff_t)centre_line_i[k]] =
        1;

  int changed = 1;
  while (changed) {
    changed = 0;

    ptrdiff_t start = -1;
    for (ptrdiff_t k = 0; k < n_centre && start < 0; k++) {
      ptrdiff_t ci = (ptrdiff_t)centre_line_i[k];
      ptrdiff_t cj = (ptrdiff_t)centre_line_j[k];
      if (!on_cl[cj * dims[0] + ci]) continue;
      int nc = 0;
      for (int n = 0; n < 8; n++) {
        ptrdiff_t ni = ci + k_di8[n], nj = cj + k_dj8[n];
        if (ni < 0 || ni >= dims[0] || nj < 0 || nj >= dims[1]) continue;
        if (on_cl[nj * dims[0] + ni]) nc++;
      }
      if (nc == 1) start = cj * dims[0] + ci;
    }
    if (start < 0)
      start =
          (ptrdiff_t)centre_line_j[0] * dims[0] + (ptrdiff_t)centre_line_i[0];

    memset(vis, 0, (size_t)total);
    ptrdiff_t n_path = 0;
    ptrdiff_t cur = start;
    while (cur >= 0 && n_path < n_centre) {
      path[n_path++] = cur;
      vis[cur] = 1;
      ptrdiff_t ci = cur % dims[0], cj = cur / dims[0];
      ptrdiff_t nxt = -1;
      for (int n = 0; n < 4 && nxt < 0; n++) {
        ptrdiff_t ni = ci + c_di[n], nj = cj + c_dj[n];
        if (ni < 0 || ni >= dims[0] || nj < 0 || nj >= dims[1]) continue;
        ptrdiff_t nidx = nj * dims[0] + ni;
        if (on_cl[nidx] && !vis[nidx]) nxt = nidx;
      }
      for (int n = 0; n < 4 && nxt < 0; n++) {
        ptrdiff_t ni = ci + d_di[n], nj = cj + d_dj[n];
        if (ni < 0 || ni >= dims[0] || nj < 0 || nj >= dims[1]) continue;
        ptrdiff_t nidx = nj * dims[0] + ni;
        if (on_cl[nidx] && !vis[nidx]) nxt = nidx;
      }
      cur = nxt;
    }

    for (ptrdiff_t i = 1; i < n_path - 1; i++) {
      ptrdiff_t ai = path[i - 1] % dims[0], aj = path[i - 1] / dims[0];
      ptrdiff_t bi = path[i + 1] % dims[0], bj = path[i + 1] / dims[0];
      ptrdiff_t dii = ai - bi, djj = aj - bj;
      if (dii >= -1 && dii <= 1 && djj >= -1 && djj <= 1) {
        on_cl[path[i]] = 0;
        changed = 1;
        i++;
      }
    }
  }

  free(vis);
  free(path);

  ptrdiff_t new_n = 0;
  for (ptrdiff_t k = 0; k < n_centre; k++) {
    ptrdiff_t ci = (ptrdiff_t)centre_line_i[k];
    ptrdiff_t cj = (ptrdiff_t)centre_line_j[k];
    if (on_cl[cj * dims[0] + ci]) {
      centre_line_i[new_n] = centre_line_i[k];
      centre_line_j[new_n] = centre_line_j[k];
      centre_width[new_n] = centre_width[k];
      new_n++;
    }
  }
  free(on_cl);
  return new_n;
}

// ============================================================================
// Bresenham rasterization of a full polyline
// ============================================================================
TOPOTOOLBOX_API
ptrdiff_t rasterize_path(ptrdiff_t *out_i, ptrdiff_t *out_j,
                         const ptrdiff_t *ref_i, const ptrdiff_t *ref_j,
                         ptrdiff_t n_refs, int close_loop, int use_d4) {
  if (n_refs < 2) return 0;

  ptrdiff_t total = 0;
  ptrdiff_t n_segments = close_loop ? n_refs : n_refs - 1;

  for (ptrdiff_t k = 0; k < n_segments; k++) {
    ptrdiff_t i0 = ref_i[k];
    ptrdiff_t j0 = ref_j[k];
    ptrdiff_t i1 = ref_i[(k + 1) % n_refs];
    ptrdiff_t j1 = ref_j[(k + 1) % n_refs];

    int skip_first = (k > 0) ? 1 : 0;

    ptrdiff_t n;
    if (use_d4) {
      n = bresenham_d4(out_i + total, out_j + total, i0, j0, i1, j1,
                       skip_first);
    } else {
      n = bresenham_d8(out_i + total, out_j + total, i0, j0, i1, j1,
                       skip_first);
    }
    total += n;
  }

  return total;
}

// ============================================================================
// Line simplification — Iterative End-Point Fit (IEF) engine
// ============================================================================

static float seg_perp_dist(float ax, float ay, float bx, float by, float px,
                           float py) {
  float dx = bx - ax, dy = by - ay;
  float len2 = dx * dx + dy * dy;
  if (len2 < 1e-12f)
    return sqrtf((px - ax) * (px - ax) + (py - ay) * (py - ay));
  float cross = (px - ax) * dy - (py - ay) * dx;
  return fabsf(cross) / sqrtf(len2);
}

static float seg_tri_area(float ax, float ay, float bx, float by, float cx,
                          float cy) {
  return 0.5f * fabsf((bx - ax) * (cy - ay) - (by - ay) * (cx - ax));
}

typedef struct {
  float neg_dev;
  float seg_rss;
  ptrdiff_t lo, hi, max_k;
} ief_seg;

static ief_seg ief_make_seg(const float *in_i, const float *in_j, ptrdiff_t lo,
                            ptrdiff_t hi) {
  ief_seg s;
  s.lo = lo;
  s.hi = hi;
  s.max_k = lo;
  s.seg_rss = 0.0f;
  s.neg_dev = 0.0f;
  float max_d = -1.0f;
  for (ptrdiff_t k = lo + 1; k < hi; k++) {
    float d =
        seg_perp_dist(in_i[lo], in_j[lo], in_i[hi], in_j[hi], in_i[k], in_j[k]);
    s.seg_rss += d * d;
    if (d > max_d) {
      max_d = d;
      s.max_k = k;
    }
  }
  if (hi - lo > 1) s.neg_dev = -max_d;
  return s;
}

static ptrdiff_t ief_build(const float *in_i, const float *in_j, ptrdiff_t n,
                           ptrdiff_t *seq, float *rmse_curve) {
  if (n <= 2) return 0;

  ptrdiff_t max_slots = 2 * n;
  ief_seg *segs = (ief_seg *)malloc((size_t)max_slots * sizeof(ief_seg));
  ptrdiff_t *pq_heap =
      (ptrdiff_t *)malloc((size_t)max_slots * sizeof(ptrdiff_t));
  ptrdiff_t *pq_back =
      (ptrdiff_t *)malloc((size_t)max_slots * sizeof(ptrdiff_t));
  float *neg_devs = (float *)malloc((size_t)max_slots * sizeof(float));
  for (ptrdiff_t i = 0; i < max_slots; i++) pq_back[i] = -1;
  PriorityQueue q = pq_create(max_slots, pq_heap, pq_back, neg_devs, 0);
  ptrdiff_t next_slot = 0;

  ief_seg init = ief_make_seg(in_i, in_j, 0, n - 1);
  double total_rss = (double)init.seg_rss;
  segs[next_slot] = init;
  pq_insert(&q, next_slot, init.neg_dev);
  next_slot++;

  rmse_curve[0] = (float)sqrt(total_rss / (double)n);

  ptrdiff_t n_steps = n - 2;
  for (ptrdiff_t t = 0; t < n_steps; t++) {
    if (pq_isempty(&q)) {
      seq[t] = 0;
      rmse_curve[t + 1] = 0.0f;
      continue;
    }
    ptrdiff_t slot = pq_deletemin(&q);
    ief_seg top = segs[slot];
    seq[t] = top.max_k;
    total_rss -= (double)top.seg_rss;

    ief_seg left = ief_make_seg(in_i, in_j, top.lo, top.max_k);
    ief_seg right = ief_make_seg(in_i, in_j, top.max_k, top.hi);
    total_rss += (double)left.seg_rss + (double)right.seg_rss;
    if (total_rss < 0.0) total_rss = 0.0;

    if (left.hi - left.lo > 1) {
      segs[next_slot] = left;
      pq_insert(&q, next_slot, left.neg_dev);
      next_slot++;
    }
    if (right.hi - right.lo > 1) {
      segs[next_slot] = right;
      pq_insert(&q, next_slot, right.neg_dev);
      next_slot++;
    }

    rmse_curve[t + 1] = (float)sqrt(total_rss / (double)n);
  }
  free(segs);
  free(pq_heap);
  free(pq_back);
  free(neg_devs);
  return n_steps;
}

static ptrdiff_t simplify_fixed_n(float *out_i, float *out_j,
                                  const float *track_i, const float *track_j,
                                  ptrdiff_t n, float tolerance) {
  ptrdiff_t *seq = (ptrdiff_t *)malloc((size_t)(n - 2) * sizeof(ptrdiff_t));
  float *rmse_curve = (float *)malloc((size_t)(n - 1) * sizeof(float));
  ief_build(track_i, track_j, n, seq, rmse_curve);

  ptrdiff_t n_target = (ptrdiff_t)tolerance;
  if (n_target < 2) n_target = 2;
  if (n_target > n) n_target = n;

  uint8_t *keep = (uint8_t *)calloc((size_t)n, 1);
  keep[0] = 1;
  keep[n - 1] = 1;
  for (ptrdiff_t k = 0; k < n_target - 2; k++) keep[seq[k]] = 1;

  ptrdiff_t cnt = 0;
  for (ptrdiff_t k = 0; k < n; k++)
    if (keep[k]) {
      out_i[cnt] = track_i[k];
      out_j[cnt] = track_j[k];
      cnt++;
    }

  free(seq);
  free(rmse_curve);
  free(keep);
  return cnt;
}

static ptrdiff_t simplify_kneedle(float *out_i, float *out_j,
                                  const float *track_i, const float *track_j,
                                  ptrdiff_t n) {
  ptrdiff_t *seq = (ptrdiff_t *)malloc((size_t)(n - 2) * sizeof(ptrdiff_t));
  float *rmse_curve = (float *)malloc((size_t)(n - 1) * sizeof(float));
  ief_build(track_i, track_j, n, seq, rmse_curve);

  float y0 = rmse_curve[0];
  ptrdiff_t best_k = 0;
  float best_d = FLT_MAX;
  for (ptrdiff_t k = 0; k < n - 1; k++) {
    float xk = (n > 2) ? (float)k / (float)(n - 2) : 0.0f;
    float yk = (y0 > 1e-12f) ? rmse_curve[k] / y0 : 0.0f;
    float d = yk + xk - 1.0f;
    if (d < best_d) {
      best_d = d;
      best_k = k;
    }
  }
  ptrdiff_t n_target = best_k + 2;
  if (n_target < 2) n_target = 2;
  if (n_target > n) n_target = n;

  uint8_t *keep = (uint8_t *)calloc((size_t)n, 1);
  keep[0] = 1;
  keep[n - 1] = 1;
  for (ptrdiff_t k = 0; k < n_target - 2; k++) keep[seq[k]] = 1;

  ptrdiff_t cnt = 0;
  for (ptrdiff_t k = 0; k < n; k++)
    if (keep[k]) {
      out_i[cnt] = track_i[k];
      out_j[cnt] = track_j[k];
      cnt++;
    }

  free(seq);
  free(rmse_curve);
  free(keep);
  return cnt;
}

static ptrdiff_t simplify_vw_area(float *out_i, float *out_j,
                                  const float *track_i, const float *track_j,
                                  ptrdiff_t n, float tolerance) {
  ptrdiff_t *prv_kept = (ptrdiff_t *)malloc((size_t)n * sizeof(ptrdiff_t));
  ptrdiff_t *nxt_kept = (ptrdiff_t *)malloc((size_t)n * sizeof(ptrdiff_t));
  uint8_t *inserted = (uint8_t *)calloc((size_t)n, sizeof(uint8_t));
  ptrdiff_t *pq_heap = (ptrdiff_t *)malloc((size_t)n * sizeof(ptrdiff_t));
  ptrdiff_t *pq_back = (ptrdiff_t *)malloc((size_t)n * sizeof(ptrdiff_t));
  float *neg_area = (float *)malloc((size_t)n * sizeof(float));
  for (ptrdiff_t i = 0; i < n; i++) pq_back[i] = -1;
  PriorityQueue pq = pq_create(n, pq_heap, pq_back, neg_area, 0);

  inserted[0] = 1;
  inserted[n - 1] = 1;
  for (ptrdiff_t k = 1; k < n - 1; k++) {
    prv_kept[k] = 0;
    nxt_kept[k] = n - 1;
    float a = seg_tri_area(track_i[0], track_j[0], track_i[k], track_j[k],
                           track_i[n - 1], track_j[n - 1]);
    pq_insert(&pq, k, -a);
  }

  while (!pq_isempty(&pq)) {
    ptrdiff_t k = pq_deletemin(&pq);
    float area = -pq.priorities[k];
    if (area < tolerance) break;

    inserted[k] = 1;
    ptrdiff_t p = prv_kept[k], nx = nxt_kept[k];

    for (ptrdiff_t qi = p + 1; qi < k; qi++) {
      if (!inserted[qi]) {
        nxt_kept[qi] = k;
        float new_a =
            seg_tri_area(track_i[prv_kept[qi]], track_j[prv_kept[qi]],
                         track_i[qi], track_j[qi], track_i[k], track_j[k]);
        float new_prio = -new_a;
        if (new_prio < pq.priorities[qi])
          pq_decrease_key(&pq, qi, new_prio);
        else if (new_prio > pq.priorities[qi])
          pq_increase_key(&pq, qi, new_prio);
      }
    }
    for (ptrdiff_t qi = k + 1; qi < nx; qi++) {
      if (!inserted[qi]) {
        prv_kept[qi] = k;
        float new_a =
            seg_tri_area(track_i[k], track_j[k], track_i[qi], track_j[qi],
                         track_i[nxt_kept[qi]], track_j[nxt_kept[qi]]);
        float new_prio = -new_a;
        if (new_prio < pq.priorities[qi])
          pq_decrease_key(&pq, qi, new_prio);
        else if (new_prio > pq.priorities[qi])
          pq_increase_key(&pq, qi, new_prio);
      }
    }
  }

  ptrdiff_t cnt = 0;
  for (ptrdiff_t k = 0; k < n; k++)
    if (inserted[k]) {
      out_i[cnt] = track_i[k];
      out_j[cnt] = track_j[k];
      cnt++;
    }
  free(prv_kept);
  free(nxt_kept);
  free(inserted);
  free(pq_heap);
  free(pq_back);
  free(neg_area);
  return cnt;
}

ptrdiff_t simplify_line(float *out_i, float *out_j, const float *track_i,
                        const float *track_j, ptrdiff_t n_points,
                        float tolerance, int method) {
  if (n_points <= 2) {
    for (ptrdiff_t k = 0; k < n_points; k++) {
      out_i[k] = track_i[k];
      out_j[k] = track_j[k];
    }
    return n_points;
  }

  ptrdiff_t n = n_points;

  if (method == 0)
    return simplify_fixed_n(out_i, out_j, track_i, track_j, n, tolerance);
  if (method == 1) return simplify_kneedle(out_i, out_j, track_i, track_j, n);
  if (method == 2)
    return simplify_vw_area(out_i, out_j, track_i, track_j, n, tolerance);
  return 0;
}
