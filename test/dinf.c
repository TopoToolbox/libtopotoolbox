#undef NDEBUG

#include <assert.h>
#include <gdal.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "flow.h"
#include "grid.h"
#include "topotoolbox.h"
#include "utils.h"

#define test(desc)                                               \
  for (int _res = 0, _break = 1; _break; _break = 0,             \
           _res ? printf("not ok %d - %s\n", ++_testcount, desc) \
                : printf("ok %d - %s\n", ++_testcount, desc),    \
           _success |= _res)
#define test_fail _res = 1
#define test_init()                                                    \
  for (int _success = 0, _testcount = (printf("TAP version 14\n"), 0), \
           _break = 1;                                                 \
       _break; _break = 0, printf("1..%d\n", _testcount), exit(_success))

ptrdiff_t ei[8] = {0, -1, -1, -1, 0, 1, 1, 1};
ptrdiff_t ej[8] = {1, 1, 0, -1, -1, -1, 0, 1};

void print_flow(Flow f) {
  printf("%ld edges\n", f.count);
  for (ptrdiff_t e = 0; e < f.count; e++) {
    printf("%ld %ld %f\n", f.source[e], f.target[e], f.fraction[e]);
  }
}

int test_node_once(Flow f) {
  // Has every node been assigned exactly once?
  Grid count = GridCreate(GridF32, NULL, f.cellsize, f.dims);
  for (ptrdiff_t j = 0; j < f.dims[1]; j++) {
    for (ptrdiff_t i = 0; i < f.dims[0]; i++) {
      ((float *)count.data)[f.stream[j * f.dims[0] + i]] += 1.0;
    }
  }

  Grid ones = GridFill(1.0, count.cellsize, count.dims);

  int res = GridEq(count, ones);

  GridFree(&count);
  GridFree(&ones);

  return res;
}

int test_tsort_node(Flow fd) {
  // Is the node list topologically sorted?
  // For every edge u => v, u comes before v in the node list

  int res = 0;

  // idxs maps every cell to its position in the node list
  Grid idxs = GridCreate(GridIdx, NULL, fd.cellsize, fd.dims);
  for (ptrdiff_t j = 0; j < fd.dims[1]; j++) {
    for (ptrdiff_t i = 0; i < fd.dims[0]; i++) {
      ((ptrdiff_t *)idxs.data)[fd.stream[j * fd.dims[0] + i]] =
          j * fd.dims[0] + i;
    }
  }

  for (ptrdiff_t e = 0; e < fd.count; e++) {
    ptrdiff_t u = fd.source[e];
    ptrdiff_t v = fd.target[e];

    res = ((ptrdiff_t *)idxs.data)[u] < ((ptrdiff_t *)idxs.data)[v] ? 1 : 0;
  }

  GridFree(&idxs);
  return res;
}

int test_tsort(Flow fd) {
  // Are the edge lists topologically sorted?
  // All incoming edges to a vertex v occur before any outgoing edges

  int tsort_success = 1;

  Grid outdegree = GridCreate(GridU8, NULL, fd.cellsize, fd.dims);

  for (ptrdiff_t e = 0; e < fd.count; e++) {
    ptrdiff_t u = fd.source[e];
    ptrdiff_t v = fd.target[e];

    ((uint8_t *)outdegree.data)[u] += 1;
    if (((uint8_t *)outdegree.data)[v] > 0) {
      tsort_success = 0;
    }
  }
  if (!tsort_success) {
    printf("# Edge list was not topologically sorted\n");
  }
  GridFree(&outdegree);
  return tsort_success;
}

int test_nanflow(Flow fd, Grid dem) {
  // Missing data in the DEM should have no flow into or out of them.
  int res = 1;

  for (ptrdiff_t e = 0; e < fd.count; e++) {
    ptrdiff_t u = fd.source[e];
    ptrdiff_t v = fd.target[e];
    if (isnan(((float *)dem.data)[u])) {
      res = 0;
      printf("# NaN pixel of DEM is a source: %ld\n", u);
    }

    if (isnan(((float *)dem.data)[v])) {
      res = 0;
      printf("# NaN pixel of DEM is a target: %ld\n", v);
    }
  }

  return res;
}

int test_no_boundary_edges(Flow fd) {
  // No edges should originate from outside the domain or leave the
  // domain
  int res = 1;
  for (ptrdiff_t e = 0; e < fd.count; e++) {
    ptrdiff_t u = fd.source[e];
    ptrdiff_t v = fd.target[e];

    ptrdiff_t ui = u % fd.dims[0];
    ptrdiff_t uj = u / fd.dims[0];

    ptrdiff_t vi = v % fd.dims[0];
    ptrdiff_t vj = v / fd.dims[0];

    if (llabs(ui - vi) > 1 || llabs(uj - vj) > 1) {
      printf("# Indices (%ld, %ld) and (%ld,%ld) are not neighbors.\n", ui, uj,
             vi, vj);
    }

    if (ui < 0 || ui >= fd.dims[0] || uj < 0 || uj >= fd.dims[1]) {
      printf("# Edge %ld => %ld originates from outside the domain.\n", u, v);
      res = 0;
    }

    if (vi < 0 || vi >= fd.dims[0] || vj < 0 || vj >= fd.dims[1]) {
      printf("# Edge %ld => %ld leaves the domain.\n", u, v);
      res = 0;
    }
  }
  return res;
}

int test_outlets(Flow fd) {
  // There should be at least one outlet, a pixel with incoming edges
  // but no outgoing edges.
  Grid outdegree = GridCreate(GridU8, NULL, fd.cellsize, fd.dims);
  Grid indegree = GridCreate(GridU8, NULL, fd.cellsize, fd.dims);

  uint8_t *outdegree_data = outdegree.data;
  uint8_t *indegree_data = indegree.data;

  for (ptrdiff_t e = 0; e < fd.count; e++) {
    ptrdiff_t u = fd.source[e];
    ptrdiff_t v = fd.target[e];

    outdegree_data[u] += 1;
    indegree_data[v] += 1;
  }

  int res = 0;
  for (ptrdiff_t j = 0; j < fd.dims[1]; j++) {
    for (ptrdiff_t i = 0; i < fd.dims[0]; i++) {
      if (indegree_data[j * fd.dims[0] + i] > 0 &&
          outdegree_data[j * fd.dims[0] + i] == 0) {
        // This is an outlet
        res = 1;
      }
    }
  }

  GridFree(&outdegree);
  GridFree(&indegree);

  return res;
}

int test_sum_edge_weights(Flow fd) {
  // Outgoing edge weights should sum to one
  int res = 1;

  Grid edge_sums = GridCreate(GridF32, NULL, fd.cellsize, fd.dims);
  for (ptrdiff_t e = 0; e < fd.count; e++) {
    ptrdiff_t u = fd.source[e];

    ((float *)edge_sums.data)[u] += fd.fraction[e];
  }

  for (ptrdiff_t e = 0; e < fd.count; e++) {
    ptrdiff_t u = fd.source[e];

    if (((float *)edge_sums.data)[u] != 1.0) {
      printf("# Outgoing edges from pixel %ld do not sum to 1.\n", u);
      res = 0;
    }
  }
  GridFree(&edge_sums);

  return res;
}

Flow dinf_flow(Grid dem, int order) {
  Flow fd = {0};

  Grid direction = GridCreate(GridU8, NULL, dem.cellsize, dem.dims);
  Grid prop = GridCreate(GridF32, NULL, dem.cellsize, dem.dims);

  flow_routing_dinf_directions(direction.data, prop.data, dem.data, dem.dims,
                               order);

  ptrdiff_t edge_count = edgeset_count(direction.data, dem.dims);

  float *weight = calloc(edge_count, sizeof *weight);

  flow_routing_dinf_weights(weight, direction.data, prop.data, direction.dims);

  // Resolve flats
  Grid rdirection = GridCreate(GridU8, NULL, dem.cellsize, dem.dims);
  Grid resolved = GridCreate(GridU8, NULL, dem.cellsize, dem.dims);

  // Initialize resolved according to the flow directions
  for (ptrdiff_t j = 0; j < dem.dims[1]; j++) {
    for (ptrdiff_t i = 0; i < dem.dims[0]; i++) {
      if (isnan(((float *)dem.data)[j * dem.dims[0] + i])) {
        // If DEM is a NaN, treat the flow direction as resolved
        ((uint8_t *)resolved.data)[j * dem.dims[0] + i] = 0xff;
      } else {
        ((uint8_t *)resolved.data)[j * dem.dims[0] + i] =
            ((uint8_t *)direction.data)[j * dem.dims[0] + i] ? 0xff : 0;
      }
    }
  }

  Grid path_distance = GridCreate(GridF32, NULL, dem.cellsize, dem.dims);
  Grid heap = GridCreate(GridIdx, NULL, dem.cellsize, dem.dims);
  Grid back = GridCreate(GridIdx, NULL, dem.cellsize, dem.dims);

  resolve_flats_shortest_path(rdirection.data, resolved.data,
                              path_distance.data, heap.data, back.data,
                              dem.data, dem.dims, order);
  GridFree(&resolved);
  GridFree(&path_distance);
  GridFree(&heap);
  GridFree(&back);

  // Compute weights for shortest path resolution
  ptrdiff_t rcount = edgeset_count(rdirection.data, rdirection.dims);
  float *rweight = calloc(rcount, sizeof *rweight);

  resolve_flats_shortest_path_weights(rweight, rcount);

  // Merge edge sets
  edge_count =
      edgeset_count_merged(direction.data, rdirection.data, direction.dims);

  float *merged_weights = calloc(edge_count, sizeof *merged_weights);
  Grid mergedscan = GridCreate(GridIdx, NULL, dem.cellsize, dem.dims);

  ptrdiff_t edge_count2 =
      edgeset_merge(merged_weights, mergedscan.data, direction.data, weight,
                    rdirection.data, rweight, direction.dims);

  if (edge_count != edge_count2) {
    printf("# Merged edge counts not equal: %ld != %ld", edge_count,
           edge_count2);
  }

  // Compute the topological sort of the merged edge set.
  Grid stream = GridCreate(GridIdx, NULL, dem.cellsize, dem.dims);
  ptrdiff_t *source = calloc(edge_count, sizeof *source);
  ptrdiff_t *target = calloc(edge_count, sizeof *target);
  float *sorted_weight = calloc(edge_count, sizeof *sorted_weight);

  Grid stack = GridCreate(GridIdx, NULL, dem.cellsize, dem.dims);
  Grid stackdir = GridCreate(GridU8, NULL, dem.cellsize, dem.dims);
  Grid visited = GridCreate(GridU8, NULL, dem.cellsize, dem.dims);

  flow_routing_tsort(stream.data, source, target, sorted_weight, stack.data,
                     stackdir.data, direction.data, merged_weights,
                     mergedscan.data, visited.data, edge_count, dem.dims,
                     order);

  fd.count = edge_count;
  fd.stream = stream.data;
  fd.source = source;
  fd.target = target;
  fd.fraction = sorted_weight;
  fd.cellsize = dem.cellsize;
  fd.dims[0] = dem.dims[0];
  fd.dims[1] = dem.dims[1];

  GridFree(&stack);
  GridFree(&stackdir);
  GridFree(&visited);
  GridFree(&mergedscan);
  GridFree(&direction);
  GridFree(&prop);
  free(weight);
  free(merged_weights);
  free(rweight);
  GridFree(&rdirection);

  return fd;
}

int main(int argc, char *argv[]) {
  GDALAllRegister();

  const char *snapshot_directory = "";

  if (argc >= 2) {
    snapshot_directory = argv[1];
  }

  test_init() {
    test("Slope") {
      ptrdiff_t dims[2] = {3, 3};
      float data[9] = {2, 2, 2, 1, 1, 1, 0, 0, 0};
      Grid dem = GridCreate(GridF32, data, 1.0, dims);
      Flow fd = dinf_flow(dem, 0);

      if (!test_node_once(fd)) test_fail;
      if (!test_tsort_node(fd)) test_fail;
      if (!test_tsort(fd)) test_fail;
      if (!test_nanflow(fd, dem)) test_fail;
      if (!test_no_boundary_edges(fd)) test_fail;
      if (!test_outlets(fd)) test_fail;
      if (!test_sum_edge_weights(fd)) test_fail;

      FlowFree(&fd);
    }

    test("Larger slope") {
      ptrdiff_t dims[2] = {5, 5};
      float data[25] = {4, 4, 4, 4, 4, 3, 3, 3, 3, 3, 2, 2, 2,
                        2, 2, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0};
      Grid dem = GridCreate(GridF32, data, 1.0, dims);
      Flow fd = dinf_flow(dem, 0);

      if (!test_node_once(fd)) test_fail;
      if (!test_tsort_node(fd)) test_fail;
      if (!test_tsort(fd)) test_fail;
      if (!test_nanflow(fd, dem)) test_fail;
      if (!test_no_boundary_edges(fd)) test_fail;
      if (!test_outlets(fd)) test_fail;
      if (!test_sum_edge_weights(fd)) test_fail;

      FlowFree(&fd);
    }

    test("Diagonal") {
      ptrdiff_t dims[2] = {3, 3};
      float data[9] = {2, 3, 4, 1, 2, 3, 0, 1, 2};
      Grid dem = GridCreate(GridF32, data, 1.0, dims);
      Flow fd = dinf_flow(dem, 0);

      if (!test_node_once(fd)) test_fail;
      if (!test_tsort_node(fd)) test_fail;
      if (!test_tsort(fd)) test_fail;
      if (!test_nanflow(fd, dem)) test_fail;
      if (!test_no_boundary_edges(fd)) test_fail;
      if (!test_outlets(fd)) test_fail;
      if (!test_sum_edge_weights(fd)) test_fail;

      FlowFree(&fd);
    }

    test("Sink") {
      ptrdiff_t dims[2] = {3, 3};
      float data[9] = {2, 2, 2, 2, 1, 2, 2, 2, 2};
      Grid dem = GridCreate(GridF32, data, 1.0, dims);
      Flow fd = dinf_flow(dem, 0);

      if (!test_node_once(fd)) test_fail;
      if (!test_tsort_node(fd)) test_fail;
      if (!test_tsort(fd)) test_fail;
      if (!test_nanflow(fd, dem)) test_fail;
      if (!test_no_boundary_edges(fd)) test_fail;
      if (!test_outlets(fd)) test_fail;
      if (!test_sum_edge_weights(fd)) test_fail;

      FlowFree(&fd);
    }

    test("Outward sloping cone") {
      ptrdiff_t dims[2] = {16, 16};
      float cellsize = 150.0 / (dims[0] - 1);
      Grid dem = GridCreate(GridF32, NULL, cellsize, dims);
      for (ptrdiff_t j = 0; j < dims[1]; j++) {
        for (ptrdiff_t i = 0; i < dims[0]; i++) {
          float x = cellsize * i - 75.0;
          float y = cellsize * j - 75.0;
          ((float *)dem.data)[j * dims[0] + i] = 200 - sqrtf(x * x + y * y);
        }
      }

      Flow fd = dinf_flow(dem, 0);

      if (!test_node_once(fd)) test_fail;
      if (!test_tsort_node(fd)) test_fail;
      if (!test_tsort(fd)) test_fail;
      if (!test_nanflow(fd, dem)) test_fail;
      if (!test_no_boundary_edges(fd)) test_fail;
      if (!test_outlets(fd)) test_fail;
      if (!test_sum_edge_weights(fd)) test_fail;

      FlowFree(&fd);
      GridFree(&dem);
    }

    test("Outward sloping cone with flats") {
      ptrdiff_t dims[2] = {16, 16};
      float cellsize = 150.0 / (dims[0] - 1);
      Grid dem = GridCreate(GridF32, NULL, cellsize, dims);
      for (ptrdiff_t j = 0; j < dims[1]; j++) {
        for (ptrdiff_t i = 0; i < dims[0]; i++) {
          float x = cellsize * i - 75.0;
          float y = cellsize * j - 75.0;
          ((float *)dem.data)[j * dims[0] + i] =
              fminf(150.0, 200 - sqrtf(x * x + y * y));
        }
      }

      Flow fd = dinf_flow(dem, 0);

      if (!test_node_once(fd)) test_fail;
      if (!test_tsort_node(fd)) test_fail;
      if (!test_tsort(fd)) test_fail;
      if (!test_nanflow(fd, dem)) test_fail;
      if (!test_no_boundary_edges(fd)) test_fail;
      if (!test_outlets(fd)) test_fail;
      if (!test_sum_edge_weights(fd)) test_fail;

      FlowFree(&fd);
      GridFree(&dem);
    }

    if (*snapshot_directory) {
      test("bigtujunga") {
        const char *filepath = "/bigtujunga_30m/dem.tif";
        char *path = calloc(strlen(snapshot_directory) + strlen(filepath) + 1,
                            sizeof *path);
        strcpy(path, snapshot_directory);
        strncat(path, filepath, strlen(filepath) + 1);
        printf("# %s\n", path);

        Grid dem = GridFromFile(path);
        free(path);
        Grid demf = GridFillsinks(dem);
        Flow fd = dinf_flow(demf, 0);

        if (!test_node_once(fd)) test_fail;
        if (!test_tsort_node(fd)) test_fail;
        if (!test_tsort(fd)) test_fail;
        if (!test_nanflow(fd, dem)) test_fail;
        if (!test_no_boundary_edges(fd)) test_fail;
        if (!test_outlets(fd)) test_fail;
        if (!test_sum_edge_weights(fd)) test_fail;

        FlowFree(&fd);
        GridFree(&demf);
        GridFree(&dem);
      }

      test("perfectworld") {
        const char *filepath = "/perfectworld/dem.tif";
        char *path = calloc(strlen(snapshot_directory) + strlen(filepath) + 1,
                            sizeof *path);
        strcpy(path, snapshot_directory);
        strncat(path, filepath, strlen(filepath) + 1);
        printf("# %s\n", path);

        Grid dem = GridFromFile(path);
        free(path);
        Grid demf = GridFillsinks(dem);
        Flow fd = dinf_flow(demf, 0);

        if (!test_node_once(fd)) test_fail;
        if (!test_tsort_node(fd)) test_fail;
        if (!test_tsort(fd)) test_fail;
        if (!test_nanflow(fd, dem)) test_fail;
        if (!test_no_boundary_edges(fd)) test_fail;
        if (!test_outlets(fd)) test_fail;
        if (!test_sum_edge_weights(fd)) test_fail;

        FlowFree(&fd);
        GridFree(&demf);
        GridFree(&dem);
      }

      test("tibet") {
        const char *filepath = "/tibet/dem.tif";
        char *path = calloc(strlen(snapshot_directory) + strlen(filepath) + 1,
                            sizeof *path);
        strcpy(path, snapshot_directory);
        strncat(path, filepath, strlen(filepath) + 1);
        printf("# %s\n", path);

        Grid dem = GridFromFile(path);
        free(path);
        Grid demf = GridFillsinks(dem);
        Flow fd = dinf_flow(demf, 0);

        if (!test_node_once(fd)) test_fail;
        if (!test_tsort_node(fd)) test_fail;
        if (!test_tsort(fd)) test_fail;
        if (!test_nanflow(fd, dem)) test_fail;
        if (!test_no_boundary_edges(fd)) test_fail;
        if (!test_outlets(fd)) test_fail;
        if (!test_sum_edge_weights(fd)) test_fail;

        FlowFree(&fd);
        GridFree(&demf);
        GridFree(&dem);
      }

      test("kedarnath") {
        const char *filepath = "/kedarnath/dem.tif";
        char *path = calloc(strlen(snapshot_directory) + strlen(filepath) + 1,
                            sizeof *path);
        strcpy(path, snapshot_directory);
        strncat(path, filepath, strlen(filepath) + 1);
        printf("# %s\n", path);

        Grid dem = GridFromFile(path);
        free(path);
        Grid demf = GridFillsinks(dem);
        Flow fd = dinf_flow(demf, 0);

        if (!test_node_once(fd)) test_fail;
        if (!test_tsort_node(fd)) test_fail;
        if (!test_tsort(fd)) test_fail;
        if (!test_nanflow(fd, dem)) test_fail;
        if (!test_no_boundary_edges(fd)) test_fail;
        if (!test_outlets(fd)) test_fail;
        if (!test_sum_edge_weights(fd)) test_fail;

        FlowFree(&fd);
        GridFree(&demf);
        GridFree(&dem);
      }

      test("kunashiri") {
        const char *filepath = "/kunashiri/dem.tif";
        char *path = calloc(strlen(snapshot_directory) + strlen(filepath) + 1,
                            sizeof *path);
        strcpy(path, snapshot_directory);
        strncat(path, filepath, strlen(filepath) + 1);
        printf("# %s\n", path);

        Grid dem = GridFromFile(path);
        free(path);
        Grid demf = GridFillsinks(dem);
        Flow fd = dinf_flow(demf, 0);

        if (!test_node_once(fd)) test_fail;
        if (!test_tsort_node(fd)) test_fail;
        if (!test_tsort(fd)) test_fail;
        if (!test_nanflow(fd, dem)) test_fail;
        if (!test_no_boundary_edges(fd)) test_fail;
        if (!test_outlets(fd)) test_fail;
        if (!test_sum_edge_weights(fd)) test_fail;

        FlowFree(&fd);
        GridFree(&demf);
        GridFree(&dem);
      }

      test("taiwan") {
        const char *filepath = "/taiwan/dem.tif";
        char *path = calloc(strlen(snapshot_directory) + strlen(filepath) + 1,
                            sizeof *path);
        strcpy(path, snapshot_directory);
        strncat(path, filepath, strlen(filepath) + 1);
        printf("# %s\n", path);

        Grid dem = GridFromFile(path);
        free(path);
        Grid demf = GridFillsinks(dem);
        Flow fd = dinf_flow(demf, 0);

        if (!test_node_once(fd)) test_fail;
        if (!test_tsort_node(fd)) test_fail;
        if (!test_tsort(fd)) test_fail;
        if (!test_nanflow(fd, dem)) test_fail;
        if (!test_no_boundary_edges(fd)) test_fail;
        if (!test_outlets(fd)) test_fail;
        if (!test_sum_edge_weights(fd)) test_fail;

        FlowFree(&fd);
        GridFree(&demf);
        GridFree(&dem);
      }

      test("taalvolcano") {
        const char *filepath = "/taalvolcano/dem.tif";
        char *path = calloc(strlen(snapshot_directory) + strlen(filepath) + 1,
                            sizeof *path);
        strcpy(path, snapshot_directory);
        strncat(path, filepath, strlen(filepath) + 1);
        printf("# %s\n", path);

        Grid dem = GridFromFile(path);
        free(path);

        Grid demf = GridFillsinks(dem);
        Flow fd = dinf_flow(demf, 0);

        if (!test_node_once(fd)) test_fail;
        if (!test_tsort_node(fd)) test_fail;
        if (!test_tsort(fd)) test_fail;
        if (!test_nanflow(fd, dem)) test_fail;
        if (!test_no_boundary_edges(fd)) test_fail;
        if (!test_outlets(fd)) test_fail;
        if (!test_sum_edge_weights(fd)) test_fail;

        FlowFree(&fd);
        GridFree(&demf);
        GridFree(&dem);
      }
    }

    test("nans") {
      ptrdiff_t dims[2] = {5, 5};
      ptrdiff_t cellsize = 1.0;
      Grid dem = GridCreate(GridF32, NULL, cellsize, dims);

      for (ptrdiff_t j = 0; j < dims[1]; j++) {
        for (ptrdiff_t i = 0; i < dims[0]; i++) {
          float x = cellsize * i;
          float y = cellsize * j;
          ((float *)dem.data)[j * dims[0] + i] = x + y;
        }
      }
      ((float *)dem.data)[2 * dims[0] + 2] = NAN;

      Flow fd = dinf_flow(dem, 0);

      if (!test_node_once(fd)) test_fail;
      if (!test_tsort_node(fd)) test_fail;
      if (!test_tsort(fd)) test_fail;
      if (!test_nanflow(fd, dem)) test_fail;
      if (!test_no_boundary_edges(fd)) test_fail;
      if (!test_outlets(fd)) test_fail;
      if (!test_sum_edge_weights(fd)) test_fail;

      FlowFree(&fd);
      GridFree(&dem);
    }

    test("Memory order slope") {
      ptrdiff_t dims[2] = {5, 7};

      Grid fdem = GridCreate(GridF32, NULL, 1.0, dims);
      for (ptrdiff_t j = 0; j < dims[1]; j++) {
        for (ptrdiff_t i = 0; i < dims[0]; i++) {
          float z = j * dims[0] + i;
          ((float *)fdem.data)[j * dims[0] + i] = z;
        }
      }

      Grid cdem = GridTranspose(fdem);

      Grid fdemf = GridFillsinks(fdem);
      GridFree(&fdem);

      Grid cdemf = GridFillsinks(cdem);
      GridFree(&cdem);

      Flow cfd = dinf_flow(cdemf, 1);
      Flow ffd = dinf_flow(fdemf, 0);

      Grid cacc = FlowAccumulation(cfd);
      Grid facc = FlowAccumulation(ffd);
      Grid caccT = GridTranspose(cacc);

      if (!GridAllClose(caccT, facc, 1e-5)) test_fail;

      GridFree(&caccT);
      GridFree(&cacc);
      GridFree(&facc);
      GridFree(&fdemf);
      GridFree(&cdemf);

      FlowFree(&ffd);
      FlowFree(&cfd);
    }

    test("Memory order saddle (ties)") {
      ptrdiff_t dims[2] = {3, 3};
      ptrdiff_t fdims[2] = {dims[0], dims[1]};

      float fdem_buffer[9] = {2, 2, 0, 2, 1, 2, 0, 2, 2};

      Grid fdem = GridCreate(GridF32, fdem_buffer, 1.0, fdims);
      Grid cdem = GridTranspose(fdem);

      Flow cfd = dinf_flow(cdem, 1);
      Flow ffd = dinf_flow(fdem, 0);

      Grid cacc = FlowAccumulation(cfd);
      Grid facc = FlowAccumulation(ffd);
      Grid caccT = GridTranspose(cacc);

      if (!GridAllClose(caccT, facc, 1e-5)) test_fail;

      GridFree(&caccT);
      GridFree(&cacc);
      GridFree(&facc);
      GridFree(&cdem);

      FlowFree(&ffd);
      FlowFree(&cfd);
    }

    test("Memory order random") {
      ptrdiff_t dims[2] = {19, 17};

      Grid fdem = GridRandom(0x746f706f7f6f6f6c, 1.0, dims);

      Grid cdem = GridTranspose(fdem);

      Grid fdemf = GridFillsinks(fdem);
      GridFree(&fdem);
      Grid cdemf = GridFillsinks(cdem);
      GridFree(&cdem);

      Flow cfd = dinf_flow(cdemf, 1);
      Flow ffd = dinf_flow(fdemf, 0);

      Grid cacc = FlowAccumulation(cfd);
      Grid facc = FlowAccumulation(ffd);

      Grid caccT = GridTranspose(cacc);

      if (!GridAllClose(caccT, facc, 1e-5)) test_fail;

      GridFree(&caccT);
      GridFree(&cacc);
      GridFree(&facc);

      FlowFree(&ffd);
      FlowFree(&cfd);

      GridFree(&fdemf);
      GridFree(&cdemf);
    }
  }
}
