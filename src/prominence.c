#define TOPOTOOLBOX_BUILD

#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>

#include "morphology/reconstruct.h"
#include "topotoolbox.h"

// linked list to handle dynamic amount of peaks
typedef struct List_Node {
  float value;
  ptrdiff_t index;
  struct List_Node *next;
} List_Node;

TOPOTOOLBOX_API
void prominence(float **result_values, ptrdiff_t **result_indices,
                ptrdiff_t *result_size, float *dem, float tolerance,
                ptrdiff_t dims[2]) {
  ptrdiff_t size = dims[0] * dims[1];

  float min_dem_val = INFINITY;
  for (ptrdiff_t i = 0; i < size; i++) {
    if (min_dem_val > dem[i]) min_dem_val = dem[i];
  }

  float *P = malloc(size * sizeof(float));
  if (!P) {
    return;
  }
  for (ptrdiff_t i = 0; i < size; i++) {
    // P = GRIDobj(DEM)+min(DEM);
    // P[i] = dem[i] + min_dem_val;
    P[i] = min_dem_val;
  }

  ptrdiff_t counter = 0;
  List_Node *p_head = NULL;
  List_Node *p_tail = NULL;

  do {
    // [p(counter),ix(counter)] = max(DEM-P);
    float max_val = -INFINITY;
    ptrdiff_t max_index = 0;
    for (ptrdiff_t i = 0; i < size; i++) {
      if (max_val < dem[i] - P[i]) {
        max_val = dem[i] - P[i];
        max_index = i;
      }
    }
    // create new node
    List_Node *new_node = (List_Node *)malloc(sizeof(List_Node));
    new_node->value = max_val;
    new_node->index = max_index;
    new_node->next = NULL;

    // append new node at end of list
    if (p_head == NULL) {
      p_head = new_node;
      p_tail = new_node;
    } else {
      p_tail->next = new_node;
      p_tail = new_node;
    }

    // P.Z(ix) = DEM.Z(ix);
    // replace new found max value with original dem value
    P[p_tail->index] = dem[p_tail->index];

    // P.Z = imreconstruct(P.Z,DEM.Z);
    reconstruct(P, dem, dims);
    counter++;

  } while (p_tail->value > tolerance);
  free(P);

  // turn linked list into arrays containing results and free the linked list
  *result_values = malloc(counter * sizeof(float));
  *result_indices = malloc(counter * sizeof(ptrdiff_t));
  if (!result_values || !result_indices) {
    return;
  }
  List_Node *current = p_head;
  List_Node *next;

  for (ptrdiff_t i = 0; i < counter; i++) {
    // save values and indices in return arrays
    (*result_values)[i] = current->value;
    (*result_indices)[i] = current->index;

    // free the current node
    next = current->next;
    free(current);
    current = next;
  }

  *result_size = counter;
}