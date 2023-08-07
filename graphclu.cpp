/*******************************************************************************
 *
 * This file is part of TODO software.
 * Copyright (C) 2015-2018 Sami Sieranoja
 * <samisi@uef.fi>, <sami.sieranoja@gmail.com>
 *
 * TODO is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version. You should have received a copy
 * of the GNU Lesser General Public License along with TODO.
 * If not, see <http://www.gnu.org/licenses/lgpl.html>.
 *******************************************************************************/

#include <stdio.h>
#include <limits.h>
#include <cstring>
#include <pthread.h>
#include "contrib/argtable3.h"

#include <execinfo.h>
#include <signal.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdio.h>
#include <iostream>
#include <iostream>
#include <cfloat>
#include <cmath>
#include <bits/stdc++.h>
#include <sys/stat.h>
#include <errno.h>

using namespace std;

#define SIMILARITY 0
#define DISTANCE 1

// #define USE_IGRAPH 1

#ifdef USE_IGRAPH
extern "C" {
#include <igraph.h>
#include <unistd.h>
#include <libgen.h>
}
#endif

#include "heap.cpp"
#include "linked_list.hpp"
#include "nngraph.hpp"

#include "gclu_options.h"

#include "timer.h"
#include "util.h"

// #include "linked_list.h"
// #include "nngraph.h"


#include "linked_list.cpp"
#include "nngraph.cpp"

#include "graphclu.h"

#ifdef USE_IGRAPH
#include "igraph_algos.h"
#endif

// double costf_w(double intsum, int size, nnGraph *graph, Clustering *clu);

int g_iter = 0;

void handler(int sig) {
  void *array[10];
  size_t size;

  // get void*'s for all entries on the stack
  size = backtrace(array, 10);

  // print out all the frames to stderr
  fprintf(stderr, "Error: signal %d:\n", sig);
  backtrace_symbols_fd(array, size, STDERR_FILENO);
  exit(1);
}

void write_ints_to_file(const char *fn, int *data, int N) {
  int i;

  FILE *fp;
  fp = fopen(fn, "w");
  for (int i = 0; i < N; i++) {
    fprintf(fp, "%d", data[i]);
    fprintf(fp, "\n");
  }
  fclose(fp);
}

void write_ints_to_fp(FILE *fp, int *data, int N) {
  int i;
  for (int i = 0; i < N; i++) {
    fprintf(fp, "%d", data[i]);
    fprintf(fp, "\n");
  }
}

// Return integers from 0 to rangeEnd in random order.
int *getRandomOrderInts(int rangeEnd) {
  int *randSample = (int *)calloc((rangeEnd + 1), sizeof(int)); // Init with zeros
  int i, tmp;
  int swapwith;

  for (i = 0; i <= rangeEnd; i++) {
    randSample[i] = i;
  }

  for (i = 0; i <= rangeEnd; i++) {
    swapwith = rand() % (rangeEnd + 1);
    tmp = randSample[swapwith];
    randSample[swapwith] = randSample[i];
    randSample[i] = tmp;
  }
  return randSample;
}

void init_Clustering(Clustering **_clu, int N, int K) {
  *_clu = (Clustering *)malloc(sizeof(Clustering));
  Clustering *clu = *_clu;
  clu->N = N;
  clu->K = K;
  clu->part = (int *)malloc(sizeof(int) * N);
  clu->clusize = (int *)malloc(sizeof(int) * K);
  clu->ntr_sums = (double *)malloc(sizeof(double) * K);
  clu->ext_sums = (double *)malloc(sizeof(double) * K);
  clu->total_sums = (double *)malloc(sizeof(double) * K);
  clu->costs = (double *)malloc(sizeof(double) * K);

  clu->node_to_clu_w = (double **)malloc(sizeof(double *) * K);
  for (int i = 0; i < K; i++) {
    clu->node_to_clu_w[i] = (double *)calloc(N, sizeof(double));
  }

  // Initialize to random partition
  for (int i = 0; i < N; i++) {
    clu->part[i] = rand() % K;
  }
}

void find_max_val(double *data, int N, int *ret_ind, double *ret_val) {
  double maxval = -DBL_MAX;
  int maxind = 0;
  for (int i = 0; i < N; i++) {
    if (maxval < data[i]) {
      maxval = data[i];
      maxind = i;
    }
  }
  *ret_ind = maxind;
  *ret_val = maxval;
}

int select_random_node_from_clu(Clustering *clu, int clu_id) {
  int *rand_order;
  int rand_node = 0;
  // TODO: this is not a bottleneck, but could be better optimized
  rand_order = getRandomOrderInts(clu->N - 1);
  for (int i = 0; i < clu->N; i++) {
    int id = rand_order[i];
    if (clu->part[id] == clu_id) {
      rand_node = id;
      break;
    }
  }
  return rand_node;
}

// Merge two clusters, then split one cluster
void merge_and_split(Clustering *clu, nnGraph *graph) {

  // Choose two clusters to merge and one cluster to split

  int *rand_order;
  int clu_to_split;
  int to_mergeA = rand() % clu->K;
  int to_mergeB;
  int seed_id, steps;
  gNode *node;

  char partfn[1000];
  if (g_opt.debug_save_intermediate_part >= 2) {
    sprintf(partfn, "part/part-%04d-%03d-ms0.txt", g_iter, 0);
    write_ints_to_file(partfn, clu->part, clu->N);
  }

  // Select another cluster to merge with to_mergeA
  // This is failsafe strategy in case the probabilistic selection doesn't work
  while (1) {
    to_mergeB = rand() % clu->K;

    if (to_mergeB != to_mergeA) {
      break;
    }
  }

  // Select two clusters to merge, probabilistically with pdf formed by external weights
  // In the case of all clusters being completely separate (total_ext_sum=0), then select
  // two random clusters (above)
  // TODO: this could be faster

  float total_ext_sum = 0.0;
  for (int i = 0; i < clu->K; i++) {
    total_ext_sum += clu->ext_sums[i];
  }
  float rext = RAND_IN_RANGE(0.00, total_ext_sum);
  float wextsum = 0.0;

  for (int i = 0; i < clu->N && wextsum < rext && total_ext_sum > 0.0; i++) {
    node = &graph->nodes[i];
    for (int j = 0; j < node->neighbors->size; j++) {
      gItem *gi = (gItem *)ll_get_item(node->neighbors, j);
      // Not in same cluster
      if (clu->part[gi->id] != clu->part[i]) {
        wextsum += gi->dist;
        if (wextsum > rext) {
          to_mergeA = clu->part[gi->id];
          to_mergeB = clu->part[i];
          break;
        }
      }
    }
  }

  // printf("Merge %d with %d\n", to_mergeA, to_mergeB);
  // Perform merge
  for (int i = 0; i < clu->N; i++) {
    if (clu->part[i] == to_mergeB) {
      clu->part[i] = to_mergeA;
    }
  }
  clu->clusize[to_mergeA] += clu->clusize[to_mergeB];
  clu->clusize[to_mergeB] = 0;
  // END of merge

  if (g_opt.debug_save_intermediate_part >= 2) {
    sprintf(partfn, "part/part-%04d-%03d-ms1.txt", g_iter, 0);
    write_ints_to_file(partfn, clu->part, clu->N);
  }

  // Select cluster to split.
  // After merging, to_mergeB will be empty, so we can't split it
  while (1) {
    clu_to_split = rand() % clu->K;
    if (clu_to_split != to_mergeB) {
      break;
    }
  }

  // Find random node from cluster clu_to_split
  seed_id = select_random_node_from_clu(clu, clu_to_split);

  // Choose new cluster size between 5% and 95% of size of cluster_to_split

  float randf = RAND_IN_RANGE(0.05, 0.95);
  steps = (int)(randf * (clu->clusize[clu_to_split]));
  // steps = 100;

  // Since all to_mergeB have been changed to to_mergeA, the label to_mergeB is free to use for new
  // cluster
  grow_Cluster(clu, graph, seed_id, to_mergeB /* newpart*/, steps);

  if (g_opt.debug_save_intermediate_part >= 2) {
    sprintf(partfn, "part/part-%04d-%03d-ms2.txt", g_iter, 0);
    write_ints_to_file(partfn, clu->part, clu->N);
  }
}

// Starting from seed node, expand cluster to other nodes, in each step maximizing total internal
// weight of cluster
void grow_Cluster(Clustering *clu, nnGraph *graph, int seed_id, int newpart, int steps) {
  double maxval = 0.0;
  int maxnode = 0;
  gNode *seed;
  gNode *node;
  seed = &graph->nodes[seed_id];
  // printf("Grow cluster on: %d\n",seed_id);

  for (int i = 0; i < clu->N; i++) {
    clu->node_to_clu_w[newpart][i] = 0.0;

    // Dissolve cluster
    if (clu->part[i] == newpart) {
      clu->part[i] = rand() % clu->K;
    }
  }

  clu->part[seed_id] = newpart;
  clu->node_to_clu_w[newpart][seed_id] = -DBL_MAX;

  for (int j = 0; j < seed->neighbors->size; j++) {
    gItem *gi = (gItem *)ll_get_item(seed->neighbors, j);
    clu->node_to_clu_w[newpart][gi->id] += gi->dist;
  }

  for (int i_step = 0; i_step < steps; i_step++) {
    find_max_val(clu->node_to_clu_w[newpart], clu->N, &maxnode, &maxval);
    if (maxval <= 0.0) {
      break;
      // printf("i_tesp=%d maxval==0.0\n",i_step);
    }
    // if(maxval < 0.0) {} //TODO
    node = &graph->nodes[maxnode];

    for (int j = 0; j < node->neighbors->size; j++) {
      gItem *gi = (gItem *)ll_get_item(node->neighbors, j);
      if (clu->part[gi->id] != newpart) {
        clu->node_to_clu_w[newpart][gi->id] += gi->dist;
      }
    }
    clu->part[maxnode] = newpart;
    clu->node_to_clu_w[newpart][maxnode] = -DBL_MAX;
  }
}

void free_Clustering(Clustering *clu) {
  free(clu->part);
  free(clu->clusize);
  free(clu->ntr_sums);
  free(clu->ext_sums);
  free(clu->total_sums);
  free(clu->costs);

  for (int i = 0; i < clu->K; i++) {
    free(clu->node_to_clu_w[i]);
  }
  free(clu->node_to_clu_w);

  free(clu);
}

void copy_Clustering(Clustering *clu, Clustering *newclu) {
  newclu->N = clu->N;
  newclu->K = clu->K;
  newclu->cost = clu->cost;
  newclu->conductance = clu->conductance;
  newclu->balance_factor = clu->balance_factor;
  newclu->min_part_size = clu->min_part_size;
  newclu->max_part_size = clu->max_part_size;
  newclu->min_max_part_ratio = clu->min_max_part_ratio;

  memcpy(newclu->part, clu->part, sizeof(int) * clu->N);
  memcpy((newclu->clusize), (clu->clusize), sizeof(int) * clu->K);
  memcpy((newclu->ntr_sums), (clu->ntr_sums), sizeof(double) * clu->K);
  memcpy((newclu->ext_sums), (clu->ext_sums), sizeof(double) * clu->K);
  memcpy((newclu->total_sums), (clu->total_sums), sizeof(double) * clu->K);
  memcpy((newclu->costs), (clu->costs), sizeof(double) * clu->K);
}

void clone_Clustering(Clustering *clu, Clustering **newclu_) {
  init_Clustering(newclu_, clu->N, clu->K);
  Clustering *newclu = *newclu_;

  copy_Clustering(clu, newclu);
}

void track_min_int(int *setv, int val) {
  if (*setv > val) {
    *setv = val;
  }
}

void track_max_int(int *setv, int val) {
  if (*setv < val) {
    *setv = val;
  }
}

// Return by modifying *clu
double costf(nnGraph *graph, Clustering *clu, double *r_conductance) {
  gNode *nodeA;
  gNode *nodeB;
  int partA, partB;
  double ntr_sum = 0.0;
  double ntr_sum2 = 0.0;
  double ext_sum = 0.0;
  double conductance = 0.0;
  double balance_factor = 0.0;

  clu->min_part_size = INT_MAX;
  clu->max_part_size = 0;

  for (int i = 0; i < clu->K; i++) {
    clu->ntr_sums[i] = 0.0;
    clu->ext_sums[i] = 0.0;
    clu->clusize[i] = 0;
    clu->total_sums[i] = 0.0;
  }

  for (int i = 0; i < graph->size; i++) {
    nodeA = &graph->nodes[i];
    partA = clu->part[i];
    clu->clusize[partA]++;
    clu->total_sums[partA] += nodeA->weight_sum;
    // printf("nodeA->weight_sum=%f\n", nodeA->weight_sum);
    for (int j = 0; j < nodeA->neighbors->size; j++) {
      gItem *gi = (gItem *)ll_get_item(nodeA->neighbors, j);
      partB = clu->part[gi->id];
      if (partA == partB) {
        clu->ntr_sums[partA] += gi->dist;
      } else {
        clu->ext_sums[partA] += gi->dist;
        clu->ext_sums[partB] += gi->dist;
      }
    }
  }

  for (int i = 0; i < clu->K; i++) {
    double clu_ext_sum = clu->total_sums[i] - clu->ntr_sums[i];
    conductance += clu_ext_sum / clu->total_sums[i];
    ntr_sum += clu->ntr_sums[i] / clu->clusize[i];
    double clu_cost = costf_w(clu->ntr_sums[i], clu->clusize[i], graph, clu);
    clu->costs[i] = clu_cost;
    ntr_sum2 += clu_cost;
    double bf_part = (clu->N / ((double)clu->K)) / clu->clusize[i];
    if (bf_part < 1.0) {
      bf_part = 1.0 / bf_part;
    }
    balance_factor += bf_part / clu->K;
    track_min_int(&(clu->min_part_size), clu->clusize[i]);
    track_max_int(&(clu->max_part_size), clu->clusize[i]);
  }
  conductance = -(conductance / clu->K);
  ntr_sum = ntr_sum / clu->K;
  ntr_sum2 = ntr_sum2 / clu->K;
  clu->min_max_part_ratio = clu->min_part_size / ((double)clu->max_part_size);

  *r_conductance = conductance;
  clu->balance_factor = balance_factor;
  clu->conductance = conductance;
  if (g_opt.costf == 0) {
    return ntr_sum;
  } else if (g_opt.costf == 1) {
    return conductance;
  } else if (g_opt.costf >= 2) {
    return ntr_sum2;
  }
}

// Scale weights
// if invert == 1, convert from distances to similarities
// Also calculate node sum of weights and graph total weight
void scale_weights(nnGraph *graph, int scale_type) {
  gNode *nodeA;
  gItem *gi;
  double maxdist = get_max_weight(graph);
  graph->total_weight = 0.0;
  for (int i = 0; i < graph->size; i++) {
    nodeA = &graph->nodes[i];
    gi = (gItem *)ll_get_item(nodeA->neighbors, nodeA->neighbors->size - 1);
    nodeA->weight_sum = 0.0;
    for (int j = 0; j < nodeA->neighbors->size; j++) {
      gi = (gItem *)ll_get_item(nodeA->neighbors, j);
      // Scale distances to similarities
      if (scale_type == 1) {
        gi->dist = (maxdist - gi->dist) / maxdist;
      }
      // Scale similarities to [0,1] range
      else if (scale_type == 2) {
        gi->dist = 1.0 - 1.0 / (1.0 + gi->dist);
      } else {
        // No scale, just calculate weight sums
      }
      nodeA->weight_sum += gi->dist;
    }
    graph->total_weight += nodeA->weight_sum;
  }
}

double costf_w(double intsum, int size, nnGraph *graph, Clustering *clu) {
  double r;

  double p = 3.0;
  if (g_opt.costf == 2) {

    // "Expected" weight divided by observed
    if (size == 0 || intsum == 0.0) {
      return -20e1;
    }
    r = -(graph->total_weight / ((double)(clu->K))) / (pow(1.0e-10 + intsum, 1.0));

  } else if (g_opt.costf == 3) {
    if (size == 0 || intsum == 0.0) {
      return -20e1;
    }
    r = -(size / (pow(1.0e-10 + intsum, 1.0)));
  } else if (g_opt.costf == 4) {
    if (size == 0) {
      return 0.0;
    }
    r = intsum / size;
  } else if (g_opt.costf == 5) {
    r = intsum / size;
  } else if (g_opt.costf == 10) {
    r = -(pow(size, 0.5) / (pow(1.0e-10 + intsum, 1.0)));
  }
  return r;
}

int choose_best_by_delta(
    // INPUT:
    int nid, nnGraph *graph, Clustering *clu,
    // OUTPUT:
    double *d_cost, double *d_sum_oldpart, double *d_sum_newpart, double *d_total_sum) {
  gNode *nodeA;
  int partA, partB;
  double *d_ntr_sums;
  d_ntr_sums = (double *)calloc(clu->K, sizeof(double));
  double costdiff, bestcost;
  int bestpart;

  nodeA = &(graph->nodes[nid]);
  partA = clu->part[nid];

  for (int j = 0; j < nodeA->neighbors->size; j++) {
    gItem *gi = (gItem *)ll_get_item(nodeA->neighbors, j);
    partB = clu->part[gi->id];
    d_ntr_sums[partB] += gi->dist * 2;
  }

  double oldcost, newcost, d_removal, d_add;

  // Calculate removal cost from current cluster
  if (g_opt.costf == 1) {
    oldcost = -(clu->total_sums[partA] - clu->ntr_sums[partA]) / clu->total_sums[partA];
    newcost = -((clu->total_sums[partA] - nodeA->weight_sum) -
                (clu->ntr_sums[partA] - d_ntr_sums[partA])) /
              (clu->total_sums[partA] - nodeA->weight_sum);

  } else if (g_opt.costf >= 2) {
    oldcost = costf_w(clu->ntr_sums[partA], clu->clusize[partA], graph, clu);
    newcost =
        costf_w(clu->ntr_sums[partA] - d_ntr_sums[partA], clu->clusize[partA] - 1, graph, clu);
  }
  d_removal = newcost - oldcost;

  // Calculate addition cost to new cluster and total delta
  bestpart = partA;
  bestcost = -9999999999.0;
  for (int i_clu = 0; i_clu < clu->K; i_clu++) {
    if (i_clu == partA) {
      continue;
    }

    if (g_opt.costf == 1) {
      // Conductance
      oldcost = -((clu->total_sums[i_clu] - clu->ntr_sums[i_clu]) / clu->total_sums[i_clu]);
      newcost = -((clu->total_sums[i_clu] + nodeA->weight_sum) -
                  (clu->ntr_sums[i_clu] + d_ntr_sums[i_clu])) /
                (clu->total_sums[i_clu] + nodeA->weight_sum);
    } else if (g_opt.costf >= 2) {
      oldcost = costf_w(clu->ntr_sums[i_clu], clu->clusize[i_clu], graph, clu);
      newcost =
          costf_w(clu->ntr_sums[i_clu] + d_ntr_sums[i_clu], clu->clusize[i_clu] + 1, graph, clu);
    }

    d_add = newcost - oldcost;
    costdiff = d_add + d_removal;

    if (costdiff > bestcost) {
      bestcost = costdiff;
      bestpart = i_clu;
    }
    // TODO: prefer current part in tie cases
  }

  bestcost = bestcost / clu->K;

  if (bestcost <= 0.0) {
    bestpart = partA;
    (*d_cost) = 0.0;
    *d_sum_oldpart = 0.0;
    *d_sum_newpart = 0.0;
    *d_total_sum = 0.0;
  } else {
    // Change in internal sums of old and new clusters
    *d_sum_oldpart = -d_ntr_sums[partA];
    *d_sum_newpart = d_ntr_sums[bestpart];
    // Change in total sums of old and new clusters (the same)
    *d_total_sum = nodeA->weight_sum;
    // Cost function delta for the best found cluster
    (*d_cost) = bestcost;
  }

  free(d_ntr_sums);
  return bestpart;
}

bool compare_nodeItem(nodeItem i1, nodeItem i2) { return (i1.density > i2.density); }

// Create an initial clustering by growing clusters in dense parts of the graph
void density_init_partition(nnGraph *graph, Clustering *clu) {
  int grow_size;
  gNode *nodeA, *nodeB;
  gItem *gi;
  nodeItem *items = (nodeItem *)malloc(sizeof(nodeItem) * graph->size);

  // Calculate densities for each node
  for (int i = 0; i < graph->size; i++) {
    nodeA = &graph->nodes[i];
    // gi = (gItem *)ll_get_item(nodeA->neighbors, nodeA->neighbors->size - 1);
    // nodeA->weight_sum = 0.0;

    items[i].id = i;
    items[i].density = 0;
    if (g_opt.density_method == 1) {
      for (int j = 0; j < nodeA->neighbors->size; j++) {
        gItem *gi = (gItem *)ll_get_item(nodeA->neighbors, j);
        nodeB = &graph->nodes[gi->id];
        // Weight to neighbor multiplied by neighbors total weight
        items[i].density += (gi->dist) * (nodeB->weight_sum);
      }
    } else if (g_opt.density_method == 2) {
      // Total weight of node
      items[i].density = nodeA->weight_sum;

    } else {
      terminal_error("No such density method");
    }

    if (g_opt.verbose >= 2.0) {
      if (i < 20) {
        printf("nid=%d dens=%f\n", nodeA->id, items[nodeA->id].density);
      }
    }
  }
  // Sort nodes by density
  sort(items, items + graph->size, compare_nodeItem);

  if (g_opt.verbose >= 2.0) {
    printf("Sorted by density:\n");
  }
  for (int i = 0; i < graph->size; i++) {
    nodeA = &graph->nodes[items[i].id];
    // gi = (gItem *)ll_get_item(nodeA->neighbors, nodeA->neighbors->size - 1);
    // nodeA->weight_sum = 0.0;

    if (g_opt.verbose >= 2.0) {
      if (i < 20) {
        printf("nid=%d dens=%f\n", nodeA->id, items[i].density);
      }
    }
  }

  // Set all partition ids to -1
  for (int i_node = 0; i_node < graph->size; i_node++) {
    clu->part[i_node] = -1;
  }

  // Grow K clusters, starting each from the node of highest density
  // that is not yet joined to any partition.
  grow_size = (int)((clu->N * g_opt.grow_factor) / ((float)clu->K));
  for (int i_clu = 0; i_clu < clu->K; i_clu++) {
    for (int i_node = 0; i_node < graph->size; i_node++) {
      int nid = items[i_node].id; // items sorted by density
      if (clu->part[nid] == -1) {
        if (g_opt.verbose >= 2.0) {
          printf("Grow clu=%d from=%d (%d'th) dens=%f\n", i_clu, items[i_node].id, i_node,
                 items[i_node].density);
        }
        grow_Cluster(clu, graph, nid, i_clu, grow_size);
        break;
      }
    }
  }

  int num_unhandled = 0;
  // Set unhandled nodes to any random partition
  for (int i_node = 0; i_node < graph->size; i_node++) {
    if (clu->part[i_node] == -1) {
      clu->part[i_node] = rand() % clu->K;
      num_unhandled++;
    }
  }

  if (g_opt.verbose >= 3.0) {
    printf("num_unhandled=%d\n", num_unhandled);
  }
}

int k_algo(nnGraph *graph, Clustering *clu) {
  Clustering *newclu;
  double cost, bestcost, costprev, curcost, curcost2, bestcost2, costerr;
  double d_cost, d_sum_oldpart, d_sum_newpart, d_total_sum, conductance;
  int newpart, oldpart, bestpart, bestpart2;
  init_Clustering(&newclu, graph->size, clu->K);
  copy_Clustering(clu, newclu);
  costprev = costf(graph, clu, &conductance);
  curcost = costprev;
  int brutef_cost = 0;
  int verbose = 2.0;
  int *rand_order;
  int i;
  int iter;

  char partfn[1000];
  d_total_sum = 0.0;

  if (g_opt.debug_save_intermediate_part >= 2) {
    sprintf(partfn, "part/part-%04d-%03d.txt", g_iter, 0);
    write_ints_to_file(partfn, clu->part, clu->N);
  }

  for (iter = 0; iter < g_opt.max_iter; iter++) {
    rand_order = getRandomOrderInts(clu->N - 1);
    for (int i_node = 0; i_node < clu->N; i_node++) {
      // i = rand_order[i_node];
      i = i_node;
      bestpart = 0;
      oldpart = clu->part[i];
      bestcost = -1.0;

      if (brutef_cost) {
        curcost = costf(graph, clu, &conductance);
        curcost2 = costf(graph, newclu, &conductance);

        for (int j = 0; j < clu->K; j++) {
          newclu->part[i] = j;
          cost = costf(graph, newclu, &conductance);
          if (cost > bestcost) {
            bestpart = j;
            bestcost = cost;
          }
        }
      }

      bestpart2 = choose_best_by_delta(i, graph, clu, &d_cost, &d_sum_oldpart, &d_sum_newpart,
                                       &d_total_sum);

      newpart = bestpart2;

      if (brutef_cost) {
        bestcost2 = curcost + d_cost;
        costerr = bestcost2 - bestcost;
        if (g_opt.verbose >= 2.0) {
          printf("i=%d bestcost=%f bestcost2=%f bestpart=%d bestpart2=%d costerr=%f d_cost=%f", i,
                 bestcost, bestcost2, bestpart, bestpart2, costerr, d_cost);
        }
      } else {
        if (newpart != oldpart) {
          curcost = curcost + d_cost;
          clu->ntr_sums[oldpart] += d_sum_oldpart;
          clu->ntr_sums[newpart] += d_sum_newpart;
          clu->clusize[oldpart]--;
          clu->clusize[newpart]++;
          clu->part[i] = newpart;
          newclu->part[i] = newpart;

          clu->total_sums[oldpart] -= d_total_sum;
          clu->total_sums[newpart] += d_total_sum;

          // cost = costf(graph, clu, &conductance);

          // printf("i=%d d_sum_oldpart=%f d_sum_newpart=%f \n", i, d_sum_oldpart, d_sum_newpart);
        } else {
          // printf("nochange i=%d d_sum_oldpart=%f d_sum_newpart=%f \n", i, d_sum_oldpart,
          // d_sum_newpart);
        }
      }

      if (brutef_cost) {
        if (bestpart == oldpart) {
          printf(" NOCHANGE2 ");
        } else {
          printf(" CHANGE2 ");
        }

        if (bestpart == bestpart2) {
          printf(" SAME\n");
        } else {
          printf(" DIFFERENT\n");
        }
        clu->part[i] = bestpart;
        newclu->part[i] = bestpart;
      }

    } // END Loop all nodes
    free(rand_order);

    cost = costf(graph, clu, &conductance);
    clu->cost = cost;

    if (g_opt.debug_save_intermediate_part >= 2) {
      sprintf(partfn, "part/part-%04d-%03d.txt", g_iter, iter + 1);
      write_ints_to_file(partfn, clu->part, clu->N);
    }

    if (g_opt.verbose >= 2.0) {
      printf("I=%d cost=%f curcost=%f conduct=%f\n", iter, cost * (g_opt.costmultip),
             curcost * (g_opt.costmultip), conductance * -1.0);
    }
    if (cost <= costprev) {
      if (g_opt.verbose >= 2.0) {
        // printf("Converged.\n");
      }
      break;
    }
    // curcost=cost;
    costprev = cost;
  }

  free_Clustering(newclu);
  return iter; // number of iterations
}

void m_algo(nnGraph *graph, Clustering *clu, int n_repeats, int n_clusters) {

  printf("Start K/swap algo\n");
  Clustering *newclu;
  clu->cost = -999999.0;
  char partfn[1000];
  int grow_size = 30;
  float mean_iterations = 0.0;
  int num_iter;

  if (g_opt.debug_save_intermediate_part) {
    int status;
    status = mkdir("part", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    if (status == -1) {
      if (errno != EEXIST) {
        terminal_error("Could not create directory 'part'");
      }
    }
  }

  // printf("clusters=%d\n", n_clusters);
  init_Clustering(&newclu, graph->size, n_clusters);
  if (g_opt.density_method >= 1) {
    density_init_partition(graph, newclu);
  }
  num_iter = k_algo(graph, newclu);
  mean_iterations += ((float)num_iter) / (n_repeats + 1);
  printf("K-algo cost=%f\n", 0, newclu->cost * (g_opt.costmultip));
  copy_Clustering(newclu, clu);

  if (g_opt.debug_save_intermediate_part >= 1) {
    sprintf(partfn, "part/part-%04d.txt", g_iter);
    write_ints_to_file(partfn, clu->part, clu->N);
  }

  for (int i_rep = 1; i_rep <= n_repeats; i_rep++) {
    g_iter = i_rep;
    copy_Clustering(clu, newclu);
    float gf = g_opt.grow_factor;
    if (g_opt.grow_factor_start != g_opt.grow_factor_end) {
      // experimental
      gf = RAND_IN_RANGE(g_opt.grow_factor_start, g_opt.grow_factor_end);
      // printf("gf=%f\n",gf);
    }
    grow_size = (int)((clu->N * gf) / ((float)clu->K));
    // grow_Cluster(newclu, graph, rand() % newclu->N /*seed_id*/, rand() % newclu->K /*newpart*/,
    // grow_size);
    merge_and_split(newclu, graph);
    num_iter = k_algo(graph, newclu);
    mean_iterations += ((float)num_iter) / (n_repeats + 1);
    printf(
        "M-algo REP=%d, time[swap%d]=%f cost=%f best=%f balance_factor=%f conductance=%f niter=%d",
        i_rep, i_rep, g_timer.get_time(), newclu->cost * (g_opt.costmultip),
        clu->cost * (g_opt.costmultip), newclu->balance_factor, newclu->conductance * -1.0,
        num_iter);

    if (clu->cost < newclu->cost) {
      copy_Clustering(newclu, clu);
      // sprintf(partfn, "part/part-%04d-999.txt", i_rep);
      // write_ints_to_file(partfn, clu->part, clu->N);
      // printf(" %s ",partfn);
      printf(" improved=1");
      if (g_opt.debug_save_intermediate_part >= 1) {
        sprintf(partfn, "part/part-%04d.txt", g_iter);
        write_ints_to_file(partfn, clu->part, clu->N);
      }
    } else {
      printf(" improved=0");
    }
    printf("\n");

    // if (i_rep % 7 == 0) {
    copy_Clustering(clu, newclu);
    // }

    // free_Clustering(newclu);
  }
  printf("FINAL=1 cost=%f balance_factor=%f conductance=%f TIME=%f costf=%d seed=%d "
         "min_part_size=%d max_part_size=%d min_max_part_ratio=%f, niter=%f\n",
         g_opt.costmultip * clu->cost, clu->balance_factor, clu->conductance * -1.0,
         g_timer.get_time(), g_opt.costf, g_opt.seed, clu->min_part_size, clu->max_part_size,
         clu->min_max_part_ratio, mean_iterations);
}

void repeated_k_algo(nnGraph *graph, Clustering *clu, int n_repeats, int n_clusters) {

  Clustering *newclu;
  // clone_Clustering(clu, &newclu);
  clu->cost = -999999.0;
  for (int i_rep = 0; i_rep < n_repeats; i_rep++) {
    init_Clustering(&newclu, graph->size, n_clusters);
    k_algo(graph, newclu);
    printf("REP = %d, cost=%f best=%f\n", i_rep, newclu->cost, clu->cost);

    // grow_Cluster(newclu, graph, rand() % newclu->N /*seed_id*/,  rand() % newclu->K
    // /*newpart*/, 30); k_algo(graph, newclu); printf("REP_S = %d, cost=%f best=%f\n", i_rep,
    // newclu->cost, clu->cost);

    if (clu->cost < newclu->cost) {
      copy_Clustering(newclu, clu);
    }

    free_Clustering(newclu);
  }
  printf("FINAL cost=%f FINAL\n", clu->cost);

  grow_Cluster(clu, graph, rand() % clu->N /*seed_id*/, rand() % clu->K /*newpart*/, 30);
  k_algo(graph, clu);

  printf("FINAL2 cost=%f FINAL\n", clu->cost);
}

Clustering *read_partition(const char *FileName, nnGraph *graph) // without using TS
{
  int i;
  FILE *f;
  int index;
  int c;
  int N;
  int k;
  int *part;
  Clustering *clu = (Clustering *)malloc(sizeof(Clustering));

  f = fopen(FileName, "r");
  do
    c = getc(f);
  while (c != 10);

  // fscanf(f, "%i\n", &(clu->K));
  // fscanf(f, "%i\n", &(clu->N));

  fscanf(f, "%i\n", &(k));
  fscanf(f, "%i\n", &(N));
  printf("k=%d N=%d\n", k, N);

  if (N < 2 || N > 1e7) {
    return NULL;
  }
  // part=(int*) malloc(sizeof(int)*N);
  init_Clustering(&clu, graph->size, k);
  part = clu->part;

  while (getc(f) == '-') {
  }

  for (i = 0; i < N; i++) {
    c = fscanf(f, "%i", &index);
    if (c == 0) {
      // ErrorMessage("ERROR reading partitioning.\n");
      // ExitProcessing(-1);
    }

    do
      c = getc(f);
    while (c != 10 && c != EOF);

    // if (index > P->PartitionCount || index < 1) {
    // ErrorMessage("ERROR: Invalid partition index (%d) found.\n", index);
    // ExitProcessing(-1);
    // }
    part[i] = index - 1;
    // printf("%d: %d\n", i, index);
  }
  fclose(f);
  return clu;
}

int main(int argc, char *argv[]) {

  signal(SIGSEGV, handler);
  int repeats = 0;
  int clusters = 50;
  int graph_type = SIMILARITY;
  int output_write_header = 1;
  int *part;
  Clustering *clu;

  g_opt.repeats = 0;
  g_opt.dissolve = 1;
  g_opt.clusters = 50;
  g_opt.graph_type = SIMILARITY;
  g_opt.debug = 0;
  g_opt.grow_factor = 0.8;
  g_opt.grow_factor_start = 1.0;
  g_opt.grow_factor_end = 1.0;
  g_opt.verbose = 1;
  g_opt.costf = 2;
  g_opt.minimize = 0; // maximize by default
  g_opt.debug_save_intermediate_part = 0;
  g_opt.density_method = 1;
  g_opt.max_iter = 200;

  struct arg_file *infn;
  struct arg_file *outfn;
  struct arg_file *partfn;
  // struct arg_str *outf;
  struct arg_str *informat;
  struct arg_str *algo;
  struct arg_int *rngSeed;
  struct arg_int *nthreads;
  struct arg_int *a_repeats;
  struct arg_int *a_maxiter;

  struct arg_int *a_clusters;
  struct arg_int *a_costf;
  struct arg_int *a_saveiterparts;
  struct arg_str *a_graphtype;
  struct arg_int *a_density;
  struct arg_str *a_algo;
  struct arg_dbl *a_grow_factor;

  struct arg_dbl *a_grow_factor_start; // Start and end of ranges
  struct arg_dbl *a_grow_factor_end;

  struct arg_int *a_write_header;

  struct arg_dbl *a_verbose;
  struct arg_str *a_scale;
  struct arg_lit *help;
  struct arg_end *end;

  void *argtable[] = {
      help = arg_litn(NULL, "help", 0, 1, "display this help and exit"),
      outfn = arg_filen("o", "out", "<file>", 0, 1, "output file"),
      infn = arg_filen(NULL, NULL, "<file>", 1, 1, "input graph file"),

      rngSeed = arg_intn(NULL, "seed", "<n>", 0, 1, "random number seed"),
      a_repeats = arg_intn("R", "repeats", "<n>", 0, 1, "Number of repeats"),
      a_clusters = arg_intn("K", "clusters", "<n>", 0, 1, "Number of clusters"),
      a_grow_factor = arg_dbln("g", "growf", "<num>", 0, 1, "grow factor [0,1]"),

      a_verbose = arg_dbln("V", "verbose", "<num>", 0, 1, "Verbose level (default = 1)"),

      a_write_header = arg_intn("H", "write-header", "<num>", 0, 1,
                                "Write header for result partition (0 =no, 1=yes (default))"),

      // nthreads = arg_intn(NULL, "threads", "<n>", 0, 1, "Number of threads"),
      a_density = arg_intn(NULL, "density", "<n>", 0, 1, "Density method, one of {0,1,2}"),

      a_maxiter =
          arg_intn("I", "maxiter", "<num>", 0, 1, "Maximum number of iterations for the k-algo"),
      a_costf =
          arg_intn(NULL, "costf", "<num>", 0, 1, "cost function {1=cond,2=inv(default),meanw}"),

      informat = arg_str0(NULL, "format", "<ascii|binary>", "Input format: ascii or binary"),

#ifdef USE_IGRAPH
      a_algo = arg_str0(NULL, "algo", "<algorithm>",
                        "One of {k,m,fastg,edgeb,leigenv,louvain,walktrap}. (default = m)"),
#else
      a_algo = arg_str0("A", "algo", "<algorithm>", "One of {k,m}. (default = m)"),
#endif
      a_graphtype =
          arg_str0(NULL, "type", "<similarity|distance>", "Are weights distances or similarities"),

      a_scale = arg_str0(NULL, "scale", "<no|yes>", "Scale weights automatically (default:yes)"),
      partfn = arg_filen(NULL, "evalpart", "<file>", 0, 1, "Evaluate partition"),

      a_saveiterparts = arg_intn(NULL, "savparts", "<num>", 0, 1,
                                 "(debug) Save intermediate result partitions. "
                                 "1 = once for each run of k-algo. "
                                 "2 = for every k-algo iteration"),
      a_grow_factor_start =
          arg_dbln(NULL, "gstart", "<num>", 0, 1, "(experimental) grow factor range start"),
      a_grow_factor_end =
          arg_dbln(NULL, "gend", "<num>", 0, 1, "(experimental) grow factor range end"),
      end = arg_end(20),
  };

  printf("Graph clustering. M-algo and K-algo (v. 0.1).\n");
  int ok = 1;
  int nerrors = arg_parse(argc, argv, argtable);
  if (nerrors > 0) {
    ok = 0;
    printf("Unable to parse command line\n");
  }

  if (a_repeats->count > 0) {
    repeats = a_repeats->ival[0];
  }

  if (a_maxiter->count > 0) {
    g_opt.max_iter = a_maxiter->ival[0];
  }

  if (a_maxiter->count > 0) {
    g_opt.max_iter = a_maxiter->ival[0];
  }

  if (a_clusters->count > 0) {
    clusters = a_clusters->ival[0];
  }
  if (a_density->count > 0) {
    g_opt.density_method = a_density->ival[0];
  }

  if (a_write_header->count > 0) {
    output_write_header = a_write_header->ival[0];
  }

  if (a_costf->count > 0) {
    g_opt.costf = a_costf->ival[0];
  }
  if (a_saveiterparts->count > 0) {
    g_opt.debug_save_intermediate_part = 1;
  }

  // These cost functions we minimize:
  if (g_opt.costf <= 3) {
    g_opt.minimize = 1;
    g_opt.costmultip = -1.0;
  } else {
    g_opt.minimize = 0;
    g_opt.costmultip = 1.0;
  }

  if (a_grow_factor->count > 0) {
    g_opt.grow_factor = a_grow_factor->dval[0];
    g_opt.grow_factor_start = g_opt.grow_factor;
    g_opt.grow_factor_end = g_opt.grow_factor;
  }

  if (a_grow_factor_start->count > 0) {
    g_opt.grow_factor_start = a_grow_factor_start->dval[0];
  }
  if (a_grow_factor_end->count > 0) {
    g_opt.grow_factor_end = a_grow_factor_end->dval[0];
  }

  if (a_verbose->count > 0) {
    g_opt.verbose = a_verbose->dval[0];
  }

  if (rngSeed->count > 0) {
    g_opt.seed = rngSeed->ival[0];
  } else {
    g_opt.seed = time(NULL);
  }
  printf("Set RNG seed: %d\n", g_opt.seed);
  srand(g_opt.seed);

  if (a_graphtype->count > 0) {
    if (strcmp(a_graphtype->sval[0], "similarity") == 0) {
      graph_type = SIMILARITY;
    }
    if (strcmp(a_graphtype->sval[0], "distance") == 0) {
      graph_type = DISTANCE;
    }
  }

  g_opt.scale = 1;
  if (strcmp(a_scale->sval[0], "no") == 0) {
    g_opt.scale = 0;
  }

  if (infn->count > 0) {

  } else {
    printf("No input filename given\n");
    ok = 0;
  }
  // printf("INFN:=%s\n", infn->filename[0]);

  if (help->count > 0 || ok == 0) {
    printf("\ngclus");
    arg_print_syntax(stdout, argtable, "\n");
    arg_print_glossary(stdout, argtable, "  %-25s %s\n");
    // printf("ok=%d\n", ok);
    return 0;
  }

  printf("infn='%s' graph_type=%d clusters=%d repeats=%d grow_factor=%f costf=%d\n",
         infn->filename[0], graph_type, clusters, repeats, g_opt.grow_factor, g_opt.costf);

  nnGraph *graph = read_ascii_graphf(infn->filename[0]);
  if (g_opt.scale == 1) {
    if (graph_type == DISTANCE) {
      scale_weights(graph, 1);
    } else {
      scale_weights(graph, 2);
    }
  } else {
    // No scale
    scale_weights(graph, -1);
  }

  if (partfn->count > 0) {
    printf("Evaluating partition in file %s\n", partfn->filename[0]);
    Clustering *clu = read_partition(partfn->filename[0], graph);
    printf("DONE\n");
    double cond, cost;
    g_opt.costf = 1;
    cost = costf(graph, clu, &cond) * -1;
    printf("cost1=%f\n", cost);
    g_opt.costf = 2;
    cost = costf(graph, clu, &cond) * -1;
    printf("cost2=%f\n", cost);
    g_opt.costf = 4;
    cost = costf(graph, clu, &cond);
    printf("cost4=%f\n", cost);

    // double costf(nnGraph *graph, Clustering *clu, double *r_conductance) {
    return 0;
  }

  init_Clustering(&clu, graph->size, clusters);

  g_timer.tick();

  // Default option
  if (a_algo->count == 0 || strcmp(a_algo->sval[0], "m") == 0) {
    m_algo(graph, clu, repeats, clusters);
  }

  if (strcmp(a_algo->sval[0], "k") == 0) {
    m_algo(graph, clu, 0, clusters);
  }

#ifdef USE_IGRAPH
  if (a_algo->count > 0) {
    if (strcmp(a_algo->sval[0], "fastg") == 0) {
      algo_fast_greedy(graph, clu);
    } else if (strcmp(a_algo->sval[0], "edgeb") == 0) {
      algo_edge_betweenness(graph, clu);
    } else if (strcmp(a_algo->sval[0], "leigenv") == 0) {
      algo_leading_eigenv(graph, clu);
    } else if (strcmp(a_algo->sval[0], "louvain") == 0) {
      algo_louvain(graph, clu);
    } else if (strcmp(a_algo->sval[0], "walktrap") == 0) {
      algo_walktrap(graph, clu);
    }
  }
#endif

  printf("Clu sizes: ");
  int maxsize = 0;
  for (int i = 0; i < clu->K; i++) {
    printf("%d(%f) ", clu->clusize[i], -clu->costs[i]);
    if (maxsize < clu->clusize[i]) {
      maxsize = clu->clusize[i];
    }
  }
  printf(" min,max sizes: [%d,%d]", clu->min_part_size, clu->max_part_size);
  printf("\n");

  if (outfn->count > 0) {
    printf("Writing results to file %s\n", outfn->filename[0]);
    FILE *fp = fopen(outfn->filename[0], "w");
    if (output_write_header) {
      fprintf(fp, "VQ PARTITIONING 2.0\n");
      fprintf(fp, "%d\n%d\n", clu->K, clu->N);
      fprintf(fp, "-------------------------------------\n");
    }
    for (int i = 0; i < clu->N; i++) {
      clu->part[i] += 1;
    }
    write_ints_to_fp(fp, clu->part, clu->N);
    fclose(fp);
  }
}
