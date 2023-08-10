
int g_iter = 0;

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
    
    
    for (auto gi : *(node->nset)) {
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


    for (auto gi : *(seed->nset)) {
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


    for (auto gi : *(node->nset)) {
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
    
    for (auto gi : *(nodeA->nset)) {
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
    nodeA->weight_sum = 0.0;
    for (auto gi : *(nodeA->nset)) {
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


  for (auto gi : *(nodeA->nset)) {
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
    // nodeA->weight_sum = 0.0;

    items[i].id = i;
    items[i].density = 0;
    if (g_opt.density_method == 1) {
    
    
    for (auto gi : *(nodeA->nset)) {
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
  clu->cost = -999999.0;
  for (int i_rep = 0; i_rep < n_repeats; i_rep++) {
    init_Clustering(&newclu, graph->size, n_clusters);
    k_algo(graph, newclu);
    printf("REP = %d, cost=%f best=%f\n", i_rep, newclu->cost, clu->cost);

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

