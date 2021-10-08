

void to_igraph(nnGraph *graph, igraph_t *g, igraph_vector_t *weights) {
  igraph_vector_t v, res, reset;
  gNode *nodeA;
  gItem *gi;

  // igraph_t g;
  int num_weights;
  num_weights = 10;
  igraph_vector_init(&v, 0);
  int edge_i = 0;
  for (int i = 0; i < graph->size; i++) {
    nodeA = &graph->nodes[i];
    for (int j = 0; j < nodeA->neighbors->size; j++) {
      gi = (gItem *)ll_get_item(nodeA->neighbors, j);

      // Add edge only once, assume all connections are two way in graph
      if (i < gi->id) {
        igraph_vector_push_back(&v, i);
        igraph_vector_push_back(&v, gi->id);
        // if(gi->dist == 0.0) {
        // printf("gi-dist:%f\n",gi->dist);
        // }
        // igraph_vector_push_back(weights, (int) gi->dist*10000);
        igraph_vector_push_back(weights, gi->dist+0.000001);
        edge_i++;
      }
    }
  }
  igraph_create(g, &v, graph->size /*number of vertices*/, 0 /*0=undirected*/);
}

void algo_fast_greedy(nnGraph *graph, Clustering *clu) {
  igraph_t g;
  igraph_vector_t modularity, weights, membership;

  igraph_matrix_t merges;
  igraph_matrix_init(&merges, 0, 0);
  igraph_vector_init(&modularity, 0);
  igraph_vector_init(&membership, 0);
  igraph_vector_init(&weights, 0);
  to_igraph(graph, &g, &weights);

  igraph_community_fastgreedy(&g, &weights, &merges, &modularity, 0 /*membership*/);

  igraph_community_to_membership(&merges, igraph_vcount(&g), igraph_vcount(&g) - clu->K,
                                 &membership, 0);
  for (int i = 0; i < igraph_vector_size(&membership); i++) {
    clu->part[i] = (int)VECTOR(membership)[i];
  }
  printf("FINAL=1 TIME=%f\n", g_timer.get_time());
}

void algo_walktrap(nnGraph *graph, Clustering *clu) {
  igraph_t g;
  igraph_vector_t modularity, weights, membership;

  igraph_matrix_t merges;
  igraph_matrix_init(&merges, 0, 0);
  igraph_vector_init(&modularity, 0);
  igraph_vector_init(&membership, 0);
  igraph_vector_init(&weights, 0);
  to_igraph(graph, &g, &weights);

  // igraph_community_fastgreedy(&g, &weights, &merges, &modularity, 0 /*membership*/);
  igraph_community_walktrap(&g, &weights, 4, &merges, &modularity, 0 /*membership*/);

  igraph_community_to_membership(&merges, igraph_vcount(&g), igraph_vcount(&g) - clu->K,
                                 &membership, 0);
  for (int i = 0; i < igraph_vector_size(&membership); i++) {
    clu->part[i] = (int)VECTOR(membership)[i];
  }
  printf("FINAL=1 TIME=%f\n", g_timer.get_time());
}


// int igraph_community_walktrap(const igraph_t *graph, 
 // const igraph_vector_t *weights,
 // int steps,
 // igraph_matrix_t *merges,
 // igraph_vector_t *modularity, 
 // igraph_vector_t *membership);
// }

void algo_leading_eigenv(nnGraph *graph, Clustering *clu) {
  igraph_t g;
  igraph_vector_t modularity, weights, membership;

  igraph_matrix_t merges;
  igraph_matrix_init(&merges, 0, 0);
  igraph_vector_init(&modularity, 0);
  igraph_vector_init(&membership, 0);
  igraph_vector_init(&weights, 0);
  to_igraph(graph, &g, &weights);

  // igraph_community_fastgreedy(&g, &weights, &merges, &modularity, 0 /*membership*/);
  igraph_arpack_options_t options;
  igraph_arpack_options_init(&options);

  igraph_community_leading_eigenvector(
      &g, &weights, &merges, &membership, 1, &options, /*modularity=*/0, /*start=*/0,
      /*eigenvalues=*/0, /*eigenvectors=*/0, /*history=*/0, /*callback=*/0, /*callback_extra=*/0);

  // int igraph_community_leading_eigenvector(const igraph_t *graph,
  // const igraph_vector_t *weights,
  // igraph_matrix_t *merges,
  // igraph_vector_t *membership,
  // igraph_integer_t steps,
  // igraph_arpack_options_t *options,
  // igraph_real_t *modularity,
  // igraph_bool_t start,
  // igraph_vector_t *eigenvalues,
  // igraph_vector_ptr_t *eigenvectors,
  // igraph_vector_t *history,
  // igraph_community_leading_eigenvector_callback_t *callback,
  // void *callback_extra);

  // TODO:
  // igraph_le_community_to_membership
  igraph_community_to_membership(&merges, igraph_vcount(&g), igraph_vcount(&g) - clu->K,
                                 &membership, 0);
  for (int i = 0; i < igraph_vector_size(&membership); i++) {
    clu->part[i] = (int)VECTOR(membership)[i];
  }
  printf("FINAL=1 TIME=%f\n", g_timer.get_time());
}

// https://igraph.org/c/doc/igraph-Community.html#idm231968320928
//  VD Blondel, J-L Guillaume, R Lambiotte and E Lefebvre: Fast unfolding of community hierarchies
//  in large networks, J Stat Mech P10008 (2008)
void algo_louvain(nnGraph *graph, Clustering *clu) {
  igraph_t g;
  igraph_vector_t modularity, weights,membership;
  // igraph_realt_t membership;
  igraph_matrix_t memberships;
  igraph_matrix_t merges;

  igraph_matrix_init(&memberships, 0, 0);
  igraph_matrix_init(&merges, 0, 0);
  igraph_vector_init(&modularity, 0);
  igraph_vector_init(&membership, 0);
  igraph_vector_init(&weights, 0);
  to_igraph(graph, &g, &weights);

  igraph_community_multilevel(&g, &weights,1, &membership, &memberships, &modularity);
  
  //TODO:
  // igraph_vector_max(&membership);

  // igraph_community_edge_betweenness(&g, NULL, NULL, &merges, NULL, &modularity, NULL, 0,
  // &weights);

  // int igraph_community_edge_betweenness(const igraph_t *graph,
  // igraph_vector_t *result,
  // igraph_vector_t *edge_betweenness,
  // igraph_matrix_t *merges,
  // igraph_vector_t *bridges,
  // igraph_vector_t *modularity,
  // igraph_vector_t *membership,
  // igraph_bool_t directed,
  // const igraph_vector_t *weights);

  // igraph_community_fastgreedy(&g, &weights, &merges, &modularity, 0 /*membership*/);

  // igraph_community_to_membership(&merges, igraph_vcount(&g), igraph_vcount(&g) - clu->K,
  // &membership, 0);

  printf("\n");
  
  // Choosing clustering with number of clusters closest to given parameter k
  int diff = INT_MAX;
  int difftmp, best_clu, best_clu_k;
  for (int i = 0; i < igraph_matrix_nrow(&memberships); i++) {
    int maxelem = 0;
    for (int j = 0; j < igraph_vcount(&g); j++) {
      int tmp = MATRIX(memberships, i, j);
      if (tmp > maxelem) {
        maxelem = tmp;
      }
    }
    difftmp = abs(clu->K - (maxelem + 1));
    if (difftmp < diff) {
      diff = difftmp;
      best_clu = i;
      best_clu_k = maxelem + 1;
    }

    printf("row=%d clusters=%d\n", i, maxelem + 1);
  }

  printf("Clustering size closest to target, k=%d i=%d diff=%d\n", best_clu_k, best_clu, diff);

  for (int i = 0; i < igraph_vcount(&g); i++) {
    clu->part[i] = MATRIX(memberships, best_clu, i);
  }
  
  // TODO, as option: Choose case of best modularity
  // clu->K = igraph_vector_max(&membership) + 1;
  
  printf("FINAL=1 TIME=%f k=%d res_k=%d k_diff=%d \n", g_timer.get_time(), clu->K, best_clu_k,abs(best_clu_k-clu->K));
  clu->K = best_clu_k;
}

// https://igraph.org/c/doc/igraph-Community.html#idm231968320928
void algo_edge_betweenness(nnGraph *graph, Clustering *clu) {
  igraph_t g;
  igraph_vector_t modularity, weights, membership;

  igraph_matrix_t merges;
  igraph_matrix_init(&merges, 0, 0);
  igraph_vector_init(&modularity, 0);
  igraph_vector_init(&membership, 0);
  igraph_vector_init(&weights, 0);
  to_igraph(graph, &g, &weights);

  igraph_community_edge_betweenness(&g, NULL, NULL, &merges, NULL, &modularity, NULL, 0, &weights);

  igraph_community_to_membership(&merges, igraph_vcount(&g), igraph_vcount(&g) - clu->K,
                                 &membership, 0);
  for (int i = 0; i < igraph_vector_size(&membership); i++) {
    clu->part[i] = (int)VECTOR(membership)[i];
  }
  printf("FINAL=1 TIME=%f\n", g_timer.get_time());
}




