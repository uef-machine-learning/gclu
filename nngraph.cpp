

/*void addNode();*/
/*void removeNode();*/

/*void addNeighbor();*/
/*void removeNeighbor();*/
/*void hasNeighbor();*/

nnGraph *read_ascii_graphf(const char *fname) {

  int N = 0;
  float buf;
  FILE *fp;
  int max_chars = 100000;
  char line[max_chars + 1];
  char *pbuf;
  int i_elem = 0;
  int dim = 0;
  int id;
  int numlinks;
  int links[30000];
  float weights[30000];

  printf("Reading ascii graph dataset from file %s\n", fname);
  fp = fopen(fname, "r");
  if (!fp) {
    terminal_error("File does not exist\n");
  }

  N = count_lines(fp);
  printf("lines=%d\n", N);
  nnGraph *graph = init_nnGraph(N);

  // Get number of elements
  char *ok = fgets(line, max_chars, fp);
  if (ok == NULL) {
    terminal_error("");
  }
  pbuf = line;
  // TODO: Delete
  for (i_elem = 0;; i_elem++) {
    if (*pbuf == '\n')
      break;
    buf = strtof(pbuf, &pbuf);
    // printf(" %f", buf);
  }
  dim = i_elem;
  printf("\nnum_vectors=%d\n", N);

  fseek(fp, 0L, SEEK_SET);
  for (int i_vector = 0; i_vector < N; i_vector++) {
    // for (int i_vector = 0; i_vector < 10; i_vector++) {
    char *ok = fgets(line, max_chars, fp);
    pbuf = line;
    if (ok == NULL) {
      terminal_error("premature end of file");
    }

    id = (int)strtof(pbuf, &pbuf);
    numlinks = (int)strtof(pbuf, &pbuf);
    // printf("id=%d nl=%d ", id, numlinks);

    // printf("links: ");
    for (int i = 0; i < numlinks; i++) {
      if (*pbuf == '\n') {
        terminal_error("Got too few elements");
      }

      buf = strtof(pbuf, &pbuf);
      // printf("%f ", buf);
      links[i] = (int)buf;
    }

    // printf("weights: ");
    for (int i = 0; i < numlinks; i++) {
      if (*pbuf == '\n') {
        terminal_error("Got too few elements");
      }

      buf = strtof(pbuf, &pbuf);
      weights[i] = buf;
      // printf("%f ", buf);
    }
    // printf("\n", buf);

    for (int i = 0; i < numlinks; i++) {
      // printf("add link: [%d] = (%f)> [%d]\n", id, weights[i], links[i]);

      if (!nng_has_neighbor(graph, id, links[i])) {
        nng_add_neighbor(graph, id, links[i], weights[i]);
      }

      if (!nng_has_neighbor(graph, links[i], id)) {
        nng_add_neighbor(graph, links[i], id, weights[i]);
      }
    }
    // printf("\n");
  }
  return graph;
}

double get_max_weight(nnGraph *graph) {
  gNode *nodeA;
  gItem *gi;
  double maxweight = -DBL_MAX;
  for (int i = 0; i < graph->size; i++) {
    nodeA = &graph->nodes[i];

    for (auto gi : *(nodeA->nset)) {

      if (gi->dist > maxweight) {
        maxweight = gi->dist;
      }
    }
  }
  return maxweight;
}

nnGraph *init_nnGraph(int numNodes) {
  int maxNeighbors = 200;
  nnGraph *g = (nnGraph *)malloc(sizeof(nnGraph));
  g->nodes = (gNode *)malloc(sizeof(gNode) * numNodes);

  auto cmp = [](gItem *a, gItem *b) { return a->id < b->id; };
  for (int i = 0; i < numNodes; i++) {
    g->nodes[i].size = 0;
    g->nodes[i].weight = 1;
    g->nodes[i].outdated = 0;
    g->nodes[i].internalSum = 0.0;

    g->nodes[i].id = i;
    g->nodes[i].visited = -1;
    g->nodes[i].maxNeighbors = maxNeighbors;

    g->nodes[i].nset = new std::set<gItem *, custom_compare>;
    g->nodes[i].nearesth = new giHeap();

    g->nodes[i].stash = new std::vector<int>;
    g->nodes[i].stash->push_back(i);
    // gItem*
  }
  g->size = numNodes;
  return g;
}

void dealloc_nnGraph(nnGraph *g) {
  // nnGraph *g = (nnGraph *)malloc(sizeof(nnGraph));
  // g->nodes = (gNode *)malloc(sizeof(gNode) * numNodes);

  for (int i = 0; i < g->size; i++) {
    for (auto gi : *(g->nodes[i].nset)) {
      free(gi);
    }
    delete g->nodes[i].stash;
    delete g->nodes[i].nearesth;
    delete g->nodes[i].nset;
  }
  free(g->nodes);
  free(g);
}

gItem *nng_add_mutual_neighbor2(nnGraph *g, int p1, int p2, float dist) {
  gNode *node = &g->nodes[p1];
  assert(p1 != p2);

  gItem *gi_p1 = (gItem *)malloc(sizeof(gItem));
  gItem *gi_p2 = (gItem *)malloc(sizeof(gItem));

  gi_p1->cost = -1;
  gi_p1->visited = -1;
  gi_p1->heapp = -1;

  gi_p2->cost = -1;
  gi_p2->visited = -1;
  gi_p2->heapp = -1;

  float cost;
  int visited;
  gItem *pair;
  int heapp; // Pointer to heap position

  gi_p1->visited = 0;
  gi_p2->visited = 0;

  gi_p1->id = p2;
  gi_p1->dist = dist;

  gi_p2->id = p1;
  gi_p2->dist = dist;

  gi_p2->pair = gi_p1;
  gi_p1->pair = gi_p2;

  auto iter1 = g->nodes[p1].nset->insert(gi_p1);
  auto iter2 = g->nodes[p2].nset->insert(gi_p2);

  // Insertion took place
  if (iter1.second) {
    gi_p1->iterO = iter2.first;
    gi_p2->iterO = iter1.first;

    g->nodes[p1].nearesth->insert((void *)gi_p1, &(gi_p1->heapp));
    g->nodes[p2].nearesth->insert((void *)gi_p2, &(gi_p2->heapp));

  } else {
    free(gi_p1);
    free(gi_p2);
  }

  return gi_p1;
}

void nng_add_neighbor_safe(nnGraph *graph, int p1, int p2, float dist) {
  gNode *node = &graph->nodes[p1];

  if (nng_has_neighbor(graph, p1, p2)) // Already exists
  {
    return;
  }

  gItem *gi = (gItem *)malloc(sizeof(gItem));
  gi->id = p2;
  gi->dist = dist;
  // printf("dist2=%f\n",gi->dist);

  graph->nodes[p1].nset->insert(gi);

  node->size++;
}

void nng_add_neighbor(nnGraph *g, int p1, int p2, float dist) {
  gNode *node = &g->nodes[p1];

  gItem *gi = (gItem *)malloc(sizeof(gItem));
  gi->id = p2;
  gi->dist = dist;
  // printf("dist2=%f\n",gi->dist);

  g->nodes[p1].nset->insert(gi);

  node->size++;
}

void nng_add_mutual_neighbor(nnGraph *g, int p1, int p2, float dist) {
  nng_add_neighbor(g, p1, p2, dist);
  nng_add_neighbor(g, p2, p1, dist);
}

void nng_remove_neighbor(nnGraph *g, int p1, int p2) {

  gItem *gi = (gItem *)malloc(sizeof(gItem));
  gi->id = p2;
  auto it = g->nodes[p1].nset->find(gi);
  if (it != g->nodes[p1].nset->end()) {
    g->nodes[p1].nset->erase(it);
    free(*it);
  }

  // free(gi); //TODO
  /*node->size--;*/ // TODO:??
}

int nng_has_neighbor(nnGraph *g, int p1, int p2) {

  gItem *gi = (gItem *)malloc(sizeof(gItem));
  gi->id = p2;
  auto it = g->nodes[p1].nset->find(gi);
  int ret = 0;
  if (it != g->nodes[p1].nset->end()) {
    ret = 1;
  }
  free(gi);
  return ret;
}

gItem *nng_get_neighbor(nnGraph *g, int p1, int p2) {

  gItem *gi = (gItem *)malloc(sizeof(gItem));
  gi->id = p2;
  auto it = g->nodes[p1].nset->find(gi);
  gItem *ret = NULL;
  if (it != g->nodes[p1].nset->end()) {
    ret = *it;
  }
  // std::cout << "decltype(it) is " << type_name<decltype(it)>() << '\n';
  // std::cout << "decltype(it) is " << type_name<decltype(*it)>() << '\n';
  free(gi);
  // return gi;
  return ret;
}

// https://thispointer.com/how-to-access-element-by-index-in-a-set-c/
template <typename T> std::pair<T, bool> getNthElement(std::set<T> &searchSet, int n) {
  std::pair<T, bool> result;
  if (searchSet.size() > n) {
    result.first = *(std::next(searchSet.begin(), n));
    result.second = true;
  } else
    result.second = false;
  return result;
}

gItem *nng_get_neighbor2(nnGraph *g, int p1, int idx) {
  // std::pair<gItem *, bool> result = getNthElement(*(g->nodes[p1].nset), 3);

  gItem *gi = *(std::next(g->nodes[p1].nset->begin(), idx));
  return gi;
}

int nng_num_neighbors(nnGraph *g, int p1) { return g->nodes[p1].nset->size(); }

void write_nngraph_to_file(nnGraph *g, const char *fn) {
  FILE *fp;
  fp = fopen(fn, "w");
  gNode *node;
  for (int i = 0; i < g->size; i++) {
    node = &g->nodes[i];
    int num_neigh = node->nset->size();
    fprintf(fp, "%d %d", node->id, num_neigh);
    for (int j = 0; j < num_neigh; j++) {
      // gItem *gi = (gItem *)ll_get_item(node->neighbors, j);
      gItem *gi = nng_get_neighbor2(g, i, j);

      fprintf(fp, " %d", gi->id);
    }
    for (int j = 0; j < num_neigh; j++) {
      // gItem *gi = (gItem *)ll_get_item(node->neighbors, j);
      gItem *gi = nng_get_neighbor2(g, i, j);
      fprintf(fp, " %f", gi->dist);
    }

    fprintf(fp, "\n");
  }
  fclose(fp);
}

void write_nngraph_to_file_old(nnGraph *g, const char *fn) {
  FILE *fp;
  fp = fopen(fn, "w");
  gNode *node;
  for (int i = 0; i < g->size; i++) {
    node = &g->nodes[i];
    fprintf(fp, "%d %d", node->id, node->neighbors->size);
    for (int j = 0; j < node->neighbors->size; j++) {
      // gItem *gi = (gItem *)ll_get_item(node->neighbors, j);
      gItem *gi = nng_get_neighbor2(g, i, j);

      fprintf(fp, " %d", gi->id);
    }
    for (int j = 0; j < node->neighbors->size; j++) {
      // gItem *gi = (gItem *)ll_get_item(node->neighbors, j);
      gItem *gi = nng_get_neighbor2(g, i, j);
      fprintf(fp, " %f", gi->dist);
    }

    fprintf(fp, "\n");
  }
  fclose(fp);
}

#ifdef LDATASET
gItem *find_greedy_path(DataSet *data, nnGraph *g, int source, int target) {
  gNode *node;
  node = &g->nodes[source];

  gItem *ti = ll_get_node_if_exist(g->nodes[source].neighbors, target);
  float d_ab;
  if (ti == NULL) {
    d_ab = distance(data, source, target);
    /*terminal_error("ll_get_node_if_exist");*/
  } else {
    d_ab = ti->dist;
  }

  for (int j = 0; j < node->neighbors->size; j++) {
    gItem *gi = (gItem *)ll_get_item(node->neighbors, j);
    if (gi->id == target) {
      continue;
    }
    gItem *neighbor = ll_get_node_if_exist(g->nodes[gi->id].neighbors, target);
    if (neighbor == NULL) {
      /*printf("nn:NULL ");*/
    }
    if (neighbor != NULL) {
      /*printf("nn:%d %f, ",neighbor->id,neighbor->dist);*/
      if (d_ab > neighbor->dist) {
        return neighbor;
      }
    }
  }

  return NULL;
}

gItem *find_greedy_path2(DataSet *data, nnGraph *g, int source, int target) {
  gNode *node;
  node = &g->nodes[source];

  gItem *ti = ll_get_node_if_exist(g->nodes[source].neighbors, target);
  float d_ab;
  if (ti == NULL) {
    d_ab = distance(data, source, target);
    /*terminal_error("ll_get_node_if_exist");*/
  } else {
    d_ab = ti->dist;
  }

  float prevdist = d_ab;
  for (int j = 0; j < node->neighbors->size; j++) {
    gItem *gi = (gItem *)ll_get_item(node->neighbors, j);
    if (gi->id == target) {
      continue;
    }

    /*float cb = distance(data,gi->id,target);*/
    /*if(cb < prevdist) {*/
    /*return gi;*/

    /*}*/

    gItem *neighbor = ll_get_node_if_exist(g->nodes[gi->id].neighbors, target);
    if (neighbor == NULL) {
      /*printf("nn:NULL ");*/
    }
    if (neighbor != NULL) {
      /*printf("nn:%d %f, ",neighbor->id,neighbor->dist);*/
      if (d_ab > neighbor->dist) {
        return neighbor;
      }
    }
  }

  return NULL;
}
#endif

void test_nn_graph() {
  nnGraph *g = init_nnGraph(100);
  nng_add_neighbor_safe(g, 0, 99, 2.0);
  nng_add_neighbor_safe(g, 0, 90, 2.0);
  nng_add_neighbor_safe(g, 0, 10, 2.0);
  nng_add_neighbor_safe(g, 0, 20, 2.0);
  nng_add_neighbor_safe(g, 0, 7, 2.0);
  nng_add_neighbor_safe(g, 0, 90, 2.0);
  nng_add_neighbor_safe(g, 0, 70, 2.0);
  nng_add_neighbor_safe(g, 0, 90, 2.0);
  nng_add_neighbor_safe(g, 0, 77, 2.0);
  nng_add_neighbor_safe(g, 0, 90, 2.0);
  nng_add_neighbor_safe(g, 3, 90, 2.0);
  nng_add_neighbor_safe(g, 3, 70, 2.0);
  nng_add_neighbor_safe(g, 3, 90, 2.0);
  nng_add_neighbor_safe(g, 3, 77, 2.0);

  // nng_remove_neighbor(g,0,90);
  nng_remove_neighbor(g, 0, 99);
  nng_remove_neighbor(g, 0, 77);
  nng_remove_neighbor(g, 0, 10);
  // for(int i=0;i<6;i++) {
  printf("n=%d\n", nng_num_neighbors(g, 0));
  fflush(stdout);
  for (int i = 0; i < nng_num_neighbors(g, 0); i++) {
    // printf(" %d", nng_get_neighbor(g, 0, i));
  }

  printf(" | ");
  for (int i = 0; i < nng_num_neighbors(g, 0); i++) {
    gItem *gi = nng_get_neighbor2(g, 0, i);
    printf(" %d", gi->id);
  }

  printf(" | ");
  for (auto gi : *(g->nodes[0].nset)) {
    printf(" %d", gi->id);
    // use f here
  }

  printf("\n");
}

nnGraph *read_ascii_graphf2(const char *fname, int clist, std::map<std::string, int> **cdmap) {

  std::string line;
  std::string delim = " ";
  std::ifstream infile(fname);

#ifdef LDATASET
  DataSet *sd = (DataSet *)malloc(sizeof(DataSet));
#endif
  int numLines = 0;
  int numItems = 0;
  int buf[1000000];

  // std::map<std::string, int> m;
  // std::map<std::string, int> code_to_id;
  std::map<std::string, int> *code_to_id = new std::map<std::string, int>;

  std::map<int, std::string> id_to_code;
  std::stringstream ss(line);
  std::string item;
  int id = -1;
  int numlines = 0;
  while (std::getline(infile, line)) {
    numlines++;
  }
  infile.clear();
  infile.seekg(0, ios::beg);
  printf("Num lines: %d\n", numlines);
  int N = numlines - 1;

  int show_sets = 1;

  nnGraph *graph = init_nnGraph(N);

  // Read first line which contains the names/codes of the set items
  // Convert to int id's
  safeGetline(infile, line);
  {
    std::stringstream ss(line);
    std::string item;
    int id = -1;
    int setSize = 0;
    while (std::getline(ss, item, ' ')) {
      if (item.size() < 1) {
        continue;
      }
      if ((*code_to_id).find(item) == (*code_to_id).end()) {
        (*code_to_id)[item] = numItems;
        id_to_code[numItems] = item;
        numItems++;
      } else {
      }
      id = (*code_to_id)[item];
      // if (show_sets) {
      // cout << "  word: '" << item << "' mapping:" << id << endl;
      // }
      buf[setSize] = id;
      setSize++;
    }
  }

  while (safeGetline(infile, line)) {
    vector<int> links;
    vector<float> weights;

    if (line.size() <= 1) {
      continue;
    }

    std::stringstream ss(line);
    std::string item;
    std::getline(ss, item, ' ');
    // if (item.size() < 1) {
    // continue;
    // }

    int id = stoi(item);

    std::getline(ss, item, ' ');

    int numlinks = stoi(item);
    // printf("id=%d numlinks=%d ", id, numlinks);
    for (int i = 0; i < numlinks; i++) {
      std::getline(ss, item, ' ');
      int neighborid = stoi(item);
      // printf("%d ", neighborid);
      // links[i] = neighborid;
      links.push_back(neighborid);
    }
    for (int i = 0; i < numlinks; i++) {
      std::getline(ss, item, ' ');
      float weight = stof(item);
      // weights[i] = weight;
      weights.push_back(weight);
      // printf("%f ", weight);
    }
    // printf("numw=%d\n",weights.size());

    for (int i = 0; i < numlinks; i++) {
      // printf("add link: [%d] = (%f)> [%d]\n", id, weights[i], links[i]);
      if (!nng_has_neighbor(graph, id, links[i])) {
        nng_add_neighbor(graph, id, links[i], weights[i]);
      }

      if (!nng_has_neighbor(graph, links[i], id)) {
        nng_add_neighbor(graph, links[i], id, weights[i]);
      }
    }
  }

  // printf("##F43:%d\n",(*cdmap)[string("F43")]);
  printf("#Z#F43:%d\n", (*code_to_id)[string("F43")]);
  if (cdmap != NULL) {
    // cdmap = &code_to_id;
    // cdmap = code_to_id;
    *cdmap = code_to_id;
  }
  // printf("##F43:%d\n",(**cdmap)[string("F43")]);
  return graph;
}
