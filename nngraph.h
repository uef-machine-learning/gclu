#ifndef NN_GRAPH_H
#define NN_GRAPH_H

#include "linked_list.h"

typedef struct gNode {
  linkedList *neighbors;
  int numNeighbors;
  int maxNeighbors;
  double weight_sum;
  int size;
  // <id> of point which nearest neighbors this represents
  int id;
} gNode;

typedef struct nnGraph {
  int size;
  /*void * content;*/
  /*linkedListNode* next;*/
  gNode *nodes;
  double total_weight;
} nnGraph;

/*void addNode();*/
/*void removeNode();*/

/*void addNeighbor();*/
/*void removeNeighbor();*/
/*void hasNeighbor();*/

/*void nng_add_neighbor(nnGraph* g, int id, int neighbor) {*/
/*add_neighbor(g->nodes[id]->neighbors, neighbor);*/
/*}*/

nnGraph *init_nnGraph(int numNodes) {
  int maxNeighbors = 200;
  nnGraph *g = (nnGraph *)malloc(sizeof(nnGraph));
  g->nodes = (gNode *)malloc(sizeof(gNode) * numNodes);
  for (int i = 0; i < numNodes; i++) {
    /*g->nodes[i].neighbors = (gItem*) malloc(sizeof(gItem)*maxNeighbors);*/
    g->nodes[i].size = 0;
    g->nodes[i].id = i;
    g->nodes[i].maxNeighbors = maxNeighbors;

    g->nodes[i].neighbors = initLinkedList();
  }
  g->size = numNodes;
  return g;
}

int nng_neighbor_exists(nnGraph *g, int p1, int p2) {
  gNode *node;
  node = &g->nodes[p1];

  for (int j = 0; j < node->neighbors->size; j++) {
    gItem *gi = (gItem *)ll_get_item(node->neighbors, j);
    if (gi->id == p2) {
      return 1;
    }
  }
  return 0;
}

void nng_add_neighbor(nnGraph *g, int p1, int p2, float dist) {
  gNode *node = &g->nodes[p1];
  /*gItem* neighbors = node->neighbors;*/

  gItem *gi = (gItem *)malloc(sizeof(gItem));
  /*gItem* gi = &node->neighbors[node->size];*/
  gi->id = p2;
  gi->dist = dist;

  // ll_remove_node(g->nodes[p1].neighbors, p2);
  ll_add_node(node->neighbors, (void *)gi);
  /*ll_add_node_if_not_exist(node->neighbors,(void*) gi);*/

  node->size++;
}

void nng_add_mutual_neighbor(nnGraph *g, int p1, int p2, float dist) {
  nng_add_neighbor(g, p1, p2, dist);
  nng_add_neighbor(g, p2, p1, dist);
}

void nng_remove_neighbor(nnGraph *g, int p1, int p2) {
  ll_remove_node2(g->nodes[p1].neighbors, p2);
  /*node->size--;*/
}

int nng_has_neighbor(nnGraph *g, int p1, int p2) {
  gItem *gi = NULL;
  gItem *neighbor = ll_get_node_if_exist(g->nodes[p1].neighbors, p2);
  if (neighbor == NULL) {
    return 0;
  }
  return 1;
}

int nng_get_neighbor(nnGraph *g, int p1, int idx) {
  gItem *gi = (gItem *)ll_get_item(g->nodes[p1].neighbors, idx);
  return gi->id;
}

float nng_get_weight(nnGraph *g, int p1, int idx) {
  gItem *gi = (gItem *)ll_get_item(g->nodes[p1].neighbors, idx);
  return gi->dist;
}

int nng_num_neighbors(nnGraph *g, int p1) { return g->nodes[p1].neighbors->size; }

void write_nngraph_to_file(nnGraph *g, const char *fn) {
  FILE *fp;
  fp = fopen(fn, "w");
  gNode *node;
  for (int i = 0; i < g->size; i++) {
    node = &g->nodes[i];
    fprintf(fp, "%d %d", node->id, node->neighbors->size);
    for (int j = 0; j < node->neighbors->size; j++) {
      gItem *gi = (gItem *)ll_get_item(node->neighbors, j);
      fprintf(fp, " %d", gi->id);
    }
    for (int j = 0; j < node->neighbors->size; j++) {
      gItem *gi = (gItem *)ll_get_item(node->neighbors, j);
      fprintf(fp, " %f", gi->dist);
    }

    fprintf(fp, "\n");
  }
  fclose(fp);
}

void debug_dump_graph(nnGraph *g) {
  printf("==============\n");
  printf("n=%d\n", nng_num_neighbors(g, 0));
  fflush(stdout);
  for (int i_node = 0; i_node < g->size; i_node++) {
    printf("[%d] ", i_node);
    for (int i = 0; i < nng_num_neighbors(g, i_node); i++) {
      printf(" %d (%f)", nng_get_neighbor(g, i_node, i), nng_get_weight(g, i_node, i));
    }
    printf("\n");
  }
  printf("==============\n");
}

double get_max_weight(nnGraph *graph) {
  gNode *nodeA;
  gItem *gi;
  double maxweight = -DBL_MAX;
  for (int i = 0; i < graph->size; i++) {
    nodeA = &graph->nodes[i];
    gi = (gItem *)ll_get_item(nodeA->neighbors, nodeA->neighbors->size - 1);
    for (int j = 0; j < nodeA->neighbors->size; j++) {
      gi = (gItem *)ll_get_item(nodeA->neighbors, j);

      if (gi->dist > maxweight) {
        maxweight = gi->dist;
      }
    }
  }
  return maxweight;
}

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
  //TODO: Delete
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
      if (!nng_neighbor_exists(graph, id, links[i])) {
        nng_add_neighbor(graph, id, links[i], weights[i]);
      }

      if (!nng_neighbor_exists(graph, links[i], id)) {
        nng_add_neighbor(graph, links[i], id, weights[i]);
      }
    }
    // printf("\n");
  }
  return graph;
}

void test_nn_graph() {
  nnGraph *g = init_nnGraph(100);
  nng_add_neighbor(g, 0, 99, 2.0);
  nng_add_neighbor(g, 0, 10, 2.0);
  nng_add_neighbor(g, 0, 20, 2.0);
  nng_add_neighbor(g, 0, 7, 2.0);
  nng_add_neighbor(g, 0, 90, 2.0);
  nng_add_neighbor(g, 0, 70, 2.0);
  nng_add_neighbor(g, 0, 77, 2.0);
  // nng_remove_neighbor(g,0,90);
  nng_remove_neighbor(g, 0, 99);
  nng_remove_neighbor(g, 0, 77);
  nng_remove_neighbor(g, 0, 10);
  // for(int i=0;i<6;i++) {
  printf("n=%d\n", nng_num_neighbors(g, 0));
  fflush(stdout);
  for (int i = 0; i < nng_num_neighbors(g, 0); i++) {
    printf(" %d", nng_get_neighbor(g, 0, i));
  }
  printf("\n");
}

#endif
