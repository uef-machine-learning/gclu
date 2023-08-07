#ifndef NN_GRAPH_H
#define NN_GRAPH_H

#include <set>
#include <vector>
#include <vector>
#include <algorithm>
#include <queue>

// #include <stdlib.h>
#include <cstdlib>
#include <cassert>

typedef int BOOL;

struct gItem;
struct gItem {
  int id;
  double dist;
  double cost;
  int visited;
  gItem *pair;
  std::_Rb_tree_const_iterator<gItem *> iterO;
  int heapp; // Pointer to heap position
};

class giHeap : public Heap {
public:

  int compare(void *_a, void *_b) {
    // printf("comp\n");

    gItem *a = (gItem *)_a;
    gItem *b = (gItem *)_b;

    if (a->dist > b->dist) {
      return 1;
    } else if (a->dist < b->dist) {
      return -1;
    } else {
      printf("EQUALS\n");
      return 0;
    }
  }
  int checkSanity() {
    if (size <= 1) {
      return 0;
    }
    int notheap = 0;
    float firstdist = ((gItem *)(data[1]))->dist;
    float mindist = firstdist;
    for (int i2 = 1; i2 <= size; i2++) {
      float tmpdist = ((gItem *)(data[i2]))->dist;
      if (i2 == 1 || tmpdist < mindist) {
        mindist = tmpdist;
      }
      printf("%d: %d %f\n", i2, ((gItem *)(data[i2]))->id, tmpdist);
      int parent = i2 / 2;
      if (i2 > 1 && compare(data[parent], data[i2]) > 0) { // smaller than parent
        notheap = 1;
        printf("ERROR, smaller than parent %d: %d %f\n", i2, ((gItem *)(data[i2]))->id, tmpdist);
      }
    }

    // if(mindist !=firstdist) {printf("HEAPERROR: %f %f\n");}
    assert(notheap == 0);
    assert(mindist == firstdist);

    printf("End sanity test\n");
    return 0;
  }
};

struct custom_compare final {
  bool operator()(gItem *a, gItem *b) const { return a->id < b->id; }
};

typedef struct gNode {
  // List of <size> number of nearest neighbors
  /*kNNItem * items;*/
  /*float max_dist;*/
  /*gItem* neighbors;*/
  linkedList *neighbors;
  int numNeighbors;
  int maxNeighbors;
  int size;
  int weight;
  // <id> of point which nearest neighbors this represents
  int id;
  unsigned char outdated;

  int nearest_id;
  int visited;
  float visited_dist;
  float nearest_dist;
  gItem *nearest;

  int heap_index;
  float *mean;

  float internalSum;

  // std::set<gItem *> *nset;
  // Neighbors as set  for fast check if is neighbor or not
  std::set<gItem *, custom_compare> *nset;
  // HEAP * nset_dist; // Neighbors with nearest on top
  std::vector<int> *stash;
  doubleHeap *h2;
  giHeap *nearesth;
  // floatHeap *h2 = new floatHeap();
  // h2->init(20);

  /*BOOL is_exact;*/
} gNode;

typedef struct nnGraph {
  int size;
  int cur_iter;
  /*void * content;*/
  /*linkedListNode* next;*/
  gNode *nodes;
  DataSet *data;
} nnGraph;

nnGraph *init_nnGraph(int numNodes);
void dealloc_nnGraph(nnGraph *g);
gItem *nng_add_mutual_neighbor2(nnGraph *g, int p1, int p2, float dist);
void nng_add_neighbor(nnGraph *g, int p1, int p2, float dist);
void nng_add_mutual_neighbor(nnGraph *g, int p1, int p2, float dist);
void nng_remove_neighbor(nnGraph *g, int p1, int p2);
int nng_has_neighbor(nnGraph *g, int p1, int p2);
gItem *nng_get_neighbor(nnGraph *g, int p1, int p2);
gItem *nng_get_neighbor2(nnGraph *g, int p1, int idx);
int nng_num_neighbors(nnGraph *g, int p1);
void write_nngraph_to_file(nnGraph *g, const char *fn);
void write_nngraph_to_file_old(nnGraph *g, const char *fn);
gItem *find_greedy_path(DataSet *data, nnGraph *g, int source, int target);
gItem *find_greedy_path2(DataSet *data, nnGraph *g, int source, int target);
void test_nn_graph();
nnGraph *read_ascii_graphf2(const char *fname, int clist, std::map<std::string, int> **cdmap);

#endif
