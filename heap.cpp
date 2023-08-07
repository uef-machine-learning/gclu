// HEAP.CPP
// Based on implementation by Timo Kaukoranta
// in https://github.com/uef-machine-learning/CBModules

#define ProgName "HEAP"
#define VersionNumber "Version 0.02"
#define LastUpdated "27.3.98"

/* ----------------------------------------------------------------- */

#include <assert.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

// #include "memctrl.h"
// #include "owntypes.h"

#define dealloc(V) free(V)
#define alloc(V) malloc(V)

// typedef enum { NO = 0, YES = 1 } YESNO;

class Heap {
public:
  void init(int _allocatedsize);
  // void free();
  Heap();
  virtual ~Heap();
  void update(int target);
  void *remove(int heapindex);
  void insert(void *elem, int *index);
  void *get_top();
  int isHeap();
  int getSize();

  // virtual int compare(void *a, void *b);
  void Heap_siftdown(int target);
  void Heap_siftup(int target);
  // void Heap_moveelem(int target, void* elem, int* index);

  void **data;
  int **index;
  int size;
  int allocatedsize;

  virtual int compare(void *_a, void *_b) {
    // int CmpDistP(void *a, void *b, void *info) {
    // double diff = (((gNode *)a)->nearest_dist - ((gNode *)b)->nearest_dist);
    int a = *((int *)_a);
    int b = *((int *)_b);
    printf("COMPARE\n");

    if (a > b) {
      return 1;
    } else if (a < b) {
      return -1;
    } else {
      return 0;
    }
  }
};

/* ----------------------------------------------------------------- */

void Heap::init(int _allocatedsize) {
  // printf("Heap::init()\n");
  assert(0 <= _allocatedsize);
  allocatedsize = _allocatedsize + 1;
  data = (void **)alloc(allocatedsize * sizeof(void *));
  index = (int **)alloc(allocatedsize * sizeof(int *));
  // compare = compf;
  size = 0;
}

/* ----------------------------------------------------------------- */

Heap::Heap() {
  // printf("Heap::Heap()\n");
  init(10);
}

Heap::~Heap() {
  // printf("Heap::~Heap()\n");
  dealloc(data);
  dealloc(index);
}

void *Heap::get_top() { return data[1]; }

/* ----------------------------------------------------------------- */

// int Heap::compare(void *_a, void *_b) {
// // int CmpDistP(void *a, void *b, void *info) {
// // double diff = (((gNode *)a)->nearest_dist - ((gNode *)b)->nearest_dist);
// int a = *((int *)_a);
// int b = *((int *)_b);
// printf("COMPARE\n");

// if (a > b) {
// return 1;
// } else if (a < b) {
// return -1;
// } else {
// return 0;
// }
// }
// down, as in to having a smaller index.
// towards the root/top of the heap.
void Heap::Heap_siftdown(int target) {
  int i;
  bool flag = true;
  int *tmpi;
  void *tmpd;

  while (flag && target > 1) {
    i = target / 2;
    // if( compare(data[target], data[i], info) >= 0 )
    if (compare(data[target], data[i]) >= 0) {
      flag = false;
    } else {
      tmpd = data[i];
      data[i] = data[target];
      data[target] = tmpd;

      tmpi = index[i];
      index[i] = index[target];
      index[target] = tmpi;

      *(index[i]) = i;
      *(index[target]) = target;

      target = i;
    }
  }
}

void Heap::insert(void *elem, int *_index) {
  int **tmpi;
  void **tmpd;

  if (_index == NULL) {
    _index = (int *)malloc(sizeof(int));
  }

  if (allocatedsize <= size + 2) {
    // printf("Heap:: Allocating more memory (%d)!\n",allocatedsize);
    tmpd = data;
    tmpi = index;
    allocatedsize = allocatedsize * 2 + 1;

    data = (void **) malloc(allocatedsize * sizeof(void *));
    index = (int **) malloc(allocatedsize * sizeof(int *));
    // int ret1,ret2;
    void *ret1, *ret2;
    ret1 = memcpy(data, tmpd, (size + 1) * sizeof(void *));
    ret2 = memcpy(index, tmpi, (size + 1) * sizeof(int *));

    if (!ret1 || !ret2) {
      printf("memcpy fail\n");
    }

    if (!tmpd || !tmpi) {
      printf("memcpy fail\n");
    }

    free(tmpd);
    free(tmpi);
  }

  size++; /* A new element */
  data[size] = elem;
  index[size] = _index;
  *(index[size]) = size;

  Heap_siftdown(size);
}

/* ----------------------------------------------------------------- */

void Heap::Heap_siftup(int target) {
  int j;
  void *tmpd;
  int *tmpi;

  while (2 * target <= size) {
    j = 2 * target;
    if (j < size) {
      /* Select smaller of two sons. */
      // if( compare(data[j], data[j+1], info) > 0 )
      if (compare(data[j], data[j + 1]) > 0) {
        j++;
      }
    }
    // if( compare(data[target], data[j], info) > 0 )
    if (compare(data[target], data[j]) > 0) {
      tmpd = data[j];
      data[j] = data[target];
      data[target] = tmpd;

      tmpi = index[j];
      index[j] = index[target];
      index[target] = tmpi;

      *(index[j]) = j;
      *(index[target]) = target;

      target = j;
    } else {
      target = size + 1;
    }
  }
}

/* ----------------------------------------------------------------- */

int Heap::getSize() {
return size;
}


void Heap::update(int target) {

  // if (target > 1 && compare(data[target / 2], data[target], info) > 0)
  if (target > 1 && compare(data[target / 2], data[target]) > 0) {
    Heap_siftdown(target);
  } else {
    Heap_siftup(target);
  }
}

void *Heap::remove(int heapindex) {
  void *elem;

  // assert(!Heap_isempty(H));
  assert(1 <= heapindex);
  assert(heapindex <= size);

  elem = data[heapindex];
  *(index[heapindex]) = 0; /* Outside of 1..size */

  data[heapindex] = data[size];
  index[heapindex] = index[size];
  *(index[heapindex]) = heapindex;

  size--;

  update(heapindex);
  // Heap_siftup(heapindex);

  return (elem);
}

  //check if Data is in heap order
  int Heap::isHeap() {
    int isheap = 1;
    if (size <= 1) {
      return 1;
    }
    for (int i2 = 1; i2 < size; i2++) {
      // printf("%d: %d %f\n", i2, ((gItem *)(data[i2]))->id, tmpdist);
      int parent = i2 / 2;
      if (i2 > 1 && compare(data[parent], data[i2]) > 0) { // smaller than parent
        isheap = 0;
        printf("ERROR, smaller than parent: %d\n", i2 );
      }
    }

    // assert(notheap == 0);
    // printf("End sanity test\n");
    return isheap;
  }


#ifdef DISABLED00

/* ----------------------------------------------------------------- */

/* ----------------------------------------------------------------- */

void *Heap_remove2(HEAP *H, int heapindex, void *info) {
  void *elem;

  assert(!Heap_isempty(H));
  assert(1 <= heapindex);
  assert(heapindex <= size);

  elem = data[heapindex];
  *(index[heapindex]) = -3; /* Outside of 1..size */

  data[heapindex] = data[size];
  index[heapindex] = index[size];
  *(index[heapindex]) = heapindex;

  size--;

  Heap_siftup(H, heapindex, info);

  return (elem);
}

void *Heap_remove(HEAP *H, int heapindex, void *info) {
  void *elem;

  assert(!Heap_isempty(H));
  assert(1 <= heapindex);
  assert(heapindex <= size);

  elem = data[heapindex];
  *(index[heapindex]) = 0; /* Outside of 1..size */

  data[heapindex] = data[size];
  index[heapindex] = index[size];
  *(index[heapindex]) = heapindex;

  size--;

  Heap_siftup(H, heapindex, info);

  return (elem);
}

/* ----------------------------------------------------------------- */

void Heap_moveelem(HEAP *H, int target, void *elem, int *index) {
  data[target] = elem;
  index[target] = index;
  *(index[target]) = target;
}
#endif

class doubleHeap : public Heap {
public:
  int compare(void *_a, void *_b) {
    double a = *((double *)_a);
    double b = *((double *)_b);

    if (a > b) {
      return 1;
    } else if (a < b) {
      return -1;
    } else {
      return 0;
    }
  }
};

void heap_test() {
  doubleHeap *h2 = new doubleHeap();
  h2->init(20);

  double a00 = 0.1;
  double a01 = 499.88;
  double a02 = 5.0003;
  double a03 = 1.0002;
  int a02ind;
  int a00ind;
  h2->insert((void *)&a00, &a00ind);
  h2->insert((void *)&a01, (int *)malloc(sizeof(int)));
  h2->insert((void *)&a02, &a02ind);
  h2->insert((void *)&a03, (int *)malloc(sizeof(int)));
  double top = *((double *)(h2->get_top()));
  printf("TOP of heap:%f (size:%d)\n", top,h2->size);
  
  a00=9.0009;
  h2->update(a00ind);
  
  a02=0.0009;
  h2->update(a02ind);
  top = *((double *)(h2->get_top()));
  printf("TOP of heap:%f (size:%d)\n", top,h2->size);
  top = *((double *)(h2->remove(1)));
  top = *((double *)(h2->get_top()));
  printf("TOP of heap:%f (size:%d)\n", top,h2->size);
  
    top = *((double *)(h2->remove(1)));
  top = *((double *)(h2->get_top()));
  printf("TOP of heap:%f (size:%d)\n", top,h2->size);

  top = *((double *)(h2->remove(1)));
  top = *((double *)(h2->get_top()));
  printf("TOP of heap:%f (size:%d)\n", top,h2->size);

  // delete h2;
}

