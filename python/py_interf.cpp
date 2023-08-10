

#include <Python.h>

#include <stdio.h>
#include <limits.h>
#include <cstring>
#include <pthread.h>
// #include "contrib/argtable3.h"

// using namespace std;

#include <execinfo.h>
#include <signal.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdio.h>
#include <iostream>
#include <string>
#include <vector>
#include <cfloat>
#include <cmath>
#include <bits/stdc++.h>
#include <sys/stat.h>
#include <errno.h>

using std::ios;
using std::sort;
using std::string;
using std::vector;

#define SIMILARITY 0
#define DISTANCE 1

#include "heap.cpp"
#include "linked_list.hpp"
#include "nngraph.hpp"
#include "util.hpp"

#include "gclu_options.h"

#include "timer.h"
#include "util.cpp"
#include "linked_list.cpp"
#include "nngraph.cpp"
#include "graphclu.h"

#include "graphclu_lib.cpp"

using namespace std;

struct stat2 {
  int num_calc_clu_dist;
  int num_pruned;
};

int g_use_heap = 0;
struct stat2 g_stat;

void print_stat() {
  printf("STAT num_calc_clu_dist=%d num_pruned=%d\n", g_stat.num_calc_clu_dist, g_stat.num_pruned);
}

#include <Python.h>
#include <cstring>
#include <cmath>
#include <fstream>
#include <csignal> // Raise
#include <fstream>
#include <vector>
#include <float.h>
#include <math.h>
#include <numpy/arrayobject.h>
// #ifdef Py_PYTHON_H
// #include "rknng_lib.h"
#include <stdio.h>
// #include "rknng/rknng_lib.h"
// #include "dencl/dencl.hpp"

// using namespace std;

#define v(x0, x1)                                                                                  \
  (*(npy_float64 *)((PyArray_DATA(py_v) + (x0)*PyArray_STRIDES(py_v)[0] +                          \
                     (x1)*PyArray_STRIDES(py_v)[1])))
#define v_shape(i) (py_v->dimensions[(i)])

PyObject *array_to_py(int *arr, int N) {
  // Convert c array to python format
  // printf("array_to_py size=%d\n", N);
  PyObject *pyarr = PyList_New(N);
  for (int i = 0; i < N; i++) {
    PyList_SetItem(pyarr, i, Py_BuildValue("i", arr[i]));
  }
  return pyarr;
}

PyObject *py_gclu(PyListObject *py_v, int num_clusters, int num_tsp, int dfunc) {

  // PyObject *ret;
  PyObject *py_labels;

  // PyObject *PyList_GetItem(PyObject *list, Py_ssize_t index)¶
  // Py_ssize_t PyList_Size(PyObject *list)¶
  int lsize = PyList_Size((PyObject *)py_v);
  printf("size=%d\n", lsize);
  int numnodes = 0;
  for (int i = 0; i < lsize; i++) {
    PyObject *x = PyList_GetItem((PyObject *)py_v, i);
    PyObject *idA = PyList_GetItem((PyObject *)x, 0);
    PyObject *idB = PyList_GetItem((PyObject *)x, 1);
    PyObject *weight = PyList_GetItem((PyObject *)x, 2);
    int a = PyLong_AsLong(idA);
    int b = PyLong_AsLong(idB);
    if (numnodes < a) {
      numnodes = a;
    }
    if (numnodes < b) {
      numnodes = b;
    }
  }
  numnodes = numnodes + 1; // Indexes are from 0 to (N-1)
  printf("nodes=%d\n", numnodes);

  nnGraph *graph = init_nnGraph(numnodes);

  for (int i = 0; i < lsize; i++) {
    PyObject *x = PyList_GetItem((PyObject *)py_v, i);
    PyObject *idA = PyList_GetItem((PyObject *)x, 0);
    PyObject *idB = PyList_GetItem((PyObject *)x, 1);
    PyObject *weight = PyList_GetItem((PyObject *)x, 2);
    int a = PyLong_AsLong(idA);
    int b = PyLong_AsLong(idB);
    float w = PyFloat_AsDouble(weight);
    // nng_add_neighbor_safe(graph, a, b, w);
    // nng_add_neighbor_safe(graph, b, a, w);

    if (!nng_has_neighbor(graph, a, b)) {
      nng_add_neighbor(graph, a, b, w);
    }

    if (!nng_has_neighbor(graph, b, a)) {
      nng_add_neighbor(graph, b, a, w);
    }
  }
  
  
  // nnGraph *graph2 = read_ascii_graphf("data/s4_knng_k30.txt");
  // graph = graph2;
  scale_weights(graph, 1); // DISTANCE
  // TODO

  // if (g_opt.scale == 1) {
  // if (graph_type == DISTANCE) {
  // scale_weights(graph, 1);
  // } else {
  // scale_weights(graph, 2);
  // }
  // } else {
  // // No scale
  // scale_weights(graph, -1);
  // }

  Clustering *clu;
  init_Clustering(&clu, graph->size, num_clusters);

  g_opt.repeats = 0;
  g_opt.dissolve = 1;
  g_opt.clusters = 50;
  g_opt.graph_type = SIMILARITY;
  // g_opt.graph_type = DISTANCE;
  g_opt.debug = 0;
  g_opt.grow_factor = 0.8;
  g_opt.grow_factor_start = 1.0;
  g_opt.grow_factor_end = 1.0;
  g_opt.verbose = 1;
  g_opt.costf = 1;    // 1=cond, 2=meanw
  g_opt.minimize = 0; // maximize by default
  g_opt.debug_save_intermediate_part = 0;
  g_opt.density_method = 1;
  g_opt.max_iter = 200;
  g_opt.scale = 1;

  if (g_opt.costf <= 3) {
    g_opt.minimize = 1;
    g_opt.costmultip = -1.0;
  } else {
    g_opt.minimize = 0;
    g_opt.costmultip = 1.0;
  }

  g_timer.tick();
  int seed = 8833;
  rand_seed(seed);

  m_algo(graph, clu, 10 /*repeats*/, num_clusters);

  for (int i = 0; i < clu->N; i++) {
    clu->part[i] += 1;
  }
  py_labels = array_to_py(clu->part, clu->N);

  return py_labels;
}

extern "C" {

static PyObject *gclu_py(PyObject *self, PyObject *args, PyObject *kwargs);

// Define python accessible methods
static PyMethodDef GCluMethods[] = {
    {"gclu", gclu_py, METH_VARARGS | METH_KEYWORDS, "Graph clustering TODO."},
    {NULL, NULL, 0, NULL}};

#define v(x0, x1)                                                                                  \
  (*(npy_float64 *)((PyArray_DATA(py_v) + (x0)*PyArray_STRIDES(py_v)[0] +                          \
                     (x1)*PyArray_STRIDES(py_v)[1])))
#define v_shape(i) (py_v->dimensions[(i)])

/* This initiates the module using the above definitions. */
// #if PY_VERSION_HEX >= 0x03000000
static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT, "gclu", NULL, -1, GCluMethods, NULL, NULL, NULL, NULL};

PyMODINIT_FUNC PyInit_gclu(void) {
  PyObject *m;
  m = PyModule_Create(&moduledef);
  if (!m) {
    return NULL;
  }
  return m;
}

static PyObject *gclu_py(PyObject *self, PyObject *args, PyObject *kwargs) {
  import_array();
  // PyArrayObject *py_v;
  // PyObject *py_v;
  PyListObject *py_v;
  int num_neighbors = 10, num_clusters = 10, num_tsp = 10, maxiter = 100;
  float nndes_start = 0.0, endcond = 0.05;
  char *type = NULL;
  char *distance = NULL;
  // int dfunc = D_L2;

  PyObject *ret;
  static char *kwlist[] = {"v", "num_clusters", "num_tsp", "dtype", "distance", NULL};

  if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O!i|iss", kwlist, &PyList_Type, &py_v,
                                   &num_clusters, &num_tsp, &type, &distance)) {
    return NULL;
  }

  if (num_tsp <= 0) {
    PyErr_SetString(PyExc_ValueError, "num_tsp <= 0 ");
  }

  // printf("v %f %f\n", v(0, 0), v(0, 1));

  // int dtype = D_L2;
  // if (distance != NULL) {
  // if (strcmp("l2", distance) == 0) {
  // dtype = D_L2;
  // } else if (strcmp("l1", distance) == 0) {
  // dtype = D_L1;
  // } else if (strcmp("cos", distance) == 0) {
  // dtype = D_COS;
  // } else {
  // PyErr_SetString(PyExc_ValueError, "Distance must be one for {l2(default),l1,cos}");
  // return NULL;
  // }
  // }

  printf("gclu_py 0001\n");

  ret = py_gclu(py_v, num_clusters, num_tsp, 0);

  return ret;
}

} // END extern "C"
