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
#include "util.hpp"

#include "gclu_options.h"

#include "timer.h"
#include "util.cpp"

// #include "linked_list.h"
// #include "nngraph.h"


#include "linked_list.cpp"
#include "nngraph.cpp"

#include "graphclu.h"

#ifdef USE_IGRAPH
#include "igraph_algos.h"
#endif

// double costf_w(double intsum, int size, nnGraph *graph, Clustering *clu);

#include "graphclu_lib.cpp"



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
  
  test_nn_graph();

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
    // g_opt.debug_save_intermediate_part = 1;
    g_opt.debug_save_intermediate_part = a_saveiterparts->ival[0];
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
    printf("\ngclu");
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
