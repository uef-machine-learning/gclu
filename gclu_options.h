#ifndef GCLU_OPTIONS_H
#define GCLU_OPTIONS_H

#include <stdlib.h>

struct gclu_options {
  double grow_factor;
  double grow_factor_start;
  double grow_factor_end;
  int repeats;
  int dissolve;
  int verbose;
  int clusters;
  int graph_type;
  int debug;
  int costf;
  int debug_save_intermediate_part;
  int seed;
  int scale;
  int density_method;
  int max_iter;
  
  int minimize; // 1 = minimize cost function, 0 = maximize
  double costmultip;

};




gclu_options g_opt;

#endif
