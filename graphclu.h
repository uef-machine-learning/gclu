
#ifndef GRAPHCLU_H
#define GRAPHCLU_H

struct nodeItem {
  float density;
  int id;
};


typedef struct {
  int *part; // Partition
  int count;
  int N;                  // Number of objects
  int K;                  // Number of clusters
  int *clusize;           // Sizes of clusters
  int *activity;
  double **node_to_clu_w; // Sum of weights from each node to each cluster
  double *ntr_sums;
  double *ext_sums;
  double *total_sums;
  double *costs;
  double cost;
  double balance_factor;
  double conductance;

  int min_part_size;
  int max_part_size;
  double min_max_part_ratio;
} Clustering;


double costf_w(double intsum, int size, nnGraph *graph, Clustering *clu);
void handler(int sig);
void write_ints_to_file(const char *fn, int *data, int N);
void write_ints_to_fp(FILE *fp, int *data, int N);
void test_nn_graph(); 
void init_Clustering(Clustering **_clu, int N, int K);
void find_max_val(double *data, int N, int *ret_ind, double *ret_val);
void split_and_merge(Clustering *clu, nnGraph *graph);
void grow_Cluster(Clustering *clu, nnGraph *graph, int seed_id, int newpart, int steps);
void free_Clustering(Clustering *clu);
void copy_Clustering(Clustering *clu, Clustering *newclu);
void clone_Clustering(Clustering *clu, Clustering **newclu_);
void track_min_int(int *setv, int val);
void track_max_int(int *setv, int val);
double costf(nnGraph *graph, Clustering *clu, double *r_conductance);
void scale_weights(nnGraph *graph, int invert);
int choose_best_clu(
    // INPUT:
    int nid, nnGraph *graph, Clustering *clu,
    // OUTPUT:
    double *d_cost, double *d_sum_oldpart, double *d_sum_newpart);
double costf_w(double intsum, int size, nnGraph *graph, Clustering *clu);
int choose_best_clu_conductance(
    // INPUT:
    int nid, nnGraph *graph, Clustering *clu,
    // OUTPUT:
    double *d_cost, double *d_sum_oldpart, double *d_sum_newpart, double *d_total_sum);
bool compare_nodeItem(nodeItem i1, nodeItem i2);
void density_init_partition(nnGraph *graph, Clustering *clu);
int k_algo(nnGraph *graph, Clustering *clu);
void swap_algo(nnGraph *graph, Clustering *clu, int n_repeats, int n_clusters);
void repeated_k_algo(nnGraph *graph, Clustering *clu, int n_repeats, int n_clusters);
int main(int argc, char *argv[]);

#endif
