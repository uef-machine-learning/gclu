
# GCLU
The GCLU software implements the K-algorithm and M-algorithm for graph clustering.

# Input data format.
One line for each item in dataset. Line format:  
 - First number: id of node (in range 0..(N-1))
 - Second number: K = Number of neighbors
 - Next K values: the id:s of the K neighbors
 - Next K values: the weights  to the K neighbors
 
Lines end with '\n' (including last line). Input file must not have \r or \t characters.

Example row: 
0 5 24 99 483 11 444 1414008.500000 1238420.750000 0.000000 0.000000 0.000000  
...  
First node id = 0, K = 5, neighbors = [24 99 483 11 444], weights = [1414008.500000 1238420.750000 .... ]  

So, input is actually a directed graph, but gets converted to undirected. 

See unb3_knn10.txt as example input file. unb3_knn10.txt is a kNN graph of original unb3.txt dataset.

# Compile (on linux)
```make```

If needed to run other algorithms than k-algo, (e.g. walktrap, louvain, fastg), need  to set location of iGraph in Makefile: OIGRAPH var. Then compile:
```make gclu_ig```
# Running the program
```gclus [--help] [-o <file>] <file> [--seed=<n>] [-R <n>] [-K <n>] [-g <num>] [-V <num>] [-H <num>] [--density=<n>] [-I <num>] [--costf=<num>] [--format=<ascii|binary>] [-A <algorithm>] [--type=<similarity|distance>] [--scale=<no|yes>] [--evalpart=<file>] [--savparts=<num>] [--gstart=<num>] [--gend=<num>]
  --help                    display this help and exit
  -o, --out=<file>          output file
  <file>                    input graph file
  --seed=<n>                random number seed
  -R, --repeats=<n>         Number of repeats
  -K, --clusters=<n>        Number of clusters
  -g, --growf=<num>         grow factor [0,1]
  -V, --verbose=<num>       Verbose level (default = 1)
  -H, --write-header=<num>  Write header for result partition (0 =no, 1=yes (default))
  --density=<n>             Density method, one of {0,1,2}
  -I, --maxiter=<num>       Maximum number of iterations for the k-algo
  --costf=<num>             cost function {1=cond,2=inv(default),meanw}
  --format=<ascii|binary>   Input format: ascii or binary
  -A, --algo=<algorithm>    One of {k,m}. (default = m)
  --type=<similarity|distance> Are weights distances or similarities
  --scale=<no|yes>          Scale weights automatically (default:yes)
  --evalpart=<file>         Evaluate partition
  --savparts=<num>          (debug) Save intermediate result partitions. 1 = once for each run of k-algo. 2 = for every k-algo iteration
  --gstart=<num>            (experimental) grow factor range start
  --gend=<num>              (experimental) grow factor range end
```

K-algo run example:
```./gclu data/s4_knng_k30.txt -o tmp/s4.part -K 15 --algo k --costf 1 --type distance```

M-algo run example (better quality):
```./gclu data/s4_knng_k30.txt -o tmp/s4.part -R 100 -K 15 --algo m --costf 1 --type distance```

Abot parameters:
 - "-K 3": Cluster data to 3 clusters
 - For better results, increase -R parameter to e.g. 1000
 - For costf, can try other values. "--costf 2" gives more balanced cluster sizes.
 - "--type similarity" larger weights in the graph mean that nodes are closer
 - "--type distance" smaller weights in the graph mean that nodes are closer
 
Output file contains number of nodes and number of clusters as header info and then the cluster label as integer for all nodes.

Save intermediate partitions:
```./gclu data/s4_knng_k30.txt -o tmp/s4.part -R 100 -K 15 --algo k --costf 1 --type distance --seed 1632488925 -g 0.8 --savparts=1```
Will create a directory 'part' where the intermediate partitions will be stored in. With savparts=1 will store only results for successful merge&split operations.

