
# GCLU
The GCLU software implements the K-algorithm and M-algorithm for graph clustering, published (accepted for publication) in the following article:

 S. Sieranoja and P. Fränti, "Adapting k-means for graph clustering", Knowledge and Information Systems 2021. DOI: 10.1007/s10115-021-01623-y 
 
Datasets can be downloaded from:
http://cs.uef.fi/ml/article/graphclu/

Contact: samisi@cs.uef.fi

# Python interface

## Install & test
```
git clone https://github.com/uef-machine-learning/gclu.git
cd gclu
pip install -r python/requirements.txt
pip install .
python/api_example_knn.py
```

## Examples

See examples how to use in:
```
python/api_example_karate.py
python/api_example_knn.py
```

```py
import networkx as nx
import matplotlib.pyplot as plt

from gclu import gclu

G = nx.karate_club_graph()
    
edges=[]
for nodeA,nodeB in G.edges():
    weight=1.0
    edges.append([nodeA,nodeB,weight])


labels=gclu(edges,graph_type="distance",num_clusters=2,repeats=2,scale="no",seed=121288,costf="inv")
#Explanation of parameters:
# edges: list of [nodeidA,nodeidB,weight] values
#        [[0, 1, 1.0], [0, 2, 1.0], ...]    
# num_clusters: How many groups should the data be divided to
# repeats: How many times to repeat split&merge process. 
           Larger values mean better optimization.
# seed: random seed value
# graph_type: "distance" if larger values mean nodes are farther, 
#             "similarity" if closer.
# costf: is one of {"cond","inv","meanw"}
# scale: set to "no" if weights are already in suitable range. 
#        Algorithm is sensitive to large outliers values.

dlabels={}
for i in range(len(labels)): #transform to dict
    dlabels[i] = labels[i]
    
nx.draw_kamada_kawai(G, with_labels=True,labels=dlabels)

plt.show()
```

# Command line interface
## Input data format.
One line for each item in dataset. Line format:  
 - First number: id of node (in range 0..(N-1))
 - Second number: K = Number of neighbors
 - Next K values: the id:s of the K neighbors
 - Next K values: the weights  to the K neighbors
 
Lines end with '\n' (including last line). Input file must not have \r or \t characters.

Example row: 
```
0 5 24 99 483 11 444 1414008.500000 1238420.750000 0.000000 0.000000 0.000000  
...  
```
First node id = 0, K = 5, neighbors = [24 99 483 11 444], weights = [1414008.500000 1238420.750000 .... ]  

So, input is actually a directed graph, but gets converted to undirected. 

See unb3_knn10.txt as example input file. unb3_knn10.txt is a kNN graph of original unb3.txt dataset.

## Compile (on linux)
```make```

If needed to run other algorithms than k-algo, (e.g. walktrap, louvain, fastg), need  to set location of iGraph in Makefile: OIGRAPH var. Then compile:
```make gclu_ig```
## Running the program
```
./gclu [--help] [-o <file>] <file> [--seed=<n>] [-R <n>] [-K <n>] [-g <num>] [-V <num>] [-H <num>] [--density=<n>] [-I <num>] [--costf=<num>] [--format=<ascii|binary>] [-A <algorithm>] [--type=<similarity|distance>] [--scale=<no|yes>] [--evalpart=<file>] [--savparts=<num>] [--gstart=<num>] [--gend=<num>]
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
  -A, --algo=<algorithm>    One of {k,m}. (default = m)
  --type=<similarity|distance> Are weights distances or similarities
  --scale=<no|yes>          Scale weights automatically (default:yes)
  --evalpart=<file>         Evaluate partition
  --savparts=<num>          (debug) Save intermediate result partitions. 1 = once for each run of k-algo. 2 = for every k-algo iteration
  --gstart=<num>            (experimental) grow factor range start
  --gend=<num>              (experimental) grow factor range end
```

K-algo run example:
```
./gclu data/s4_knng_k30.txt -o tmp/s4.part -K 15 --algo k --costf 1 --type distance
```

M-algo run example (better quality):
```
./gclu data/s4_knng_k30.txt -o tmp/s4.part -R 100 -K 15 --algo m --costf 1 --type distance
```

About parameters:
 - "-K 3": Cluster data to 3 clusters
 - For better results, increase -R parameter to e.g. 1000
 - For costf, can try other values. "--costf 2" gives more balanced cluster sizes.
 - "--type similarity" larger weights in the graph mean that nodes are closer
 - "--type distance" smaller weights in the graph mean that nodes are closer
 
Output file contains number of nodes and number of clusters as header info and then the cluster label as integer for all nodes.

Save intermediate partitions:
```
./gclu data/s4_knng_k30.txt -o tmp/s4.part -R 10 -K 15 --algo m --costf 1 --type distance --seed 1632488925 -g 0.8 --savparts=1
```
Will create a directory 'part' where the intermediate partitions will be stored in. With savparts=1 will store only results for successful merge&split operations.

