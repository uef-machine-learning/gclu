#!/usr/bin/python3

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
#          Larger values mean better optimization.
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
