#!/usr/bin/python3

import networkx as nx

# import numpy as np
# import matplotlib.pyplot as plt
import random

from gclu import gclu

infname ="data/s4_knng_k30.txt"
G = nx.Graph()
edges = []
# edges= np.ndarray(shape=(0,3),dtype=float)
with open(infname) as f:
    # x = [int(x) for x in next(f).split()] # read first line
    
    for line in f: # read rest of lines
    # x = next(f).split() # read first line
        # G.add_edge(i, i+1, weight=0.000001)
        x = line.split() # read first line
        i = int(x[0])
        num = int(x[1])
        # print(x)
        # print(i)
        # print(num)
        for c in range(0, num):
            # print(c)
            nid = int(x[2+c])
            w = float(x[2+num+c])
            # print("from=%d to=%d w=%f" % (i,nid,w)) 
            # G.add_edge(i, nid, weight=w)
            edges.append([i,nid,w])
            
# labels=gclu(edges,15)
labels=gclu(edges,graph_type="distance",num_clusters=13,repeats=10)
# labels=gclu(edges,graph_type="similarity",num_clusters=13,repeats=100)
print(labels)

G = nx.karate_club_graph()

