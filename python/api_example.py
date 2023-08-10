#!/usr/bin/python3

import networkx as nx

import numpy as np
import matplotlib.pyplot as plt
import random

from gclu import gclu

#Only needed in case of string distance:
# pip install rapidfuzz
# from rapidfuzz.distance import Levenshtein

# x=np.loadtxt('data/s4.txt')
# a=gclu(x,3)
# print(a)

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
            G.add_edge(i, nid, weight=w)
            edges.append([i,nid,w])
            # if i > maxid: maxid=i
            # if nid > maxid: maxid=nid
    # array = []
    # for line in f: # read rest of lines
        # array.append([int(x) for x in line.split()])
a=gclu(edges,15)
print(a)
       

# import pdb; pdb.set_trace()

# def show_clusters_2d(x,labels,numclu):
	# colormap = plt.cm.gist_ncar
	# colorst = [colormap(i) for i in np.linspace(0, 0.9,numclu)]
	# # print(colorst)
	
	# u_labels = np.unique(labels)
	
	# for i in u_labels:
		# plt.scatter(x[labels == i , 0] , x[labels == i , 1] , label = i, color = colorst[i-1])
	# plt.show()


# Fast version using built in distance functions written in C:
# def example_vec(ds,numclu):
	# # For higher quality:
	# #  - increase number of tsp paths (num_tsp), (in range [2,100])
	
	# labels = tspg(ds,numclu,distance="l2",num_tsp=5,dtype="vec")
	# # print(labels)
	# show_clusters_2d(ds,labels,numclu)

# x=np.loadtxt('data/s1.txt')
# example_vec(x,15)

# Slower version using distance function provided by python:
# Can work with any kind of distance.
# Recommended only when no suitable distance function implemented in C++.
# Takes around two minutes for 100k (2D) dataset

# Class implementing distance measure needs to have following properties:
# - attribute 'size' that reflects the number of data objects
# - function distance(a,b) where parameters a,b are integers between 0..(size-1)
# - distance(a,b) must return a float value
# Name of class is arbitrary
# See examples DistanceMeasureL2 and EditDistance


