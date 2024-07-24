#!/bin/python
#
#
# SNA_plotting.py
#
#
# LABELS
# author: Patrik.bey@bih-charite.de
# 
#
#
# CONTENT
# This script contains the plotting functions for
# Figure 4 (G) of 1.
#
#
#
# REFERENCES
# 1. Rivalan et al. (2024), https://doi.org/10.1101/2024.02.02.578690



#####################################
#                                   #
#        IMPORT LIBRARIES           #
#                                   #
#####################################

import os, numpy, matplotlib.pyplot as plt, scipy.stats, matplotlib.lines as mlines, pyreadr

#####################################
#                                   #
#             FUNCTIONS             #
#                                   #
#####################################

def load_data(filename):
    '''
    load batch csv file
    '''
    array=numpy.genfromtxt(filename, delimiter=',', dtype=str)
    return(array[1:,:])

def get_subset(_array, day=1,behav = 'STR.G', phase = 'Dark'):
    '''
    get day specific subset of behavioral scores for given behvaiour and phase
    '''
    if isinstance(day,int):
        day = str(day)
    behav_idx = numpy.where(_array[:,2]==behav)
    tmp = _array[behav_idx[0],:]
    day_idx = numpy.where(tmp[:,4] == day)
    tmp = tmp[day_idx[0],:]
    phase_idx = numpy.where(tmp[:,0] == phase)
    return(tmp[phase_idx[0],:])

def get_conf_matrix(animals, behav):
    '''
    compute confusion matrix for given behaviour array
    '''
    conf = numpy.zeros([len(animals), len(animals)])
    for a in animals:
        a_idx = numpy.where(behav[:,1] == a)[0]
        for b in animals:
            conf[animals.index(a),animals.index(b)] = len(numpy.where(behav[a_idx,3] == b)[0])
    conf = conf - numpy.eye(len(animals))*conf
    return(conf)

#####################################
#                                   #
#         DEFINE VARIABLES          #
#                                   #
#####################################

#---- get study folder ----#
Path = f'{os.getcwd()}/Data'

#---- get animal batch order ----#
Batches=(5,2,6,10,7,1,8,9,11,12) # ordered (control / knock out)

#---- define example batch for visualization ----#
bex = [2,8]

#---- define behavor ----#
STRG = 'allogrooming'
STRF = 'Struggle at Feeder'
#####################################
#                                   #
#          DO COMPUTATIONS          #
#                                   #
#####################################

for b in bex:
#---- load data ----#
# Data= load_data(os.path.join(Path,f'Batch-B{b}.csv'))
Data = numpy.genfromtxt(os.path.join(Path,'B2_date_behav.csv'), delimiter = ';',dtype=str)
animals = list(numpy.unique(Data[:,1]))
days = list(numpy.unique(Data[:,-1]))
# behav = Data[:,7:10].copy()


for d in days:
    behav = get_subset(Data, day = d)
    plt.subplot(1,5,int(d))
    c = get_conf_matrix(animals, behav)

import networkx


numpy.array([[0,0,0],[0,1,0],[1,0,0],[1,1,0]])


G = networkx.from_numpy_array(c, numpy.array([[0,0],[0,1],[1,0],[1,1]]))
networkx.draw_networkx(G, with_labels=True,arrows = True)
fig, axs = plt.subplots(1,1, sharey = False, figsize=(25,10))




D = networkx.DiGraph()
c[1,0] = 1
nonzero = numpy.where(c!=0)
for i in numpy.arange(len(nonzero[0])):
    D.add_edge(nonzero[0][i],nonzero[1][i], weight = c[nonzero[0][i],nonzero[1][i]])

import numpy as np
np.random.seed(5)

edges = D.edges()
weights = [D[u][v]['weight'] for u,v in edges]
plt.subplot(212)
networkx.draw_circular(D, connectionstyle="arc3,rad=0.2", width = weights)
plt.show()