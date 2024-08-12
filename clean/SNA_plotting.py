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

import os, numpy, matplotlib.pyplot as plt, networkx, matplotlib.lines as mlines, pyreadr

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

def get_subset(_array, day=1,action = ' STR.G', phase = ' Dark phase'):
    '''
    get day specific subset of behavioral scores for given behvaiour and phase
    '''
    if isinstance(day,int):
        day = str(day)
    behav_idx = numpy.where(_array[:,8]==action)
    tmp = _array[behav_idx[0],:]
    day_idx = numpy.where(tmp[:,-1] == day)
    tmp = tmp[day_idx[0],:]
    phase_idx = numpy.where(tmp[:,3] == phase)
    return(tmp[phase_idx[0],:])

def get_conf_matrix(animals, behav):
    '''
    compute confusion matrix for given behaviour array
    '''
    conf = numpy.zeros([len(animals), len(animals)])
    for a in animals:
        a_idx = numpy.where(behav[:,7] == a)[0]
        for b in animals:
            conf[animals.index(a),animals.index(b)] = len(numpy.where(behav[a_idx,9] == b)[0])
    conf = conf - numpy.eye(len(animals))*conf
    return(conf)

def date2day(array):
    '''
    transform calendar dates to experimental days
    '''
    start_time = array[0,5]
    dates = numpy.unique(array[:,4]).astype(int)
    dates.sort()
    days = numpy.arange(1,6)
    new_date = numpy.zeros([array.shape[0],1])
    idx = numpy.zeros(len(days)-1)
    day = 1
    for d in numpy.arange(len(days)):
        if day == len(days):
            new_date[int(idx[-1]):] = day
        else:
            dateidx = numpy.where(array[:,4] == str(dates[d+1]))[0]
            tmp = array[dateidx,5].astype(int)
            if tmp.min() >= int(start_time):
                timeidx = numpy.where(tmp >= int(start_time))[0]
                idx[day-1] = dateidx[timeidx[0]]
            else:
                timeidx = numpy.where(array[dateidx,5].astype(int) < int(start_time))[0]
                idx[day-1] = dateidx[timeidx[-1]]
            if day==1:
                new_date[:int(idx[day-1])] = day
            else:
                new_date[numpy.arange(idx[day-2],idx[day-1]).astype(int)] = day
        day = day+ 1
    return(numpy.concatenate([array, new_date], axis = 1))

def get_digraph(conf):
    '''
    create directed graph object for given interaction matrix
    '''
    import networkx
    digraph = networkx.DiGraph()
    positions = numpy.array([[0,0],[0,1],[1,0],[1,1]])
    for i in numpy.arange(conf.shape[1]):
        digraph.add_node(i, pos = positions[i,:])
    nonzero = numpy.where(conf!=0)
    for i in numpy.arange(len(nonzero[0])):
        digraph.add_edge(nonzero[0][i],nonzero[1][i], weight = conf[nonzero[0][i],nonzero[1][i]])
    return(digraph)


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
# Example batches: 2,8

#####################################
#                                   #
#          DO COMPUTATIONS          #
#                                   #
#####################################


#---- visualize all batches indvidually ----#

# for b in Batches:
# #---- load data ----#
#     Data= load_data(os.path.join(Path,f'Batch-B{b}.csv'))
#     DataDays = date2day(Data)
#     animals = list(numpy.unique(DataDays[:,7]))
#     days = list(numpy.unique(DataDays[:,-1]))
#     Graphs = dict()
#     plt.figure(figsize=[15,3])
#     for d in days:
#         behav = get_subset(DataDays, day = d)
#         c = get_conf_matrix(animals, behav)
#         Graphs[d] = get_digraph(c)
#     for d in days:
#         plt.subplot(150+days.index(d)+1)
#         plt.ylabel('asd')
#         if d == days[0]:
#             plt.xlabel(f'B{b}')
#         edges = Graphs[d].edges
#         weights = [Graphs[d][u][v]['weight'] for u,v in edges]
#         pos = networkx.get_node_attributes(Graphs[d],'pos')
#         if Batches.index(b) <= 4:
#             networkx.draw_circular(Graphs[d], connectionstyle="arc3,rad=0.2", width = weights*1000, with_labels = False, edgecolors = 'black', node_color = 'Black' , node_size = numpy.array(Graphs[d].degree)[:,1]*200+100)
#         else:
#             networkx.draw_circular(Graphs[d], connectionstyle="arc3,rad=0.2", width = weights*1000, with_labels = False, edgecolors = 'black', node_color = 'none' , node_size = numpy.array(Graphs[d].degree)[:,1]*200+100)
#     plt.tight_layout()
#     plt.show()


#####################################
#                                   #
#        CREATE FIGURE 4(G)         #
#                                   #
#####################################

#---- STRUGGLE AT FEEDER ----#

fig, all_axes = plt.subplots(2, 5, figsize=(25, 10))
#---- EXAMPLE BATCHES----#
Data= load_data(os.path.join(Path,f'Batch-B2.csv'))
DataDays = date2day(Data)
animals = list(numpy.unique(DataDays[:,7]))
days = list(numpy.unique(DataDays[:,-1]))
Graphs = dict()
for d in days:
    behav = get_subset(DataDays, day = d, action = ' STR.F')
    c = get_conf_matrix(animals, behav)
    Graphs[d] = get_digraph(c)

for d in days:
    edges = Graphs[d].edges
    weights = [Graphs[d][u][v]['weight'] for u,v in edges]
    pos = networkx.get_node_attributes(Graphs[d],'pos')
    networkx.draw_circular(Graphs[d], connectionstyle="arc3,rad=0.2", width = numpy.array(weights)/3, arrowstyle = '-', edge_color = 'black', node_color = '#707070', node_size = 1000, ax = all_axes[0,int(float(d))-1])

Data= load_data(os.path.join(Path,f'Batch-B8.csv'))
DataDays = date2day(Data)
animals = list(numpy.unique(DataDays[:,7]))
days = list(numpy.unique(DataDays[:,-1]))
Graphs = dict()
for d in days:
    behav = get_subset(DataDays, day = d, action = ' STR.F')
    c = get_conf_matrix(animals, behav)
    Graphs[d] = get_digraph(c)

for d in days:
    edges = Graphs[d].edges
    weights = [Graphs[d][u][v]['weight'] for u,v in edges]
    pos = networkx.get_node_attributes(Graphs[d],'pos')
    networkx.draw_circular(Graphs[d], connectionstyle="arc3,rad=0.2", arrowstyle = '-', width = numpy.array(weights)/3, edgecolors = 'black', node_color = 'black', linewidths = 7, node_size = 1000, ax = all_axes[1,int(float(d))-1])

all_axes[0,0].set_ylabel('asd',size = 10)
all_axes[1,0].set_ylabel('asd',size = 10)
all_axes[0,0].tick_params(labelleft=True)
plt.tight_layout()
# plt.suptitle('Struggle at feeder', fontsize = 20, x=0.05, y=1)#horizontalalignment='left')
plt.show()


#---- ALLOGROOMING ----#

fig, all_axes = plt.subplots(2, 5, figsize=(29, 12))
#---- EXAMPLE BATCHES----#
Data= load_data(os.path.join(Path,f'Batch-B2.csv'))
DataDays = date2day(Data)
animals = list(numpy.unique(DataDays[:,7]))
days = list(numpy.unique(DataDays[:,-1]))
Graphs = dict()
for d in days:
    behav = get_subset(DataDays, day = d, action = ' STR.G')
    c = get_conf_matrix(animals, behav)
    Graphs[d] = get_digraph(c)

for d in days:
    edges = Graphs[d].edges
    weights = [Graphs[d][u][v]['weight'] for u,v in edges]
    pos = networkx.get_node_attributes(Graphs[d],'pos')
    networkx.draw_circular(Graphs[d], connectionstyle="arc3,rad=0.2", width = numpy.array(weights)/3, arrowstyle = '-', edge_color = 'black', node_color = '#707070', node_size = 1000, ax = all_axes[0,int(float(d))-1])

Data= load_data(os.path.join(Path,f'Batch-B8.csv'))
DataDays = date2day(Data)
animals = list(numpy.unique(DataDays[:,7]))
days = list(numpy.unique(DataDays[:,-1]))
Graphs = dict()
for d in days:
    behav = get_subset(DataDays, day = d, action = ' STR.G')
    c = get_conf_matrix(animals, behav)
    Graphs[d] = get_digraph(c)

for d in days:
    edges = Graphs[d].edges
    weights = [Graphs[d][u][v]['weight'] for u,v in edges]
    pos = networkx.get_node_attributes(Graphs[d],'pos')
    networkx.draw_circular(Graphs[d], connectionstyle="arc3,rad=0.2", arrowstyle = '-', width = numpy.array(weights)/3, edgecolors = 'black', node_color = 'black', linewidths = 7, node_size = 1000, ax = all_axes[1,int(float(d))-1])

all_axes[0,0].set_ylabel('asd',size = 10)
all_axes[1,0].set_ylabel('asd',size = 10)
all_axes[0,0].tick_params(labelleft=True)
plt.tight_layout()
# plt.suptitle('Allogrooming', fontsize = 20, x=0.05, y=1)#horizontalalignment='left')
plt.show()



