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
bex = 2
#####################################
#                                   #
#          DO COMPUTATIONS          #
#                                   #
#####################################

#---- load data ----#
Data = dict()

Data= load_data(os.path.join(Path,f'Batch-B{bex}.csv'))


behav = Data[:,7:10].copy()

STR.G = 'allogrooming'
STR.F = 'Struggle at Feeder'

days = '31,1,2,3,4,5'

def get_subset(_array, day=1,behav = ' STR.G'):
    '''
    get day specific subset of behavioral scores for given behvaiour
    '''
    if isinstance(day,int):
        day = str(day)
    behav_idx = numpy.where(_array[:,8]==behav)
    tmp = _array[behav_idx[0],:]
    day_idx = numpy.where(tmp[:,4] == day)