#!/bin/python
#
#
# glicko_analysis.py
#
#
# LABELS
# author: Patrik.bey@bih-charite.de
# 
#
#
# CONTENT
# This script contains the analysis of glick ranking history
# performed in 1.
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

import os, numpy, matplotlib.pyplot as plt, scipy.stats, matplotlib.lines as mlines


#####################################
#                                   #
#             FUNCTIONS             #
#                                   #
#####################################


def load_data(batch):
    '''
    load batch file
    '''
    temp = numpy.genfromtxt(os.path.join(Path,f'GlickoHistoryB{str(batch)}.csv'), delimiter=',')
    temp=temp[1:,1:]
    return(temp)

def plot_single_animals(_ax, _data):
    '''
    plot each animal according to its final rank
    with corresponding line type
    '''
    _style = [':','-.','--','-']
    _order = numpy.argsort(_data[:,-1], kind='stable')
    if abs(_data).max() in _data[:,-1]:
        _order = list(_order[::-1])
    else:
        _order = list(_order)
    _data = _data / abs(_data).max()
    for i in numpy.arange(4):
        _ax.plot(_data[_order.index(i),:], ls = _style[i], color = 'black')

def plot_single_animals(_ax, _data):
    '''
    plot each animal according to its final rank
    with corresponding line type
    '''
    _style = [':','-.','--','-']
    _order = sort_ranks(_data[:,-1])
    _data = _data / abs(_data).max()
    for i in numpy.arange(4):
        _ax.plot(_data[_order[i],:], ls = _style[i], color = 'black', linewidth = 3)

def sort_ranks(ranks):
    '''
    sort final glicko ranks
    '''
    order=[]
    order.append(numpy.where(ranks==ranks.min())[0])
    mranks = numpy.ma.array(ranks, mask = False)
    maxid = numpy.where(ranks==ranks.max())
    mranks.mask[order[0]] = True
    mranks.mask[maxid[0]] = True
    order.append(numpy.where(mranks==mranks.min())[0])
    order.append(numpy.where(mranks==mranks.max())[0])
    order.append(maxid[0])
    order=numpy.asarray(order).reshape(-1)
    return(order)

def get_bins(_vector, _bincount = 20):
    '''
    subdivide vector into 5% bins for plotting as boxplots
    '''
    _dim = len(_vector)
    bin_size = round((1/_bincount)*_dim)
    mismatch = numpy.mod(_dim, bin_size)
    if mismatch > 0:
        extension_size = bin_size - mismatch
        _vector = numpy.concatenate([numpy.repeat(numpy.nan,extension_size),_vector])
    elif mismatch == 0:
        _bincount = _dim/bin_size
    return(numpy.split(_vector, _bincount))

def plot_clean():
    '''
    removing borders of plot to fix previous plots
    '''
    plt.gca().spines['right'].set_visible(False)
    plt.gca().spines['top'].set_visible(False)
    plt.gca().yaxis.set_ticks_position('left')
    plt.gca().xaxis.set_ticks_position('bottom')

def plot_corr_boxes(_list):
    '''
    plot boxplot for all 5% bins in _list object returned by
    get_bins()
    '''
    import statsmodels.api as sm
    n_bins = len(_list)
    # plt.axes([0, 0, 1, 1], frameon=True)
    for b in numpy.arange(n_bins):
        temp = _list[b]
        plt.boxplot(temp[~numpy.isnan(temp)], positions=[b])
    plt.gca().spines['right'].set_visible(False)
    plt.gca().spines['top'].set_visible(False)
    plt.gca().yaxis.set_ticks_position('left')
    plt.gca().xaxis.set_ticks_position('bottom')
    smoothing = sm.nonparametric.lowess(numpy.mean(numpy.array(_list), axis = 1),numpy.arange(n_bins), frac = 1./2)
    plt.plot(smoothing[:,0],smoothing[:,1], color = 'darkgrey')


# Plot alpha male number of interactions until ranking does not change anymore
def get_dominance_point(_array):
    '''
    extract the last time first time point where final dominant animal doesn't change anymore
    '''
    id = numpy.where(_array[:,-1]==_array[:,-1].max())[0][0]
    if _array[id,-2] != _array[:,-2].max():
        return(_array.shape[1])
    order = []
    test_index = numpy.sort(numpy.arange(_array.shape[1]-1))[::-1]
    for i in test_index:
        if _array[id,i] == _array[:,i].max():
            order.append(i)
        else:
            break;
    return(numpy.array(order).min())

#####################################
#                                   #
#         DEFINE VARIABLES          #
#                                   #
#####################################

#---- get study folder ----#
Path = f'{os.getcwd()}/Data'

#---- get animal batch order ----#
Batches=(5,2,6,10,7,1,8,9,11,12) # ordered (control / knock out)



#####################################
#                                   #
#          DO COMPUTATIONS          #
#                                   #
#####################################

#---- extract final ranking ----#

LastRating = numpy.zeros([4,len(Batches)])
LastRatingRescale = numpy.zeros([4,len(Batches)])
MiddleRatingRescale = numpy.zeros([4,len(Batches)])
LastRank = numpy.zeros([4,len(Batches)])
MiddleRank = numpy.zeros([4,len(Batches)])

for b in Batches:
    temp = numpy.genfromtxt(os.path.join(Path,f'GlickoHistoryB{str(b)}.csv'), delimiter=',')
    temp=temp[1:,1:]
    rescale = temp / abs(temp).max()
    middle = int(temp.shape[1]/2)
    LastRating[:,Batches.index(b)] = temp[:,-1]
    LastRatingRescale[:,Batches.index(b)] = rescale[:,-1]
    MiddleRatingRescale[:,Batches.index(b)] = rescale[:,middle]
    LastRank[:,Batches.index(b)] = sort_ranks(rescale[:,-1])
    MiddleRank[:,Batches.index(b)] = sort_ranks(rescale[:,middle])



#---- compure rank correlations ----#

RankCorrs = dict()
for b in Batches:
    temp = numpy.genfromtxt(os.path.join(Path,f'GlickoHistoryB{str(b)}.csv'), delimiter=',')
    temp=temp[1:,1:]
    Corr = numpy.zeros(temp.shape[1])
    LastRank = temp[:,-1]
    for i in numpy.arange(len(Corr)):
        if numpy.sum(temp[:,i])==0:
            Corr[i] = 0
        else:
            # Corr[i] = scipy.stats.spearmanr(temp[:,i], LastRank)[0]
            Corr[i] = scipy.stats.pearsonr(temp[:,i], LastRank)[0]
    RankCorrs[b] = Corr


#---- compute despotism score ----#
# Despotism was defined as the ratio of exhibited alpha power (absolute difference in glicko rating of alpha over beta animal) over the total power exhibited (absolute difference in glicko rating of alpha over delta animal).

Alpha = numpy.zeros([len(Batches)])
Delta = numpy.zeros([len(Batches)])
for b in Batches:
    temp = LastRating[:,Batches.index(b)]
    if temp.min() <0:
        Delta[Batches.index(b)] = temp.max() + abs(temp.min())
    else:
        Delta[Batches.index(b)] = temp.max() - temp.min()
    alpha_idx = numpy.where(temp == temp.max())[0]
    subset = temp[numpy.arange(len(temp))!=alpha_idx]
    beta_idx = numpy.where(temp == subset.max())[0]
    if temp[beta_idx] < 0:
        Alpha[Batches.index(b)] = temp[alpha_idx] + abs(temp[beta_idx])
    else:
        Alpha[Batches.index(b)] = temp[alpha_idx] - temp[beta_idx]


Despotism = Alpha / Delta