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

#---- define phenotypes ----#
phenotypes = list(['+/+','-/-'])

#---- get animal batch order ----#
Batches=(5,2,6,10,7,1,8,9,11,12) # ordered (control / knock out)

WTidx = numpy.arange(5)
KOidx = numpy.arange(5,10)

#---- example batches ----#
FinalBatches = list([6,12])

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

# Alpha = numpy.zeros([len(Batches)])
# Delta = numpy.zeros([len(Batches)])
# for b in Batches:
#     temp = LastRating[:,Batches.index(b)]
#     if temp.min() <0:
#         Delta[Batches.index(b)] = temp.max() + abs(temp.min())
#     else:
#         Delta[Batches.index(b)] = temp.max() - temp.min()
#     alpha_idx = numpy.where(temp == temp.max())[0]
#     subset = temp[numpy.arange(len(temp))!=alpha_idx]
#     beta_idx = numpy.where(temp == subset.max())[0]
#     if temp[beta_idx] < 0:
#         Alpha[Batches.index(b)] = temp[alpha_idx] + abs(temp[beta_idx])
#     else:
#         Alpha[Batches.index(b)] = temp[alpha_idx] - temp[beta_idx]

# Despotism = Alpha / Delta

Despotism_old = (0.67,0.24,0.34,0.35,0.52,0.67,0.36,0.08,0.29,0.1)
'''
Using initial R-based despotism ranks in manuscript:
despotism:  (0.67,0.24,0.34,0.35,0.52,0.67,0.36,0.08,0.29,0.1)
'''


#---- compute dominance score ----#

Dominance = numpy.zeros([1,len(Batches)])
Interactions = numpy.zeros([1,len(Batches)])
for b in Batches:
    temp = numpy.genfromtxt(os.path.join(Path,f'GlickoHistoryB{str(b)}.csv'), delimiter=',')
    temp=temp[1:,1:]
    Dominance[0,Batches.index(b)] = int(get_dominance_point(temp))
    Interactions[0,Batches.index(b)] = int(temp.shape[1])


# DominancePoints
# TpH2+/+
DPWT = Dominance[:,WTidx].astype(int)
ICWT = Interactions[:,WTidx].astype(int)
DPWT.mean() #85.8
DPWT.std() # 54.13
numpy.median(DPWT) #61
# RATIOS
RDPWT = DPWT / ICWT
# TpH2-/-
DPKO = Dominance[:,KOidx].astype(int)
ICKO = Interactions[:,KOidx].astype(int)
DPKO.mean() # 269.2
DPKO.std() # 122.53
numpy.median(DPKO) # 243
# RATIOS
RDPKO = DPKO / ICKO

print(f'median ratio dominance point WT: {numpy.median(RDPWT)}') #68.75%
print(f'median ratio dominance point KO: {numpy.median(RDPKO)}') # 80.08%
# RDPWT.mean() # 64.46%
# RDPWT.std() # 34.08%
# RDPKO.mean() # 80.3%
# RDPKO.std() # 12.43%


#####################################
#                                   #
#           VISUALIZATION           #
#                                   #
#####################################


fig, axs = plt.subplots(2,3, sharey = False, figsize=(25,10), gridspec_kw={'width_ratios':[3.5,1,1]})


fontsize = 20

for i in FinalBatches:
    temp = load_data(i)
    plot_single_animals(axs[FinalBatches.index(i),0],temp)
    axs[FinalBatches.index(i),0].set_title(f'A. Glicko history \nTph2 {phenotypes[FinalBatches.index(i)]}' , loc = 'left', fontsize = fontsize+5)
    axs[FinalBatches.index(i),0].spines['right'].set_visible(False)
    axs[FinalBatches.index(i),0].spines['top'].set_visible(False)
    axs[FinalBatches.index(i),0].yaxis.set_ticks_position('left')
    axs[FinalBatches.index(i),0].xaxis.set_ticks_position('bottom')
    axs[FinalBatches.index(i),0].tick_params(width=3, length=7)
    axs[FinalBatches.index(i),0].set_ylabel('Glicko Rating')
    axs[FinalBatches.index(i),0].spines['left'].set_linewidth(4)
    axs[FinalBatches.index(i),0].spines['bottom'].set_linewidth(4)
    axs[FinalBatches.index(i),0].yaxis.label.set_fontsize(fontsize)
    axs[FinalBatches.index(i),0].set_ylim(-1.05,1.05)
    axs[FinalBatches.index(i),0].set_xlim(-1,round(Interactions[:,Batches.index(i)][0],-1).astype(int))

axs[1,0].tick_params(width=3, length=7)
axs[1,0].set_xlabel('Interactions', fontsize=fontsize)

for label in ([axs[0,0].title,axs[1,0].title] + axs[0,0].get_xticklabels() + axs[0,0].get_yticklabels() + axs[1,0].get_xticklabels() + axs[1,0].get_yticklabels()):
    label.set_fontsize(fontsize)


### FINAL GLICKO RATING

for i in WTidx:
    temp = LastRating[:,i] / abs(LastRating[:,WTidx]).max()
    axs[0,1].scatter(numpy.repeat(i+1,4),temp,label=f'B{Batches[i]}', s = 100, color='black')

for i in KOidx:
    temp = LastRating[:,i] / abs(LastRating[:,KOidx]).max()
    axs[1,1].scatter(numpy.repeat(i-4,4),temp,label=f'B{Batches[i]}', s=100, color='black',  marker = 'o', facecolor='none')

for i in [0,1]:
    axs[i,1].set_ylim(-1.1,1.1)
    axs[i,1].set_title(f'B. Final Glicko \nTph2 {phenotypes[i]}', loc = 'left',fontsize = fontsize+5)
    axs[i,1].spines['right'].set_visible(False)
    axs[i,1].spines['top'].set_visible(False)
    axs[i,1].yaxis.set_ticks_position('left')
    axs[i,1].xaxis.set_ticks_position('bottom')
    axs[i,1].set_ylabel('Glicko Rating')
    axs[i,1].yaxis.label.set_fontsize(fontsize)
    axs[i,1].spines['left'].set_linewidth(4)
    axs[i,1].spines['bottom'].set_linewidth(4)
    axs[i,1].set_xticks(numpy.arange(1,6), labels=numpy.arange(1,6)+i*5)
    axs[i,1].axhline(0, ls ='--', color='darkgray')
    axs[i,1].tick_params(width=3, length=7)

axs[1,1].set_xlabel('Groups ', fontsize=fontsize)

for label in ([axs[0,1].title,axs[1,1].title] + axs[0,1].get_xticklabels() + axs[0,1].get_yticklabels() + axs[1,1].get_xticklabels() + axs[1,1].get_yticklabels()):
    label.set_fontsize(fontsize)

### DOMINANCE POINTS

#---- random jitter array for spread of data points
x_jitter1 = numpy.array([0.95534285, 1.04809788, 0.99010567, 1.0181319 , 1.0422938 , 1.00470313, 1.00175763, 0.96126289, 1.00909865, 0.97047873])

for i in WTidx:
    axs[0,2].scatter(x_jitter1[i], Dominance[0][i], marker = 'o', facecolor='black', color = 'black', s = 100, label = f'Tphs {phenotypes[0]}')

for i in KOidx:
    axs[0,2].scatter(x_jitter1[i], Dominance[0][i], marker = 'o', facecolor='none', color = 'black', s = 100, label = f'Tphs {phenotypes[1]}')

axs[0,2].set_title('C. Dominance \n', loc = 'left', fontsize = fontsize+5)
axs[0,2].set_xlim(0.75,1.25)
axs[0,2].set_ylim(0,599)
axs[0,2].spines['right'].set_visible(False)
axs[0,2].spines['top'].set_visible(False)
axs[0,2].yaxis.set_ticks_position('left')
axs[0,2].xaxis.set_ticks_position('bottom')
axs[0,2].set_ylabel('Interaction Count')
axs[0,2].yaxis.label.set_fontsize(fontsize)
axs[0,2].spines['left'].set_linewidth(4)
axs[0,2].spines['bottom'].set_linewidth(4)
axs[0,2].set_xticks([])
axs[0,2].legend(loc='best', frameon = False)
axs[0,2].tick_params(width=3, length=7)

#### DESPOTISM
#---- random jitter array for spread of data points
x_jitter2 = numpy.array([1.04197265, 0.98218399, 0.98294308, 1.02805604, 0.97595451, 1.02746764, 0.97661856, 0.97030259, 0.99442258, 1.02175054])

WT = mlines.Line2D([], [], color='black', marker='o', linestyle='None', markersize=10, label='Tph2 +/+')
KO = mlines.Line2D([], [], color='black',fillstyle = 'none',   marker='o', linestyle='None', markersize=10, label='Tph2 -/-')

for i in numpy.arange(10):
    if i < 5:
        filled = 'black'
    else:
        filled = 'none'
    axs[1,2].scatter(x_jitter2[i], Despotism_old[i], marker = 'o', facecolor=filled, color = 'black', s = 100)

axs[1,2].axhline(0.33, ls ='--', color='darkgray')
axs[1,2].set_title('D. Despotism \n', loc = 'left', fontsize = fontsize+5)
axs[1,2].set_xlim(0.75,1.25)
axs[1,2].spines['right'].set_visible(False)
axs[1,2].spines['top'].set_visible(False)
axs[1,2].yaxis.set_ticks_position('left')
axs[1,2].xaxis.set_ticks_position('bottom')
axs[1,2].set_ylabel('Alpha Power')
axs[1,2].yaxis.label.set_fontsize(fontsize)
axs[1,2].spines['left'].set_linewidth(4)
axs[1,2].spines['bottom'].set_linewidth(4)
axs[1,2].set_xticks([])

axs[1,2].tick_params(width=3, length=7)

for label in ([axs[0,2].title,axs[1,2].title] + axs[0,2].get_xticklabels() + axs[0,2].get_yticklabels() + axs[1,2].get_xticklabels() + axs[1,2].get_yticklabels()):
    label.set_fontsize(fontsize)


plt.tight_layout()

plt.savefig(os.path.join(Path,'Figure5.png'))
plt.close()

