# Docker container run snippets

# docker build . -t vbs:dev



# docker run -it -v C:\\Users\\me\\GitHub\\VBSX\\Data:/data vbs:dev bash

# python3

# file='/data/Glicko_results_2018-10-24.Rdata'




# # R
# R

# load(file)

# temp = glicko$batch$B1$history
# temp[,,1]

# Batches=names(glicko$batch)

# for (b in Batches) {
#     temp=glicko$batch[[b]]$history
#     write.csv(temp[,,1], sprintf("/data/GlickoHistory%s.csv", b), row.names=TRUE)
# }

# write.csv(temp[,,1], "/data/test.csv", row.names=TRUE)



# python
# python3
import os, numpy, matplotlib.pyplot as plt, scipy.stats, matplotlib.lines as mlines
# Path='/data' 
Path = f'{os.getcwd()}/Data'


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
        _ax.plot(_data[_order[i],:], ls = _style[i], color = 'black')

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


# Batches=(1,10,11,12,2,5,6,7,8,9)
Batches=(5,2,6,10,7,1,8,9,11,12) # ordered (control / knock out)
 


fig, axs = plt.subplots(2, 5, sharey = True, figsize=(30,15))
for b in Batches:
    id = Batches.index(b)
    if id < 5:
        x=0
    else:
        x=1
        id=id-5
    print(f'{b} - {x}')
    temp = load_data(b)
    plot_single_animals(axs[x,id],temp)
    axs[x,id].set_title(str(Batches.index(b)+1), loc = 'left')
    if id == 0:
        if Batches.index(b)<5:
           axs[x,id].set(ylabel='Tph2 +/+ \n glicko rating')
        else:
            axs[x,id].set(ylabel='Tph2 -/- \n glicko rating')

plt.tight_layout()
plt.savefig(f'/{Path}/Plots/GlickoHistory-classic.png')
plt.close()


    
# Last Rating per animal per group
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


# # plot ratings
# plt.subplot(1,2,1)
# for i in numpy.arange(5):
#     plt.scatter(MiddleRatingRescale[:,i],LastRatingRescale[:,i])

# plt.subplot(1,2,2)
# for i in numpy.arange(5,10):
# #     plt.scatter(MiddleRatingRescale[:,i],LastRatingRescale[:,i])

# # last rating rescaled
# plt.subplot(1,2,1)
# for i in numpy.arange(5):
#     plt.scatter(numpy.repeat(i+1,4),LastRatingRescale[:,i],label=f'B{Batches[i]}', color='black')

# plt.plot(numpy.arange(1,6),numpy.repeat(0,5),ls='--',color = 'darkgray')
# plt.title('wildtype')
# plt.subplot(1,2,2)
# for i in numpy.arange(5,10):
#     plt.scatter(numpy.repeat(i-4,4),LastRatingRescale[:,i],label=f'B{Batches[i]}', color='black')

# plt.yticks(ticks = numpy.arange(-1,1.25,0.25), labels=[])
# plt.plot(numpy.arange(1,6),numpy.repeat(0,5),ls='--',color = 'darkgray')
# plt.title('KO')
# plt.show()

# # last rating
# plt.subplot(1,2,1)
# for i in numpy.arange(5):
#     plt.scatter(numpy.repeat(i+1,4),LastRating[:,i],label=f'B{Batches[i]}', color='black')

# # plt.plot(numpy.mean(LastRating[:,0:5], axis = 0), ls='--', label='mean')
# plt.plot(numpy.arange(1,6),numpy.repeat(0,5),ls='--',color = 'darkgray')
# plt.title('wildtype')
# plt.subplot(1,2,2)
# for i in numpy.arange(5,10):
#     plt.scatter(numpy.repeat(i-4,4),LastRating[:,i],label=f'B{Batches[i]}', color='black')

# # plt.plot(numpy.mean(LastRating[:,5:10], axis = 0), ls='--', label='mean')
# plt.plot(numpy.arange(1,6),numpy.repeat(0,5),ls='--',color = 'darkgray')
# plt.title('KO')
# plt.show()

# #########################################################
#
# last rating rescaled within genotype at last timepoint
#
#########################################################

fig, axs = plt.subplots(2, 1, figsize=(5,10))

plt.subplot(2,1,1)
for i in numpy.arange(5):
    temp = LastRating[:,i] / abs(LastRating[:,0:5]).max()
    plt.scatter(numpy.repeat(i+1,4),temp,label=f'B{Batches[i]}', color='black')
    plt.ylim(-1.1,1.1)

plt.ylabel('Final Glicko rating (rescaled)')
plt.xticks(ticks = numpy.arange(1,6), labels=numpy.arange(1,6))
plt.plot(numpy.arange(1,6),numpy.repeat(0,5),ls='--',color = 'darkgray')
# plt.plot(numpy.arange(1,6),numpy.mean(LastRating[:,0:5], axis = 0), ls='--', label='mean')
plt.title('Tph2 +/+')

plt.subplot(2,1,2)
for i in numpy.arange(5,10):
    temp = LastRating[:,i] / abs(LastRating[:,5:10]).max()
    plt.scatter(numpy.repeat(i-4,4),temp,label=f'B{Batches[i]}', color='black')
    plt.ylim(-1.1,1.1)

plt.ylabel('Final Glicko rating (rescaled)')
# plt.yticks(ticks = numpy.arange(-1,1.25,0.25), labels=[])
plt.plot(numpy.arange(1,6),numpy.repeat(0,5),ls='--',color = 'darkgray')
plt.xticks(ticks = numpy.arange(1,6), labels=numpy.arange(6,11))
# plt.plot(numpy.arange(1,6),numpy.mean(LastRating[:,5:10], axis = 0), ls='--', label='mean')
plt.title('Tph2 -/-')
plt.xlabel('Batches')

# plt.tick_params(which = 'both', top=False, bottom=False, left=False, right=False)
# plt.xlabel("Batches")

plt.tight_layout()
plt.show()


# last rating rescaled at last timepoint
# plt.subplot(1,2,1)
# for i in numpy.arange(5):
#     temp = LastRating[:,i] / abs(LastRating).max()
#     plt.scatter(numpy.repeat(i+1,4),temp,label=f'B{Batches[i]}', color='black')
#     plt.ylim(-1.1,1.1)

# plt.plot(numpy.arange(1,6),numpy.mean(LastRating[:,0:5], axis = 0), ls='--', label='mean')
# plt.title('wildtype')

# plt.subplot(1,2,2)
# for i in numpy.arange(5,10):
#     temp = LastRating[:,i] / abs(LastRating).max()
#     plt.scatter(numpy.repeat(i-4,4),temp,label=f'B{Batches[i]}', color='black')
#     plt.ylim(-1.1,1.1)

# plt.plot(numpy.arange(1,6),numpy.mean(LastRating[:,5:10], axis = 0), ls='--', label='mean')
# plt.title('KO')

# plt.show()


# # plot ranks
# plt.subplot(1,2,1)
# for i in numpy.arange(5):
#     plt.scatter(MiddleRank[:,i],LastRank[:,i])

# plt.subplot(1,2,2)
# for i in numpy.arange(5,10):
#     plt.scatter(MiddleRank[:,i],LastRank[:,i])




#########################################################
#                                                       #
#                  RANK CORRELATIONS                    #
#                                                       #
#########################################################

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


# for b in Batches:
#     plt.plot(RankCorrs[b], label=f'B{b}')

# plt.legend()


BatchesList = list(Batches)
for b in BatchesList:
    plt.subplot(2,5,BatchesList.index(b)+1)
    plot_corr_boxes(get_bins(RankCorrs[b]))
    plt.ylim([-1,1.1])
    plt.title(str(BatchesList.index(b)+1), loc = 'left')


plt.show()
# VarPlot

# plt.figure(figsize=(40,20))

# for b in Batches:
#     temp = numpy.genfromtxt(os.path.join(Path,f'GlickoHistoryB{str(b)}.csv'), delimiter=',')
#     temp=temp[1:,1:]
#     vardynamic = numpy.var(temp, axis=0)
#     meandynamic = numpy.mean(temp, axis=0)
#     plt.subplot(2,5,Batches.index(b)+1)
#     # plt.plot(meandynamic, color = 'black')
#     maxid = numpy.where(temp[:,-1]==numpy.max(temp[:,-1]))[0]
#     plt.plot(temp[maxid,:].T, color='purple', label='final alpha male')
#     plt.errorbar(numpy.arange(len(vardynamic)),meandynamic,vardynamic, alpha = .5, color='black', label='rank dynamic')
#     plt.title(f'Batch{str(b)}')

# plt.legend()


# plt.savefig(f'/data/GlickoHistory-Variance.png')
# plt.tight_layout()
# plt.close()



# plt.plot(meandynamic)
# plt.errorbar(numpy.arange(len(vardynamic)),meandynamic,vardynamic)
# plt.title(f'Batch{str(b)}-histor-error')
# plt.savefig(f'/data/Batch{str(b)}-history-error.png')
# plt.close()

#########################################################
#                                                       #
#                        DESPOTISM                      #
#                                                       #
#########################################################


# Plot alpha power for each batch
# alpha power defined as the ratio of power (difference in glicko rating) expressed from alpha male over beta over the full power spectrum from alpha to delta.

Alpha = numpy.zeros([len(Batches)])
Delta = numpy.zeros([len(Batches)])
for b in Batches:
    temp = numpy.genfromtxt(os.path.join(Path,f'GlickoHistoryB{str(b)}.csv'), delimiter=',')
    temp=temp[1:,-1]
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

fig, axs = plt.subplots(1, figsize=(4,7))
x_jitter = numpy.repeat(1.0,10) + numpy.random.uniform(low=-.05, high=0.05, size=(10,))
WT = mlines.Line2D([], [], color='black', marker='o', linestyle='None', markersize=10, label='Tph2 +/+')
KO = mlines.Line2D([], [], color='black',fillstyle = 'none',   marker='o', linestyle='None', markersize=10, label='Tph2 -/-')

for i in numpy.arange(10):
    if i < 5:
        filled = 'black'
    else:
        filled = 'none'
    plt.scatter(x_jitter[i], Despotism_old[i], marker = 'o', facecolor=filled, color = 'black', s = 50)

plt.xlim(0.75,1.25)
plt.ylabel('Alpha power')
plt.xticks([])
plt.legend(handles = [WT,KO])
# plt.legend(['*','+'],['Tph2 +/+','Tph2 -/-'])
plt.show()

numpy.mean(Despotism[:5])
numpy.mean(Despotism[5:])
numpy.var(Despotism[:5])
numpy.var(Despotism[5:])
scipy.stats.ttest_ind(Despotism[:5],Despotism[5:])

Despotism_old = (0.67,0.24,0.34,0.35,0.52,0.67,0.36,0.08,0.29,0.1)
#########################################################
#                                                       #
#                   STABLE DOMINANCE                    #
#                                                       #
#########################################################

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


Dominance = numpy.zeros([1,len(Batches)])
Interactions = numpy.zeros([1,len(Batches)])
for b in Batches:
    temp = numpy.genfromtxt(os.path.join(Path,f'GlickoHistoryB{str(b)}.csv'), delimiter=',')
    temp=temp[1:,1:]
    Dominance[0,Batches.index(b)] = int(get_dominance_point(temp))
    Interactions[0,Batches.index(b)] = int(temp.shape[1])
    # print(f'{get_dominance_point(temp)} for Batch B{b}')



fig, axs = plt.subplots(1, figsize=(4,7))
x_jitter = numpy.repeat(1.0,10) + numpy.random.uniform(low=-.05, high=0.05, size=(10,))
WT = mlines.Line2D([], [], color='black', marker='o', linestyle='None', markersize=10, label='Tph2 +/+')
KO = mlines.Line2D([], [], color='black',fillstyle = 'none',   marker='o', linestyle='None', markersize=10, label='Tph2 -/-')

for i in numpy.arange(10):
    if i < 5:
        filled = 'black'
    else:
        filled = 'none'
    plt.scatter(x_jitter[i], Dominance[0][i], marker = 'o', facecolor=filled, color = 'black', s = 50)

plt.xlim(0.75,1.25)
plt.ylabel('Count')
plt.xticks([])
plt.legend(handles = [WT,KO])
# plt.legend(['*','+'],['Tph2 +/+','Tph2 -/-'])
plt.show()





