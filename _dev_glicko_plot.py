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
import os, numpy, matplotlib.pyplot as plt, scipy.stats
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
    axs[x,id].set_title(str(Batches.index(b)+1))
    if id == 0:
        axs[x,id].set(ylabel='glicko rating')

plt.tight_layout()
plt.savefig(f'/{Path}/GlickoHistory-classic.png')
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
#     plt.scatter(MiddleRatingRescale[:,i],LastRatingRescale[:,i])

# last rating rescaled
plt.subplot(1,2,1)
for i in numpy.arange(5):
    plt.scatter(numpy.repeat(i+1,4),LastRatingRescale[:,i],label=f'B{Batches[i]}', color='black')

plt.plot(numpy.arange(1,6),numpy.repeat(0,5),ls='--',color = 'darkgray')
plt.title('wildtype')
plt.subplot(1,2,2)
for i in numpy.arange(5,10):
    plt.scatter(numpy.repeat(i-4,4),LastRatingRescale[:,i],label=f'B{Batches[i]}', color='black')

plt.yticks(ticks = numpy.arange(-1,1.25,0.25), labels=[])
plt.plot(numpy.arange(1,6),numpy.repeat(0,5),ls='--',color = 'darkgray')
plt.title('KO')
plt.show()

# last rating
plt.subplot(1,2,1)
for i in numpy.arange(5):
    plt.scatter(numpy.repeat(i+1,4),LastRating[:,i],label=f'B{Batches[i]}', color='black')

# plt.plot(numpy.mean(LastRating[:,0:5], axis = 0), ls='--', label='mean')
plt.plot(numpy.arange(1,6),numpy.repeat(0,5),ls='--',color = 'darkgray')
plt.title('wildtype')
plt.subplot(1,2,2)
for i in numpy.arange(5,10):
    plt.scatter(numpy.repeat(i-4,4),LastRating[:,i],label=f'B{Batches[i]}', color='black')

# plt.plot(numpy.mean(LastRating[:,5:10], axis = 0), ls='--', label='mean')
plt.plot(numpy.arange(1,6),numpy.repeat(0,5),ls='--',color = 'darkgray')
plt.title('KO')
plt.show()

# last rating rescaled within genotype at last timepoint
plt.subplot(1,2,1)
for i in numpy.arange(5):
    temp = LastRating[:,i] / abs(LastRating[:,0:5]).max()
    plt.scatter(numpy.repeat(i+1,4),temp,label=f'B{Batches[i]}', color='black')
    plt.ylim(-1.1,1.1)

plt.plot(numpy.arange(1,6),numpy.mean(LastRating[:,0:5], axis = 0), ls='--', label='mean')
plt.title('wildtype')

plt.subplot(1,2,2)
for i in numpy.arange(5,10):
    temp = LastRating[:,i] / abs(LastRating[:,5:10]).max()
    plt.scatter(numpy.repeat(i-4,4),temp,label=f'B{Batches[i]}', color='black')
    plt.ylim(-1.1,1.1)

plt.plot(numpy.arange(1,6),numpy.mean(LastRating[:,5:10], axis = 0), ls='--', label='mean')
plt.title('KO')

plt.show()


# last rating rescaled at last timepoint
plt.subplot(1,2,1)
for i in numpy.arange(5):
    temp = LastRating[:,i] / abs(LastRating).max()
    plt.scatter(numpy.repeat(i+1,4),temp,label=f'B{Batches[i]}', color='black')
    plt.ylim(-1.1,1.1)

plt.plot(numpy.arange(1,6),numpy.mean(LastRating[:,0:5], axis = 0), ls='--', label='mean')
plt.title('wildtype')

plt.subplot(1,2,2)
for i in numpy.arange(5,10):
    temp = LastRating[:,i] / abs(LastRating).max()
    plt.scatter(numpy.repeat(i-4,4),temp,label=f'B{Batches[i]}', color='black')
    plt.ylim(-1.1,1.1)

plt.plot(numpy.arange(1,6),numpy.mean(LastRating[:,5:10], axis = 0), ls='--', label='mean')
plt.title('KO')

plt.show()


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


scipy.stats.entropy(temp[:,0])
