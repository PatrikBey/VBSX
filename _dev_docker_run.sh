# Docker container run snippets

docker build . -t vbs:dev



docker run -it -v C:\\Users\\me\\GitHub\\VBSX\\Data:/data vbs:dev bash

python3

file='/data/Glicko_results_2018-10-24.Rdata'




# R
R

load(file)

temp = glicko$batch$B1$history
temp[,,1]

Batches=names(glicko$batch)

for (b in Batches) {
    temp=glicko$batch[[b]]$history
    write.csv(temp[,,1], sprintf("/data/GlickoHistory%s.csv", b), row.names=TRUE)
}

write.csv(temp[,,1], "/data/test.csv", row.names=TRUE)



# python
python3
import os, numpy, matplotlib.pyplot as plt, scipy
Path='/data' 

Batches=(1,10,11,12,2,5,6,7,8,9)
 
plt.figure(figsize=(40,20))

for b in Batches:
    temp = numpy.genfromtxt(os.path.join(Path,f'GlickoHistoryB{str(b)}.csv'), delimiter=',')
    temp=temp[1:,1:]
    plt.subplot(2,5,Batches.index(b)+1)
    plt.plot(temp.T)
    plt.title(f'Batch{str(b)}')


plt.savefig(f'/data/GlickoHistory-classic.png')
plt.tight_layout()
plt.close()





# VarPlot

plt.figure(figsize=(40,20))

for b in Batches:
    temp = numpy.genfromtxt(os.path.join(Path,f'GlickoHistoryB{str(b)}.csv'), delimiter=',')
    temp=temp[1:,1:]
    vardynamic = numpy.var(temp, axis=0)
    meandynamic = numpy.mean(temp, axis=0)
    plt.subplot(2,5,Batches.index(b)+1)
    # plt.plot(meandynamic, color = 'black')
    maxid = numpy.where(temp[:,-1]==numpy.max(temp[:,-1]))[0]
    plt.plot(temp[maxid,:].T, color='purple', label='final alpha male')
    plt.errorbar(numpy.arange(len(vardynamic)),meandynamic,vardynamic, alpha = .5, color='black', label='rank dynamic')
    plt.title(f'Batch{str(b)}')

plt.legend()


plt.savefig(f'/data/GlickoHistory-Variance.png')
plt.tight_layout()
plt.close()



plt.plot(meandynamic)
plt.errorbar(numpy.arange(len(vardynamic)),meandynamic,vardynamic)
plt.title(f'Batch{str(b)}-histor-error')
plt.savefig(f'/data/Batch{str(b)}-history-error.png')
plt.close()


scipy.stats.entropy(temp[:,0])
