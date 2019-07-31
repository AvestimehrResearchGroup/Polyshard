import numpy as np
import matplotlib.pyplot as plt
import pickle
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
MS = 1000


def plots(fileName):
    with open(fileName, 'rb') as handle:
        result = pickle.load(handle)
    data = result['data']
    numShardsRange = np.array(result['numShardsRange'])
    numNodesRange = np.array(result['numNodesRange'])
    numEpochsRange = np.array(result['numEpochsRange'])
    schemes = list(data.keys())

    # plot 3D
    for s in schemes:
        plot3D(numShardsRange, numNodesRange, numEpochsRange, data, s,
               fileName)
    # plot throughput v.s. epoches for the largest # of shards
    plotEpoch(numShardsRange, numNodesRange, numEpochsRange, data, schemes,
              fileName)
    # plot throughput v.s. # of nodes for the last epoch
    plotNode(numShardsRange, numNodesRange, numEpochsRange, data, schemes,
             fileName)


def plot3D(numShardsRange, numNodesRange, numEpochsRange, data, scheme,
           fileName):
    for method in ['Max', 'Median', 'Mean']:
        t = data[scheme]['tVer' + method]
        thpt = np.repeat(np.transpose(np.matrix(numShardsRange)),
                         len(numEpochsRange), axis=1) / t / MS
        X, Y = np.meshgrid(numEpochsRange, numNodesRange)
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        surf = ax.plot_surface(X, Y, thpt, cmap=cm.coolwarm,
                               linewidth=0, antialiased=False)

        fig.colorbar(surf, shrink=0.5, aspect=5)
        plt.xlabel('Epoch index')
        plt.ylabel('Number of ndoes')
        ax.set_zlabel('Throughput (# blocks/ms)')
        title = 'Data souce: \n' + fileName + '\n' + 'scheme: ' + scheme + \
                ', struggler method: ' + method
        plt.title(title)
        plt.show()


def plotEpoch(numShardsRange, numNodesRange, numEpochsRange, data, schemes,
              fileName):
    for method in ['Max', 'Median', 'Mean']:
        for s in schemes:
            t = data[s]['tVer' + method]
            thpt = numShardsRange[-1] / t[-1, :] / MS
            plt.plot(numEpochsRange, thpt, label=s)
            plt.xlabel('Epoch index')
            plt.ylabel('Throughput (# blocks/ms)')
            title = 'Data souce: \n' + fileName + '\n' + \
                    'scheme: ' + s + ', # of shards: ' + \
                    str(numShardsRange[-1]) + ', struggler method: ' + method
            plt.title(title)
            plt.grid()
            plt.legend(loc='best')
        plt.show()


def plotNode(numShardsRange, numNodesRange, numEpochsRange, data, schemes,
             fileName):
    for method in ['Max', 'Median', 'Mean']:
        for s in schemes:
            t = data[s]['tVer' + method]
            thpt = numShardsRange / t[:, -1] / MS
            plt.plot(numNodesRange, thpt, label=s)
            plt.xlabel('Number of nodes')
            plt.ylabel('Throughput (# blocks/ms)')
            title = 'Data souce: \n' + fileName + '\n' + \
                    'scheme: ' + s + ', # of shards: ' + \
                    str(numShardsRange[-1]) + ', struggler method: ' + method
            plt.title(title)
            plt.grid()
            plt.legend(loc='best')
        plt.show()


fileName = 'all_schemes_dense_K=[5,50]_M=2000_r=3_epoch=100,1000]_s=0.01_1537330.pickle'

plots((fileName))
