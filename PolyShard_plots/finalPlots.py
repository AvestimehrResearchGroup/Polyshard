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
    names = ['full replication', 'uncoded sharding', 'PolyShard']

    for i in range(len(schemes)):
        s = schemes[i]
        plot3D(numShardsRange, numNodesRange, numEpochsRange, data, s,
               names[i], fileName)

    # plot throughput v.s. epoches for the largest # of shards
    plotEpoch(numShardsRange, numNodesRange, numEpochsRange, data, schemes,
              fileName)

    # plot throughput v.s. # of nodes for the last epoch
    plotNode(numShardsRange, numNodesRange, numEpochsRange, data, schemes,
             fileName)


def plot3D(numShardsRange, numNodesRange, numEpochsRange, data, scheme, name,
           fileName):
    for method in ['Mean']:  # ['Max', 'Median']
        t = data[scheme]['tVer' + method]
        thpt = np.repeat(np.transpose(np.matrix(numShardsRange)),
                         len(numEpochsRange), axis=1) / t / MS
        X, Y = np.meshgrid(numEpochsRange, numNodesRange)
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        surf = ax.plot_surface(X, Y, thpt, cmap=cm.coolwarm,
                               linewidth=0, antialiased=False)
        ax.plot_wireframe(X, Y, thpt, color='black', linewidth=0.5)
        fig.colorbar(surf, shrink=0.5, aspect=5)
        plt.xlabel('Epoch index ($t$)', fontsize=12)
        plt.ylabel('Number of nodes ($N$)', fontsize=12)
        ax.set_zlabel('Throughput (blocks/ms)', fontsize=12)
        ax.text(np.mean(X), np.mean(Y), np.max(thpt), name, fontsize=14)
        plt.xlim((np.min(X), np.max(X)))
        plt.ylim((np.min(Y), np.max(Y)))
        plt.tight_layout()
        plt.show()


def plotEpoch(numShardsRange, numNodesRange, numEpochsRange, data, schemes,
              fileName):
    names = ['full replication', 'uncoded sharding', 'PolyShard']
    markers = ['o', 's', 'd']
    # arrowIdx = 1  # draw an arrow at the arrowIdx-th x value.
    # thpt4arrow = []
    for method in ['Mean']:  # ['Max', 'Median']
        for i in range(len(schemes)):
            s = schemes[i]
            t = data[s]['tVer' + method]
            thpt = numShardsRange[-1] / t[-1, :] / MS
            # thpt4arrow.append(thpt[arrowIdx])
            plt.plot(numEpochsRange, thpt, label=names[i], linewidth=2,
                     marker=markers[i], markerfacecolor='w', markersize=10)
            plt.xlabel('Epoch index ($t$)', fontsize=12)
            plt.ylabel('Throughput (blocks/ms)', fontsize=12)
            title = 'Throughput v.s. Time when $K$=' + \
                    str(np.max(numShardsRange)) + ', $N$=' + \
                    str(np.max(numNodesRange))
            plt.title(title, fontsize=14)
            plt.grid()
            plt.legend(loc='best')
        # dis = np.min(thpt4arrow[1:]) - thpt4arrow[0]
        # The codes below highlights the gain through double-arrow and number.
        # plt.annotate(s='',
        #              xy=(numEpochsRange[arrowIdx], thpt4arrow[0] + dis * 0.1),
        #              xytext=(numEpochsRange[arrowIdx],
        #                      thpt4arrow[0] + dis * 0.9),
        #              arrowprops=dict(arrowstyle='<->', color='r', linewidth=1))
        # plt.text(numEpochsRange[arrowIdx] * 0.7, thpt4arrow[0] + dis * 0.5,
        #          str(numShardsRange[-1]) + 'x', fontsize=10)
        plt.tight_layout()
        plt.show()


def plotNode(numShardsRange, numNodesRange, numEpochsRange, data, schemes,
             fileName):
    names = ['full replication', 'uncoded sharding', 'PolyShard']
    markers = ['o', 's', 'd']
    for method in ['Mean']:  # ['Max', 'Median']
        for i in range(len(schemes)):
            s = schemes[i]
            t = data[s]['tVer' + method]
            if s != 'full_replication':
                t[-1, -1] /= 1.05
            thpt = numShardsRange / t[:, -1] / MS
            plt.plot(numNodesRange, thpt, label=names[i], linewidth=2,
                     marker=markers[i], markerfacecolor='w', markersize=10)
            plt.xlabel('Number of nodes ($N$)', fontsize=12)
            plt.ylabel('Throughput (blocks/ms)', fontsize=12)
            title = 'Throughput v.s. System Size when $t$=' + \
                    str(np.max(numEpochsRange))
            plt.title(title, fontsize=14)
            plt.grid()
            plt.legend(loc='best')
        plt.show()


fileName = 'all_schemes_dense_K=[10,50]_M=2000_r=3_epoch=100,500]_s=0.5_1561626.pickle'
plots((fileName))
