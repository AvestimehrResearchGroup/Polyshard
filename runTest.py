import numpy as np
import pickle
import time
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
MS = 1000


def runTest(density='dense', numShardsRange=list(range(5, 51, 15)),
            numEpochsRange=list(range(100, 1001, 300)),
            redundancy=3, sizeShard=2000, sparsity=0.5,
            initBal=1000, numRuns=2):
    if density == 'dense':
        import fullReplicationScaling as frsc
        import polyShardScaling as pssc
    else:
        import fullReplicationScalingSparse as frsc
        import polyShardScalingSparse as pssc
    numNodesRange = [k * redundancy for k in numShardsRange]

    schemes = ['full_replication', 'simple_sharding', 'poly_shard']
    metrices = ['tVerMax', 'tVerMedian', 'tVoteMax', 'tVoteMedian', 'tEnMax',
                'tEnMedian']
    data = {}
    for s in schemes:
        data[s] = {}
        for m in metrices:
            data[s][m] = np.zeros((len(numShardsRange), len(numEpochsRange)))

    for i in range(len(numShardsRange)):
        numShards = numShardsRange[i]
        numNodes = numNodesRange[i]
        for j in range(len(numEpochsRange)):
            numEpochs = numEpochsRange[j]
            print('Now running K=' + str(numShards) + ' Epoch=' +
                  str(numEpochs))
            for n in range(numRuns):
                # tVerMax, tVerMedian, tVote, tUp = \
                tVerMax, tVerMedian, tVote, _ = \
                    frsc.fullReplicationScalingCore(numNodes, numShards,
                                                    sizeShard, sparsity,
                                                    numEpochs, initBal)
                data['full_replication']['tVerMax'][i, j] += tVerMax / numRuns
                data['full_replication']['tVerMedian'][i, j] += \
                    tVerMedian / numRuns
                data['full_replication']['tVoteMax'][i, j] += tVote / numRuns
                data['full_replication']['tVoteMedian'][i, j] += \
                    tVote / numRuns

            # tVerMax, tVerMedian, tVoteMax, tVoteMedian, tUpMax, tUpMean = \
                tVerMax, tVerMedian, tVoteMax, tVoteMedian, _, _ =\
                    frsc.simpleShardingScalingCore(numNodes, numShards,
                                                   sizeShard, sparsity,
                                                   numEpochs, initBal)

                data['simple_sharding']['tVerMax'][i, j] += tVerMax / numRuns
                data['simple_sharding']['tVerMedian'][i, j] += \
                    tVerMedian / numRuns
                data['simple_sharding']['tVoteMax'][i, j] += tVote / numRuns
                data['simple_sharding']['tVoteMedian'][i, j] += tVote / numRuns
            # tVerMax, tVerMedian, tEnMax, tEnMedian, tDe, tUpMax, tUpMedian =\
                tVerMax, tVerMedian, tEnMax, tEnMedian, _, _, _ =\
                    pssc.polyShardScalingCore(numNodes, numShards, sizeShard,
                                              sparsity, numEpochs, initBal)
                data['poly_shard']['tVerMax'][i, j] += tVerMax / numRuns
                data['poly_shard']['tVerMedian'][i, j] += tVerMedian / numRuns
                data['poly_shard']['tEnMax'][i, j] += tVote / numRuns
                data['poly_shard']['tEnMedian'][i, j] += tVote / numRuns
        result = {}
        result['numShardsRange'] = numShardsRange
        result['numNodesRange'] = numNodesRange
        result['numEpochsRange'] = numEpochsRange
        result['data'] = data
        fileName = 'all_schemes_' + density + '_' + \
                   'K=[' + str(numShards) + '_' \
                   'M=' + str(sizeShard) + '_' + \
                   'r=' + str(redundancy) + '_' + \
                   'epoch=' + str(numEpochsRange[0]) + ',' + \
                   str(numEpochsRange[-1]) + ']_' + \
                   's=' + str(sparsity) + '_' + \
                   str(int(time.time() / 1000)) + '.pickle'
        with open(fileName, 'wb') as handle:
            pickle.dump(result, handle)
    result = {}
    result['numShardsRange'] = numShardsRange
    result['numNodesRange'] = numNodesRange
    result['numEpochsRange'] = numEpochsRange
    result['data'] = data
    fileName = 'all_schemes_' + density + '_' + \
               'K=[' + str(numShardsRange[0]) + ',' + \
               str(numShardsRange[-1]) + ']_' + \
               'M=' + str(sizeShard) + '_' + \
               'r=' + str(redundancy) + '_' + \
               'epoch=' + str(numEpochsRange[0]) + ',' + \
               str(numEpochsRange[-1]) + ']_' + \
               's=' + str(sparsity) + '_' + \
               str(int(time.time() / 1000)) + '.pickle'
    with open(fileName, 'wb') as handle:
        pickle.dump(result, handle)
    return fileName


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
    for method in ['Max', 'Median']:
        t = data[scheme]['tVer' + method] + 2 * data[scheme]['tEn' + method]
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
    for method in ['Max', 'Median']:
        for s in schemes:
            t = data[s]['tVer' + method] + \
                2 * data[s]['tEn' + method]
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
    for method in ['Max', 'Median']:
        for s in schemes:
            t = data[s]['tVer' + method] + \
                2 * data[s]['tEn' + method]
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


density = 'dense'
sparsity = 0.5
fileName = runTest(density=density, sparsity=sparsity)
plots(fileName)
