import numpy as np
import scipy.sparse as ss
import pickle
import time
import matplotlib.pyplot as plt
import fullReplicationSparse as frs


def scaleFullReplication(numShards, redundancy, sizeShard, sparsity,
                         numEpochs, initBal):
    numNodes = numShards * redundancy
    scaling = numShards.size

    tVer = np.zeros(scaling)
    tUp = np.zeros(scaling)
    tVote = np.zeros(scaling)
    throughput = np.zeros(scaling)

    for j in range(scaling):
        tVer[j], tVote[j], tUp[j] = \
            fullReplicationScalingCore(numNodes[j], numShards[j], sizeShard,
                                       sparsity, numEpochs, initBal)
        throughput[j] = numShards[j] / (tVer[j] + tVote[j])

    # save the results
    result = {}
    result['verification time'] = tVer
    result['voting time'] = tVote
    result['updating time'] = tUp
    result['throughput'] = throughput
    fileName = 'full_replication_redundancy_' + str(redundancy) +\
               '_M_' + str(sizeShard) + '_s_' +\
               str(sparsity) + '_numEpochs_' + str(numEpochs) + '_' +\
               str(int(time.time())) + '.pickle'
    with open(fileName, 'wb') as handle:
        pickle.dump(result, handle)

    return tVer, tVote, tUp, throughput, fileName


def fullReplicationScalingCore(numNodes, numShards, sizeShard, sparsity,
                               numEpochs, initBal):
    '''
    This function simulates the verification process of a full replication
    scheme in a single epoch, with numEpochs blocks verified in each chain
    '''

    # initialize the system after numEpochs epoches
    chainLength = numEpochs - 1
    chains = chainInit(numShards, sizeShard, sparsity, chainLength, initBal)
    tVer = np.zeros(numEpochs)
    tUp = np.zeros(numEpochs)
    tVote = np.zeros(numEpochs)

    tVerMax, tVerMedian, tVote, tUp = \
        frs.simEpoch(chains, tVer, tUp, tVote, numNodes, numShards,
                     sizeShard, sparsity, numEpochs - 1, initBal / numEpochs)

    return tVerMax, tVerMedian, tVote, tUp


def chainInit(numShards, sizeShard, sparsity, chainLength, initBal):
    chains = frs.chainInit(numShards, sizeShard, initBal)
    txCap = initBal / chainLength / 2
    for k in range(numShards):
        for e in range(chainLength):
            chains[k] = ss.vstack([chains[k],
                                   frs.blockGenCore(sizeShard, sparsity,
                                                    txCap)])
    return chains


def scaleSimpleSharding(numShards, redundancy, sizeShard, sparsity, numEpochs,
                        initBal):
    numNodes = numShards * redundancy
    scaling = numShards.size

    tVer = np.zeros(scaling)
    tUp = np.zeros(scaling)
    tVote = np.zeros(scaling)
    throughput = np.zeros(scaling)

    for j in range(scaling):
        tVer[j], tVote[j], tUp[j] = \
            simpleShardingScalingCore(numNodes[j], numShards[j], sizeShard,
                                      sparsity, numEpochs, initBal)
        throughput[j] = numShards[j] / (tVer[j] + tVote[j])

    # save the results
    result = {}
    result['verification time'] = tVer
    result['voting time'] = tVote
    result['updating time'] = tUp
    result['throughput'] = throughput
    fileName = 'simple_sharding_redundancy_' + str(redundancy) +\
               '_M_' + str(sizeShard) + '_s_' +\
               str(sparsity) + '_numEpochs_' + str(numEpochs) + '_' +\
               str(int(time.time())) + '.pickle'
    with open(fileName, 'wb') as handle:
        pickle.dump(result, handle)

    return tVer, tVote, tUp, throughput, fileName


def simpleShardingScalingCore(numNodes, numShards, sizeShard, sparsity,
                              numEpochs, initBal):
    assert numNodes % numShards == 0

    # the number of nodes that repeats a shard.
    numRep = int(numNodes / numShards)
    tVerMax = []
    tVerMedian = []
    tVote = []
    tUp = []
    for k in range(numShards):
        t1, t2, t3, t4 = \
            fullReplicationScalingCore(numNodes=numRep, numShards=1,
                                       sizeShard=sizeShard, sparsity=sparsity,
                                       numEpochs=numEpochs, initBal=initBal)

        tVerMax.append(t1)
        tVerMedian.append(t2)
        tVote.append(t3)
        tUp.append(t4)
    return np.max(tVerMax), np.median(tVerMedian), np.max(tVote), \
        np.median(tVote), np.max(tUp), np.mean(tUp)


def plots(fileName, numNodes):
    print(fileName)
    with open(fileName, 'rb') as handle:
        result = pickle.load(handle)
    for key in result.keys():
        if key == 'throughput':
            plt.plot(numNodes, result[key], label=key)
    plt.xlabel('Number of nodes')
    plt.ylabel('Throughput (# of blocks per sec)')
    plt.grid()
    plt.legend(loc='best')
    plt.title('Data source:\n' + fileName)
    plt.xlim((0, None))
    plt.ylim((0, None))
    plt.tight_layout()
    plt.show()
