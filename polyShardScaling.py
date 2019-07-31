# V2.1: fixed the bug that codedChain is not updated in codedChainUpdate()
import numpy as np
import pickle
import time
import fullReplicationScaling as frsc
import polyShard as ps


def scalePolyShard(numShards, redundancy, sizeShard, sparsity, numEpochs,
                   initBal):
    numNodes = numShards * redundancy
    scaling = numShards.size

    tVer = np.zeros(scaling)
    tUp = np.zeros(scaling)
    tEn = np.zeros(scaling)
    tDe = np.zeros(scaling)
    throughput = np.zeros(scaling)

    for j in range(scaling):
        tVer[j], tUp[j], tEn[j], tDe[j] = \
            polyShardScalingCore(numNodes[j], numShards[j], sizeShard,
                                 sparsity, numEpochs, initBal)
        throughput[j] = numShards[j] / (tVer[j] + tEn[j] + tDe[j])
    # save the results
    result = {}
    result['verification time'] = tVer
    result['encoding time'] = tEn
    result['decoding time'] = tDe
    result['updating time'] = tUp
    result['throughput'] = throughput
    fileName = 'PolyShard_redundancy_' + str(redundancy) +\
               '_M_' + str(sizeShard) + '_s_' +\
               str(sparsity) + '_numEpochs_' + str(numEpochs) + '_' + \
               str(int(time.time())) + '.pickle'
    with open(fileName, 'wb') as handle:
        pickle.dump(result, handle)

    return tVer, tEn, tDe, tUp, throughput, fileName


def polyShardScalingCore(numNodes, numShards, sizeShard, sparsity, numEpochs,
                         initBal):
    chainLength = numEpochs - 1
    chains = frsc.chainInit(numShards, sizeShard, sparsity, chainLength,
                            initBal)

    beta = [x + 1 for x in range(numShards)]
    alpha = [x + 1 for x in range(numNodes)]

    coefficients = ps.coeGen(numNodes, numShards, beta, alpha)

    codedChains = codedInit(chains, coefficients, sizeShard,
                            numNodes, numShards)

    tVer = np.zeros(numEpochs)
    tUp = np.zeros(numEpochs)
    tEn = np.zeros(numEpochs)
    tDe = np.zeros(numEpochs)

    # print("processing epoch:", numEpochs - 1)
    tVerMax, tVerMedian, tEnMax, tEnMedian, tDe, tUpMax, tUpMedian = \
        ps.simEpoch(codedChains, coefficients, beta, alpha, tVer, tUp, tEn,
                    tDe, numNodes, numShards, sizeShard, sparsity,
                    numEpochs - 1, initBal / numEpochs)

    return tVerMax, tVerMedian, tEnMax, tEnMedian, tDe, tUpMax, tUpMedian


def codedInit(chains, C, sizeShard, numNodes, numShards):
    codedChains = [np.zeros(chains[0].shape) for n in range(numNodes)]
    for n in range(numNodes):
        for s in range(numShards):
            codedChains[n] = codedChains[n] + C[n, s] * chains[s]
    return codedChains
