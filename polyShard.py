# V2.1: fixed the bug that codedChain is not updated in codedChainUpdate()
import numpy as np
import pickle
import time
import fullReplication as fr


def polyShard(numNodes, numShards, sizeShard, sparsity, numEpoches,
              initBal):

    # initialize system
    chains = fr.chainInit(numNodes, numShards, sizeShard, initBal, numEpoches)
    # Generate parameters for Lagrange coding
    beta = [x + 1 for x in range(numShards)]
    alpha = [x + 1 for x in range(numNodes)]

    # coefficients is a numNodes * numShards matrix whose ith row is
    # the coefficients for ith node
    coefficients = coeGen(numNodes, numShards, beta, alpha)

    # generate numNodes coded chains as linear combinations of initial chains
    codedChains = codedInit(chains, coefficients, sizeShard,
                            numNodes, numShards)

    tVer = np.zeros(numEpoches)
    tUp = np.zeros(numEpoches)
    # time to encode incoming blocks before verification
    tEn = np.zeros(numEpoches)
    # time to decode verification results
    tDe = np.zeros(numEpoches)

    # run the system
    for idxEpoch in range(numEpoches):
        print("processing epoch:", idxEpoch)
        simEpoch(codedChains, coefficients, beta, alpha, tVer, tUp, tEn, tDe,
                 numNodes, numShards, sizeShard, sparsity, idxEpoch)

    # save the results
    result = {}
    result['verification time'] = tVer
    result['updating time'] = tUp
    result['encoding time'] = tEn
    result['decoding time'] = tDe
    fileName = 'poly_shard__N_' + str(numNodes) + '_K_' +\
               str(numShards) + '_M_' + str(sizeShard) + '_s_' +\
               str(sparsity) + '_' + str(int(time.time())) + '.pickle'
    with open(fileName, 'wb') as handle:
        pickle.dump(result, handle)

    return tVer, tUp, tEn, tDe, fileName


def coeGen(numNodes, numShards, beta, alpha):
    C = np.ones((numNodes, numShards))

    for i in range(numNodes):  # generate the coefficients for ith node
        for j in range(numShards):
            multiply = list(range(numShards))
            multiply.remove(j)
            for l in multiply:
                C[i][j] = C[i][j] * (alpha[i] - beta[l]) / (beta[j] - beta[l])
    return C


def codedInit(chains, C, sizeShard, numNodes, numShards):
    codedChains = [np.zeros([2, sizeShard]) for n in range(numNodes)]
    for n in range(numNodes):
        for s in range(numShards):
            codedChains[n] = codedChains[n] + C[n, s] * chains[s]
    return codedChains


def simEpoch(codedChains, coefficients, beta, alpha, tVer, tUp, tEn, tDe,
             numNodes, numShards, sizeShard, sparsity, idxEpoch):
    '''
    This function simulates and measures the verification, encoding/decoding,
    and updating happened in each epoch.
    Inputs:
        see above
        idxEpoch: the idx of the current epoch
    Outputs:
        None. It updates codedChains, tVer, tUp, tEncode, and tDecode.
    '''

    # generate blocks
    blocks = fr.blockGen(numShards, sizeShard, sparsity)

    # fetch only the sender vector in the block needed for verification
    senderBlocks = [blocks[k][0, :] for k in range(numShards)]

    senderBlocks_mat = np.zeros((numShards, sizeShard))
    for k in range(numShards):
        senderBlocks_mat[k, :] = senderBlocks[k]

    # generate coded blocks for verification
    codedBlocks = [np.zeros([1, sizeShard]) for k in range(numNodes)]
    tEncode = 0
    for n in range(numNodes):
        tStart = time.time()
        codedBlocks[n] = encode(senderBlocks_mat, coefficients[n])
        tPassed = time.time() - tStart
        tEncode = np.max([tPassed, tEncode])
    tEn[idxEpoch] = tEncode

    # verification
    compResults = np.zeros((numNodes, sizeShard))
    tVerification = 0
    for n in range(numNodes):
        tStart = time.time()
        compResults[n, :] = verifyByNode(codedChains[n], codedBlocks[n], n)
        tPassed = time.time() - tStart
        # verification time across the nodes is the maximum of each node
        tVerification = np.max([tPassed, tVerification])
    tVer[idxEpoch] = tVerification

    # decode
    tStart = time.time()
    results = decode(compResults, numShards)
    decisions = [(results[k, :] >= 0).all() for k in range(numShards)]
    tDe[idxEpoch] = time.time() - tStart

    # update
    receiverBlocks = [blocks[k][1, :] for k in range(numShards)]
    receiverBlocks_mat = np.zeros((numShards, sizeShard))
    for k in range(numShards):
        receiverBlocks_mat[k, :] = receiverBlocks[k]
    decisions_mat = np.array(decisions)
    tUpdate = 0
    for n in range(numNodes):
        tStart = time.time()
        codedChains[n] = codedChainUpdate(codedChains[n], senderBlocks_mat,
                                          receiverBlocks_mat,
                                          decisions_mat, coefficients[n])
        tPassed = time.time() - tStart
        tUpdate = np.max([tPassed, tUpdate])
    tUp[idxEpoch] = tUpdate

    pass


def encode(blocks, coefficient):
    return np.dot(coefficient, blocks)


def verifyByNode(codedChain, codedBlock, n):
    codedBalance = np.sum(codedChain, axis=0)
    codedResult = codedBalance + codedBlock
    temperOpinion(codedResult, n)
    return codedResult


def temperOpinion(codedResult, n):
    '''
    This function decides whether and how node-n tempers its opinion.
    We currently assume all nodes are honest.
    '''
    return codedResult


def decode(compResults, numShards):
    # this happens when the first numShards nodes are honest
    results = compResults[:numShards, :]
    return results


def codedChainUpdate(codedChain, senderBlocks, receiverBlocks, decisions,
                     coefficient):
    codedChain = np.vstack([codedChain,
                            encode(senderBlocks,
                                   np.multiply(decisions, coefficient))])
    codedChain = np.vstack([codedChain,
                            encode(receiverBlocks,
                                   np.multiply(decisions, coefficient))])
    return codedChain
