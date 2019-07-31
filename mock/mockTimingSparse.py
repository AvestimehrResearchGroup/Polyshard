import numpy as np
import scipy.sparse as scs
import time


# Block chain bulit on top of sparse matrix
def frEpoch(numShards, numNodes, sizeShard, sparsity, chainLength, initBal):
    initChain = np.vstack([np.zeros((1, sizeShard)),
                           initBal * np.ones((1, sizeShard))])
    block = blockGenCore(sizeShard, sparsity, initBal / chainLength / 2)
    chain = np.vstack([initChain, np.repeat(block, chainLength, axis=0)])

    chain = scs.coo_matrix(chain)
    senderBlock = scs.coo_matrix(block[0, :])

    tVer = np.zeros(numNodes)
    for n in range(numNodes):
        # measure the verification time of each node
        for k in range(numShards):
            tVer[n] += verifyCore(chain, senderBlock)
    return np.max(tVer), np.median(tVer), np.mean(tVer)


def ssEpoch(numShards, numNodes, sizeShard, sparsity, chainLength, initBal):
    tVerMax = []
    tVerMedian = []
    tVerMean = []
    numRep = int(numNodes / numShards)
    for k in range(numShards):
        tMax, tMedian, tMean = frEpoch(1, numRep, sizeShard, sparsity,
                                       chainLength, initBal)
        tVerMax.append(tMax)
        tVerMedian.append(tMedian)
        tVerMean.append(tMean)
    return np.max(tVerMax), np.median(tVerMedian), np.mean(tVerMean)


def rowwiseShuffle(arr):
    # shuffle each row of a matrix independently
    x, y = arr.shape
    rows = np.indices((x, y))[0]
    cols = [np.random.permutation(y) for _ in range(x)]
    return arr[rows, cols]


def psEpoch(numShards, numNodes, sizeShard, sparsity, chainLength, initBal):
    initChain = np.vstack([np.zeros((1, sizeShard)),
                           initBal * np.ones((1, sizeShard))])
    block = blockGenCore(sizeShard, sparsity, initBal / chainLength / 2)
    blocks = np.tile(block, (chainLength, 1))
    blocks = rowwiseShuffle(blocks)

    chain = np.vstack([initChain, blocks])

    senderBlocks = blocks[0:numShards * 2:2, ]
    senderBlocks = rowwiseShuffle(senderBlocks)
    beta = np.array(range(numShards)) + 1
    alpha = np.array(range(numNodes)) + 1

    coeff = coeGen(numNodes, numShards, beta, alpha)
    codedChain = chainLength
    for i in range(1, numShards):
        codedChain += coeff[-1, i] * rowwiseShuffle(chain)

    codedChain = scs.coo_matrix(codedChain)
    senderBlocks = scs.coo_matrix(senderBlocks)

    tVer = []
    for n in range(numShards, numNodes):  # the first numShards nodes are fast.
        # measure the verification time of each node
        start = time.time()
        # run twice to mimic encoding + decoding
        codedSenderBlock = scs.coo_matrix.dot(coeff[n, :], senderBlocks)
        codedSenderBlock = scs.coo_matrix.dot(coeff[n, :], senderBlocks)
        tVer.append(time.time() - start)
        tVer[-1] += verifyCore(codedChain, codedSenderBlock)
    return np.max(tVer), np.median(tVer), np.mean(tVer)


def blockGenCore(sizeShard, sparsity, txCap):
    '''
    This function creates a block that contains sizeShard * sparisity Trans.
    Inputs:
        see above
    Output:
        block, see above
    '''
    numTrans = int(sizeShard * sparsity)
    userShuffle = np.random.permutation(sizeShard)
    idxSenders = userShuffle[:numTrans]
    userShuffle = np.random.permutation(sizeShard)
    idxRecivers = userShuffle[:numTrans]
    block = np.zeros((2, sizeShard))
    block[0, idxSenders] = -txCap
    block[1, idxRecivers] = txCap
    return block


def coeGen(numNodes, numShards, beta, alpha):
    C = np.ones((numNodes, numShards))

    for i in range(numNodes):  # generate the coefficients for ith node
        for j in range(numShards):
            multiply = list(range(numShards))
            multiply.remove(j)
            for l in multiply:
                C[i][j] = C[i][j] * (alpha[i] - beta[l]) / (beta[j] - beta[l])
    return C


def verifyCore(chain, block):
    start = time.time()
    bal = np.sum(chain, axis=0)
    newBal = bal + block
    ignore = (newBal > 0).all()
    t = time.time() - start
    return t

# def oneRound(numShards, numNodes, sizeShard, sparsity, chainLength, initBal):
#     tVerMax, tVerMedian = frEpoch(numShards, numNodes, sizeShard, sparsity,
#                                   chainLength, initBal)
#     tVerMax, tVerMedian = ssEpoch(numShards, numNodes, sizeShard, sparsity,
#                                   chainLength, initBal)
#     tVerMax, tVerMedian = psEpoch(numShards, numNodes, sizeShard, sparsity,
#                                   chainLength, initBal)


# numShards = 10
# numNodes = 30
# sizeShard = 100
# sparsity = 0.5
# chainLength = 100
# initBal = 1000
# oneRound(numShards, numNodes, sizeShard, sparsity, chainLength, initBal)
