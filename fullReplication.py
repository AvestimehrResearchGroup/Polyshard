# V1.0
import numpy as np
import pickle
import time
import matplotlib.pyplot as plt


def fullReplication(numNodes, numShards, sizeShard, sparsity, numEpoches,
                    initBal):
    '''
    This function simulates the verification, voting, and updating time of
    full replication block chain.
    Inputs:
        numNodes: the number of nodes in the system. Every node stores all the
                  transactions.
        numShards: the number of shards the users are partitioned into.
                   Users from different shards cannot do TXs with each other.
        sizeShard: the number of users in each shard.
        sparsity: the number of transactions per block per shard as a fraction
                  of the number of users (sizeShard). We assume every user
                  makes at most one transaction per block.
        numEpoch: the number of epochs we simulate. In each epoch, every shard
                  submits a new block.
        initFund: the initial balance of each user
    Outputs:
        tVer: a vector of length numEpoch recording the verification time of
              each epoch.
        tVote: a vector of length numEpoch recording the voting time of
              each epoch.
        tUp: a vector of length numEpoch recording the chain logging time of
              each epoch.
        fileName: the name of the pickle file that stores the above 3 results.
    '''

    # initialize the system
    chains = chainInit(numNodes, numShards, sizeShard, initBal, numEpoches)
    tVer = np.zeros(numEpoches)
    tUp = np.zeros(numEpoches)
    tVote = np.zeros(numEpoches)

    # run the system
    for idxEpoch in range(numEpoches):
        print("processing epoch:", idxEpoch)
        simEpoch(chains, tVer, tUp, tVote, numNodes, numShards,
                 sizeShard, sparsity, idxEpoch)

    # save the results
    result = {}
    result['verification time'] = tVer
    result['voting time'] = tUp
    result['updating time'] = tVote
    fileName = 'full_replication_N_' + str(numNodes) + '_K_' +\
               str(numShards) + '_M_' + str(sizeShard) + '_s_' +\
               str(sparsity) + '_' + str(int(time.time())) + '.pickle'
    with open(fileName, 'wb') as handle:
        pickle.dump(result, handle)

    return tVer, tVote, tUp, fileName


def chainInit(numNodes, numShards, sizeShard, initBal, numEpoch):
    '''
    This function creates the numShards chains stored at each node.
    Inputs:
        see above
    Output:
        chains: a list of numShards chains.

                chain: Each chain is a matrix with sizeShard columns.
                Every two rows constitute a block.
                The first row records the balance the senders will pay.
                The second row records the balance the receivers will receive.
                Its number of rows increases over time to up to numEpoch * 2.

                NB: In practice we should create a separate copy of chains
                for each node. But due to Python memory limitation,
                we only create one. This is fine for simulation because all
                the nodes maintain the same chains.
    '''
    chain = np.zeros([2, sizeShard])
    chain[1, :] = initBal  # the 0th + 1st row constitute a genesis block
    chains = [np.copy(chain) for k in range(numShards)]
    return chains


def simEpoch(chains, tVer, tUp, tVote, numNodes, numShards,
             sizeShard, sparsity, idxEpoch):
    '''
    This function simulates and measures the verification, voting, and updating
    happened in each epoch.
    Inputs:
        see above
        idxEpoch: the idx of the current epoch
    Outputs:
        None. It updates chains, tVer, tLog, and tVote.
    '''

    # generate blocks
    blocks = blockGen(numShards, sizeShard, sparsity)

    # verification
    # each of the numNodes nodes' opinion on each of the numShards blocks
    opinions = np.zeros((numNodes, numShards))
    tVerification = 0
    for n in range(numNodes):
        tStart = time.time()
        opinions[n, :] = verifyByNode(chains, blocks, n)
        tPassed = time.time() - tStart
        # verification time across the nodes is the maximum of each node
        tVerification = np.max([tPassed, tVerification])
    tVer[idxEpoch] = tVerification

    # voting
    tStart = time.time()
    decisions = voting(opinions)
    tVote[idxEpoch] = time.time() - tStart

    # update
    tStart = time.time()
    chainUpdate(chains, blocks, decisions)
    tUp[idxEpoch] = time.time() - tStart
    pass


def blockGen(numShards, sizeShard, sparsity):
    '''
    This function generates numShards blocks
    Inputs:
        see above
    Output:
        blocks: a list of numShards blocks.

                block: Each block is a 2 * sizeShard matrix.
                The 1st vector of a block shows the amount of money all the
                sizeShard users will send. All entries are non-positive.
                The 2nd vector of a block shows the amount of money all the
                sizeShard users will receive. All entries are non-positive.
                Every user can both send and receive money.
    '''
    blocks = [blockGenCore(sizeShard, sparsity) for n in range(numShards)]
    return blocks


def blockGenCore(sizeShard, sparsity):
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
    amount = np.random.randint(1, 100, size=numTrans)
    block = np.zeros((2, sizeShard))
    block[0, idxSenders] = -amount
    block[1, idxRecivers] = amount
    return block


def verifyByNode(chains, blocks, n):
    '''
    This function verifies all blocks against its respective chain based on
    balance check at node-n.
    Inputs:
        chains: see above
        blocks: see above
        n: the idx of the node which perform this verification.
           This is a place holder for futher implementation of malicious nodes.
    Output:
        opinion: a list of numShards True/False opinions
    '''
    numShards = len(chains)
    opinion = np.zeros(numShards)
    for k in range(numShards):
        opinion[k] = verifyBlock(chains[k], blocks[k])
    temperOpinion(opinion, n)  # temper the opinion
    return opinion


def verifyBlock(chain, block):
    '''
    This function verifies a block against its chain based on the balance
    Inputs:
        see above
    output:
        a True/False decision
    '''
    balance = np.sum(chain, axis=0)
    return (balance + block[0] >= 0).all()


def temperOpinion(opinion, n):
    '''
    This function decies whetehr and how node-n tempers its opinion.
    We currently assume all nodes are honest.
    '''
    return opinion


def voting(opinions):
    '''
    This function makes acceptance decisions for all the numShards blocks.
    Input:
        opinions: a numNodes * numShards binary matrix.
                  The [n, k]-th entry is the opinion of node-n on block-k.
    Output:
        decisions: a length-numShards binary vector.
                   The k-th entry is the final decision on block-k.
    '''
    numNodes, _ = opinions.shape
    votes = np.sum(opinions, axis=0)
    decisions = votes > numNodes / 2  # pass if more than half vote up
    return decisions


def chainUpdate(chains, blocks, decisions):  # incomplete
    '''
    This function updates each chain using the respective block based on the
    respective decision.
    Inputs:
        see above
    Outputs:
        None. It updates chains.
    '''
    for i in range(len(chains)):
        chains[i] = np.vstack([chains[i], blocks[i] * decisions[i]])


def plots(fileName):
    print(fileName)
    with open(fileName, 'rb') as handle:
        result = pickle.load(handle)
    for key in result.keys():
        plt.plot(range(len(result[key])), result[key] * 1000, label=key)
    plt.xlabel('Epoch index')
    plt.ylabel('Time (ms)')
    plt.grid()
    plt.legend(loc='best')
    plt.title('Data source:\n' + fileName)
    plt.xlim((0, None))
    plt.ylim((0, None))
    plt.tight_layout()
    plt.show()


def uniTests():
    # test verifyBlock()
    chain = np.matrix([[0, 0, 0], [100, 100, 100]])
    block = np.matrix([[-1, -2, -3], [3, 2, 1]])
    assert verifyBlock(chain, block)
    block = np.matrix([[-101, -2, -3], [3, 2, 1]])
    assert not verifyBlock(chain, block)

    # test voting()
    opinions = np.matrix([[True, False, False],
                          [True, True, True],
                          [False, False, False],
                          [True, True, False]])
    assert (voting(opinions) == [True, False, False]).all()

    # test chainUpdate()
    chains = [np.matrix([[0, 0, 0], [100, 100, 100]]) for i in range(2)]
    blocks = [np.matrix([[-1, -2, -3], [3, 2, 1]]) for i in range(2)]
    decisions = [True, False]
    chainUpdate(chains, blocks, decisions)
    assert (chains[0][2, :] == [-1, -2, -3]).all()
    assert (chains[0][3, :] == [3, 2, 1]).all()
    assert (chains[1][2, :] == [0, 0, 0]).all()
    assert (chains[1][3, :] == [0, 0, 0]).all()
    print('All unit tests passed!')


def simpleSharding(numNodes, numShards, sizeShard, sparsity,
                   numEpoches, initBal):
    '''
    This function simulates the verification, voting, and updating time of
    block chain with simple sharding. In this system, each shard is repeated
    by a cluster of numNodes / numShards times. The handling of the numShards
    shards are independent of each otehr.

    Simulation method:
    In such a simple sharding system, each shard/cluster is indeed equivalently
    to a small full replication system with numNodes / numShards nodes and
    only one shard. Thus, it is sufficient to simulate simple sharding by
    running small a full replication system numShards times.
    Inputs:
        see above
    Outputs:
        tVer: see above
        tVote: see above
        tUp: see above
        fileName: see above
    '''

    # numNodes must be a multiple of numShards
    assert numNodes % numShards == 0

    # the number of nodes that repeats a shard.
    numRep = int(numNodes / numShards)
    tVer = []
    tVote = []
    tUp = []
    for k in range(numShards):
        print("processing shard:", str(k))
        t1, t2, t3, _ = fullReplication(numNodes=numRep, numShards=1,
                                        sizeShard=sizeShard, sparsity=sparsity,
                                        numEpoches=numEpoches, initBal=initBal)

        tVer.append(t1)
        tVote.append(t2)
        tUp.append(t3)

    # each time at each epoch should be the maximum across all shards.
    tVer = np.max(tVer, axis=0)
    tVote = np.max(tVote, axis=0)
    tUp = np.max(tUp, axis=0)
    result = {}
    result['verification time'] = tVer
    result['voting time'] = tUp
    result['updating time'] = tVote
    fileName = 'simple_sharding_N_' + str(numNodes) + '_K_' +\
               str(numShards) + '_M_' + str(sizeShard) + '_s_' +\
               str(sparsity) + '_' + str(int(time.time())) + '.pickle'
    with open(fileName, 'wb') as handle:
        pickle.dump(result, handle)
    return tVer, tVote, tUp, fileName


# uniTests()
