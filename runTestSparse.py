import fullReplicationScalingSparse as frscs
import polyShardScalingSparse as psscs


numNodes = 10
numShards = 5
sizeShard = 1000
sparsity = 1
numEpochs = 100
initBal = 100
result = {}
schemes = ['full_replication', 'simple_sharding', 'poly_shard']
result['fr'] = {}
tVerMax, tVerMedian, tVote, tUp = \
    frscs.fullReplicationScalingCore(numNodes, numShards, sizeShard, sparsity,
                                     numEpochs, initBal)
tVerMax, tVerMedian, tVoteMax, tVoteMedian, tUpMax, tUpMean =\
    frscs.simpleShardingScalingCore(numNodes, numShards, sizeShard, sparsity,
                                    numEpochs, initBal)
tVerMax, tVerMedian, tEnMax, tEnMedian, tDe, tUpMax, tUpMedian =\
    psscs.polyShardScalingCore(numNodes, numShards, sizeShard, sparsity,
                               numEpochs, initBal)
