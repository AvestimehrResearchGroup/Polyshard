import fullReplicationScaling as frs
import polyShardScaling as pss
import time
import numpy as np

numShards = np.arange(10, 105, 10)
redundancy = 3
sizeShard = 50
sparsity = 0.2
numEpochs = 2000
initBal = 100


# Run full replication
# start = time.time()
# out1 = frs.scaleFullReplication(numShards, redundancy, sizeShard, sparsity, numEpochs, initBal)
# print("FR: ", time.time() - start)
# frs.plots(out1[-1],numShards * redundancy)


# Run simple sharding
start = time.time()
out2 = frs.scaleSimpleSharding(numShards, redundancy, sizeShard, sparsity, numEpochs, initBal)
print("SS: ", time.time() - start)
frs.plots(out2[-1],numShards * redundancy)


# Run polynomial sharding
start = time.time()
out3 = pss.scalePolyShard(numShards, redundancy, sizeShard, sparsity, numEpochs, initBal)
print("PS: ", time.time() - start)
frs.plots(out3[-1],numShards * redundancy)
