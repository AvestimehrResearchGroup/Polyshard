import fullReplicationSparse as frs
import polyShardSparse as pss
import time


numNodes = 10
numShards = 5
sizeShard = 1000
sparsity = 1
numEpoches = 1000
initBal = 100

# Run full replication
start = time.time()
out1 = frs.fullReplication(numNodes, numShards, sizeShard,
                           sparsity, numEpoches, initBal)
print("FRS: ", time.time() - start)
frs.plots(out1[-1])


# # Run simple sharding
start = time.time()
out2 = frs.simpleSharding(numNodes, numShards, sizeShard,
                          sparsity, numEpoches, initBal)
print("SSS: ", time.time() - start)
frs.plots(out2[-1])

# Run polynomial sharding
start = time.time()
out3 = pss.polyShard(numNodes, numShards, sizeShard, sparsity,
                     numEpoches, initBal)
print("PSS: ", time.time() - start)
frs.plots(out3[-1])
