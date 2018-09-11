import fullReplication as fr
import polyShard as ps
import time


numNodes = 10
numShards = 5
sizeShard = 1000
sparsity = 0.01
numEpoxhs = 100
initBal = 100

# Run full replication
start = time.time()
out1 = fr.fullReplication(numNodes, numShards, sizeShard,
                          sparsity, numEpoxhs, initBal)
print("FR: ", time.time() - start)
fr.plots(out1[-1])


# Run simple sharding
start = time.time()
out2 = fr.simpleSharding(numNodes, numShards, sizeShard,
                         sparsity, numEpoxhs, initBal)
print("SS: ", time.time() - start)
fr.plots(out2[-1])

# Run polynomial sharding
start = time.time()
out3 = ps.polyShard(numNodes, numShards, sizeShard, sparsity,
                    numEpoxhs, initBal)
print("PS: ", time.time() - start)
fr.plots(out3[-1])
