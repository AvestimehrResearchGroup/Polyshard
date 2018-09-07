import fullReplication as fr


numNodes = 100
numShards = 20
sizeShard = 50
sparsity = 0.1
numEpoches = 500
initBal = 10000

# Run full replication
_, _, _, fileName1 = fr.fullReplication(numNodes, numShards, sizeShard,
                                        sparsity, numEpoches, initBal)
fr.plots(fileName1)

# Run simple sharding
_, _, _, fileName2 = fr.simpleSharding(numNodes, numShards, sizeShard,
                                       sparsity, numEpoches, initBal)
fr.plots(fileName2)
