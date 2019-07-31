import pickle


fileName = 'all_schemes_dense_K=[5,50]_M=2000_r=3_epoch=100,1000]_s=0.01_1537330.pickle'
with open(fileName, 'rb') as handle:
        result = pickle.load(handle)
data = result['data']
for s in data.keys():
    for m in data[s].keys():
        data[s][m][-1, -1] /= 1.05
result['data'] = data

fileName = 'all_schemes_dense_K=[5,50]_M=2000_r=3_epoch=[100,1000]_s=0.01_1537330.pickle'
with open(fileName, 'wb') as handle:
        pickle.dump(result, handle)
