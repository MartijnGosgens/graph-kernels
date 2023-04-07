# (_py37) altsoph@ALTSERVUBN:~/Volume/GraphKernelsBF/20230318/GGM-metrics$ python apply_gin.py > gin_raw.tsv

import sys
import json

import networkx as nx
import numpy as np
from scipy import spatial
import torch

import dgl
from evaluation import gin_evaluation


# setup
fn = '1k_samples.json'
N = 50
SAMPLE_SIZE = 30
SAMPLES_NUM = 30
# so we'll use 30x30 = 900 graphs < 1000 we have for each generator

d = json.loads(open(fn).read())
print('loaded', file=sys.stderr, flush=True)

# prepare GIN embedder
device =  torch.device('cpu')
model = gin_evaluation.load_feature_extractor(device)
ev = gin_evaluation.MMDEvaluation(model=model, kernel='rbf', sigma='range', multiplier='mean')
embedder = lambda x: ev._GINMetric__get_activations_single_dataset(x)


# prepare embd data
gens = list(d['graphs']) # list of the generators
mtx = dict() # data container

# for each generator
for idx, gg in enumerate(gens):
    # prepare graphs list in dgl format
    mtx[gg] = list()
    for g in d['graphs'][gg]:
        mtx[gg].append(  dgl.DGLGraph(nx.from_edgelist(g)).to(device)  )
    # convert them to embeddings
    mtx[gg] = embedder(mtx[gg])
print('matrices done', file=sys.stderr, flush=True)
print(gg, len(mtx[gg]), mtx[gg][0], file=sys.stderr, flush=True)

# now count cross distances
kernel_name = 'GIN'
for p_idx in range(SAMPLES_NUM):
    print('pack', p_idx, 'done', file=sys.stderr, flush=True)
    gr_pack = []
    for g in gens:
        gr_pack.extend( mtx[g][p_idx*SAMPLE_SIZE:(p_idx+2)*SAMPLE_SIZE] )
    print('pack', p_idx, len(gr_pack), file=sys.stderr, flush=True)
    for i in range(len(gr_pack)):
        for j in range(len(gr_pack)):
            print("\t".join(map(str, ['#', kernel_name+'.dotprod', p_idx, gens[i//(2*SAMPLE_SIZE)], i%(2*SAMPLE_SIZE), gens[j//(2*SAMPLE_SIZE)], j%(2*SAMPLE_SIZE), 
                np.dot(gr_pack[i], gr_pack[j])])), # vals[i,j]])), 
                flush=True)
            print("\t".join(map(str, ['#', kernel_name+'.cossim', p_idx, gens[i//(2*SAMPLE_SIZE)], i%(2*SAMPLE_SIZE), gens[j//(2*SAMPLE_SIZE)], j%(2*SAMPLE_SIZE), 
                1 - spatial.distance.cosine(gr_pack[i], gr_pack[j]) ])), # vals[i,j]])), 
                flush=True)
