# python apply_kernels_v2.py > kernels_raw.v2.tsv

import sys
import json
import networkx as nx
import numpy as np
from grakel import Graph
from grakel import GraphKernel
from grakel.kernels import RandomWalkLabeled, ShortestPathAttr, RandomWalk, PyramidMatch, NeighborhoodHash, ShortestPath, GraphletSampling, SubgraphMatching, WeisfeilerLehman, HadamardCode, NeighborhoodSubgraphPairwiseDistance, SvmTheta, Propagation, PropagationAttr, OddSth, MultiscaleLaplacian, HadamardCode, VertexHistogram, EdgeHistogram, GraphHopper, CoreFramework, WeisfeilerLehmanOptimalAssignment
from collections import defaultdict
from time import time

# setup
fn = '1k_samples.json'
N = 50
SAMPLE_SIZE = 30
SAMPLES_NUM = 30
# so we'll use 30x30 = 900 graphs < 1000 we have for each generator


d = json.loads(open(fn).read())
print('loaded', file=sys.stderr, flush=True)

mtx = dict() # graphs in matrix form
gens = list(d['graphs']) # list of the generators

for idx, gg in enumerate(gens):
    mtx[gg] = list()
    for g in d['graphs'][gg]:
        mtx[gg].append( Graph(list(map(tuple,g)), node_labels={i: 'A' for i in range(N)}, edge_labels={e: 'B' for e in map(tuple,g)}) )
print('matrices done', file=sys.stderr, flush=True)
print(gg, len(mtx[gg]), mtx[gg][0], file=sys.stderr, flush=True)


kernel_params = defaultdict(dict)
kernel_params[GraphletSampling] = {
    'sampling': {'n_samples': 5000}
}
kernel_params[RandomWalk] = {
    'lamda': 0.1, # not 'lambda' (typo in grakel?)
    'p': 5
}

ks = []
kn = []
for k in (
     #RandomWalk, # ERRRs,
     GraphletSampling,
     PyramidMatch,
     NeighborhoodHash,
     ShortestPath,
     WeisfeilerLehman,
     Propagation,
     OddSth,
     WeisfeilerLehmanOptimalAssignment,
     NeighborhoodSubgraphPairwiseDistance,
     # SvmTheta, # ERRR
    ):

    kernel_name = str(k).split("'")[1].split(".")[-1]
    try:
        kstart_time = time()
        for p_idx in range(SAMPLES_NUM):
            gr_pack = []
            for g in gens:
                gr_pack.extend( mtx[g][p_idx*SAMPLE_SIZE:(p_idx+2)*SAMPLE_SIZE] )
                # from each generator we take current pack + next pack to fit on
            pstart_time = time()
            vals = k(normalize=True,**kernel_params[k]).fit_transform(gr_pack)
            print(k.__name__,'pack', p_idx,'took',time()-pstart_time,'seconds', file=sys.stderr, flush=True)
            for i in range(len(vals)):
                for j in range(len(vals)):
                    print("\t".join(map(str, ['#', kernel_name, p_idx, gens[i//(2*SAMPLE_SIZE)], i%(2*SAMPLE_SIZE), gens[j//(2*SAMPLE_SIZE)], j%(2*SAMPLE_SIZE), vals[i,j]])), flush=True)
        print(k.__name__,'took',time()-kstart_time,'seconds in total', file=sys.stderr, flush=True)
    except:
        print(kernel_name, 'err')
        pass

