# (_py37) altsoph@ALTSERVUBN:~/Volume/GraphKernelsBF/20230126$ python apply_netlsd.py > netlsd_raw.tsv

import sys
import json
import netlsd
import networkx as nx
import numpy as np
# from grakel import Graph
# from grakel import GraphKernel
# from grakel.kernels import RandomWalkLabeled, ShortestPathAttr, RandomWalk, PyramidMatch, NeighborhoodHash, ShortestPath, GraphletSampling, SubgraphMatching, WeisfeilerLehman, HadamardCode, NeighborhoodSubgraphPairwiseDistance, SvmTheta, Propagation, PropagationAttr, OddSth, MultiscaleLaplacian, HadamardCode, VertexHistogram, EdgeHistogram, GraphHopper, CoreFramework, WeisfeilerLehmanOptimalAssignment

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
        mtx[gg].append( netlsd.heat(nx.from_edgelist(g)) )
print('matrices done', file=sys.stderr, flush=True)
print(gg, len(mtx[gg]), mtx[gg][0], file=sys.stderr, flush=True)


ks = []
kn = []
for k in ('netlsd',):
    kernel_name = k
    if 1:
#    try:
        for p_idx in range(SAMPLES_NUM):
            print('pack', p_idx, 'done', file=sys.stderr, flush=True)
            gr_pack = []
            for g in gens:
                gr_pack.extend( mtx[g][p_idx*SAMPLE_SIZE:(p_idx+2)*SAMPLE_SIZE] )
            #     # from each generator we take current pack + next pack to fit on
            # vals = k(normalize=True).fit_transform(gr_pack)
            print('pack', p_idx, len(gr_pack), file=sys.stderr, flush=True)
            for i in range(len(gr_pack)):
                for j in range(len(gr_pack)):
                    # print("\t".join(map(str, ['#', kernel_name, p_idx, gens[i//(2*SAMPLE_SIZE)], i%(2*SAMPLE_SIZE), gens[j//(2*SAMPLE_SIZE)], j%(2*SAMPLE_SIZE), vals[i,j]])), flush=True)
                    print("\t".join(map(str, ['#', kernel_name, p_idx, gens[i//(2*SAMPLE_SIZE)], i%(2*SAMPLE_SIZE), gens[j//(2*SAMPLE_SIZE)], j%(2*SAMPLE_SIZE), 
                        netlsd.compare(gr_pack[i], gr_pack[j])])), # vals[i,j]])), 
                        flush=True)
#    except:
#        print('err')
#        pass

