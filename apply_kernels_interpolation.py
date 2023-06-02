'''
Usage:
    python apply_kernels_interpolation.py > kernels_raw.interpolation.tsv

Output format: tsv with 8 columns of the format:
    # kernel_name    interpolator_name    pack_idx    step    idx_in_step    other_step    idx_in_other_step    kernel_value
other_step is either 0 (the beginning of the interpolation) or 1 (the end of the interpolation)
Thus, this column shows the kernel value between d['interpolators'][interpolator_name][step][SAMPLE_SIZE*pack_idx+idx_in_step]
and d['interpolators'][interpolator_name][other_step][SAMPLE_SIZE*pack_idx+idx_in_other_step]
'''
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
fn = 'interpolation_graphs.json'
N = 50
SAMPLE_SIZE = 10
SAMPLES_NUM = 9


d = json.loads(open(fn).read())
print('loaded', file=sys.stderr, flush=True)

mtx = dict() # graphs in matrix form
interpolators = list(d['interpolators'])
steps = list(list(d['interpolators'].values())[0].keys())
first_step = steps[0]
last_step = steps[-1]

for idx, interpolator in enumerate(interpolators):
    mtx[interpolator] = dict()
    for step in steps:
        mtx[interpolator][step] = []
        for g in d['interpolators'][interpolator][step]:
            mtx[interpolator][step].append( Graph(list(map(tuple,g)), node_labels={i: 'A' for i in range(N)}, edge_labels={e: 'B' for e in map(tuple,g)}) )
print('matrices done', file=sys.stderr, flush=True)


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
for k in [ShortestPath,WeisfeilerLehman]:
    '''(
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
    ):'''
    kernel_name = str(k).split("'")[1].split(".")[-1]
    
    try:
        kstart_time = time()
        for interpolator in interpolators:
            for p_idx in range(SAMPLES_NUM):
                for step in steps:
                    print('starting',kernel_name,interpolator, file=sys.stderr, flush=True)
                    pack_first = mtx[interpolator][step][p_idx*SAMPLE_SIZE:(p_idx+2)*SAMPLE_SIZE]
                    pack_first.extend(mtx[interpolator][first_step][p_idx*SAMPLE_SIZE:(p_idx+2)*SAMPLE_SIZE])
                    pack_last = mtx[interpolator][step][p_idx*SAMPLE_SIZE:(p_idx+2)*SAMPLE_SIZE]
                    pack_last.extend(mtx[interpolator][last_step][p_idx*SAMPLE_SIZE:(p_idx+2)*SAMPLE_SIZE])
                    # from each step we take current pack + next pack to fit on
                    pstart_time = time()
                    vals_first = k(normalize=True,**kernel_params[k]).fit_transform(pack_first)
                    vals_last = k(normalize=True,**kernel_params[k]).fit_transform(pack_last)
                    print(k.__name__,interpolator,'step',step,'pack', p_idx,'took',time()-pstart_time,'seconds', file=sys.stderr, flush=True)
                    for i in range(2*SAMPLE_SIZE):
                        for j in range(2*SAMPLE_SIZE,4*SAMPLE_SIZE):
                            # kernel_name    interpolator_name    pack_idx    step    idx_in_step    other_step    idx_in_other_step    kernel_value
                            print("\t".join(map(str, ['#', kernel_name, interpolator, p_idx, step, i, first_step, j-2*SAMPLE_SIZE, vals_first[i,j]])), flush=True)
                    for i in range(2*SAMPLE_SIZE):
                        for j in range(2*SAMPLE_SIZE,4*SAMPLE_SIZE):
                            # kernel_name    interpolator_name    pack_idx    step    idx_in_step    other_step    idx_in_other_step    kernel_value
                            print("\t".join(map(str, ['#', kernel_name, interpolator, p_idx, step, i, last_step, j-2*SAMPLE_SIZE, vals_last[i,j]])), flush=True)
    except:
        print(kernel_name, 'err')
        pass

