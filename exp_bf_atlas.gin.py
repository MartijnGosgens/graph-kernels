import json
import itertools as it
from collections import defaultdict
from time import time

import numpy as np
import networkx as nx

from generate_graphs import ig2edges,edges2grakel
from grakel.kernels import (RandomWalk,
                            GraphletSampling,
                            PyramidMatch,
                            NeighborhoodHash,
                            ShortestPath,
                            WeisfeilerLehman,
                            Propagation,
                            OddSth,
                            WeisfeilerLehmanOptimalAssignment,
                            NeighborhoodSubgraphPairwiseDistance)
from other_kernels import NetLSD,Gin
from grakel import Graph
from grakel.utils import graph_from_networkx

selected_kernels = (
     Gin,
    )


slices = defaultdict(set)
names = []

gk_graphs = []
for idx, g in enumerate(nx.graph_atlas_g()):
    if not idx:
      continue
    if not len(g.edges()):
      continue
    v_n = len(g.nodes())
    if not v_n:
      continue
    conn = nx.is_connected(g)
    gn = str(g).split("'")[1]
    names.append( gn )
    slices['full'].add(idx)
    slices[f'size_{v_n}'].add(idx)
    if conn:
        slices['connected'].add(idx)
        slices[f'connected_size_{v_n}'].add(idx)
    else:
        slices['non_connected'].add(idx)
        slices[f'non_connected_size_{v_n}'].add(idx)

    nx.set_node_attributes(g, dict([(i,'A') for i in range(v_n)]), 'label')
    nx.set_edge_attributes(g, {tuple(ee): {'elabel':'B'} for ee in g.edges()})
    gg = list(graph_from_networkx([g,], node_labels_tag='label', edge_labels_tag='elabel', as_Graph=True))[0] # , as_Graph=True)
    gk_graphs.append( gg )

kernel_params = defaultdict(dict)
kernel_params[GraphletSampling] = {
    'sampling': {'n_samples': 500}
}
kernel_params[RandomWalk] = {
    'lamda': 0.1, # not 'lambda' (typo in grakel?)
    'p': 5
}

kernel_vals = dict()
for k in selected_kernels:
    print(k.__name__)
    vals = k(normalize=True,**kernel_params[k]).fit_transform(gk_graphs)
    print((vals == 1.).sum(), len(gk_graphs))
    kernel_vals[k.__name__] = vals

pairs = defaultdict(list)
for k in kernel_vals:
  for i, r in enumerate(kernel_vals[k].tolist()):
    for j, c in enumerate(r):
      if i>j and c == 1.:
        pairs[k].append((names[i],names[j]))
with open('pairs_gin.json', 'w') as ofh:
  print(json.dumps(pairs), file=ofh)

print(pairs('Gin'))