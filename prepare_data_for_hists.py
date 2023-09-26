import json

from glob import glob
from collections import defaultdict, Counter

import networkx as nx
import numpy as np


def get_data(fn):
    return json.loads(open(fn).read())

fns = [fn for fn in glob('interpolation_graphs/*')]

def process_graph(el):
    G = nx.Graph(el)
    n = len(G.nodes)
    degree_cnt = Counter([d for n, d in G.degree()])
    sp_cnt = Counter()
    p = dict(nx.shortest_path_length(G))
    for i in range(n):
        for j in range(n):
            # if i==j: continue
            if i in p and j in p[i]:
                sp_cnt[p[i][j]]+=1
            else:
                sp_cnt[-1]+=1
    return degree_cnt, sp_cnt


# from matplotlib import pyplot as plt

def process_transition(fn):
    data = get_data(fn)
    for k in data:
        print(k)
        degrees = []
        spaths = []
        _, k1, k2 = k.replace('_torus', '-torus').split('_')
        print(k1,'->',k2)
        for kk, v in data[k].items():
            print(kk)
            d_cnt = Counter()
            p_cnt = Counter()
            for idx, pack in enumerate(v):
                print(idx, end=' ')
                for g in pack:
                    d, p = process_graph(g)
                    d_cnt.update(d)
                    p_cnt.update(p)
            degrees.append( d_cnt )
            spaths.append( p_cnt )

        with open(f'graphstats/{k}_undir_degrees.json','w',encoding='utf-8') as ofh:
            print(json.dumps(degrees, indent=2), file=ofh)
        with open(f'graphstats/{k}_undir_spaths.json','w',encoding='utf-8') as ofh:
            print(json.dumps(spaths, indent=2), file=ofh)    
    
        # draw_endpoints(degrees, 'Degrees', k1, k2)
        # draw_endpoints(spaths, 'Shortest paths', k1, k2)
    return degrees, spaths

for tr_fn in fns:
	_degrees, _spaths = process_transition(tr_fn)
