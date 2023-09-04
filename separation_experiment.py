
from experiment import Experiment,selected_kernels,tuple2str
from scipy.stats import spearmanr
from collections import defaultdict
from generate_graphs import (interpolate_ER_triangular,
                             interpolate_ER_PPM,
                             interpolate_ER_GRG_torus,
                             interpolate_ER_inhomogeneous,
                             interpolate_GRG_torus_circle,
                             interpolate_ER_GCG)
from other_kernels import NetLSD,Gin
import numpy as np
import itertools as it
import json

load = True
compute_mmds = False
g_file = 'separation_experiment_graphs.json'
vals_file = 'separation_vals.tsv'
mmds_file = 'separation_mmds.json'

#interpolators = [interpolate_ER_PPM, interpolate_ER_GCG, interpolate_ER_inhomogeneous, interpolate_ER_triangular, interpolate_ER_GRG_torus, interpolate_GRG_torus_circle]

def start_generator(interpolator,name):
    gen = lambda: interpolator(step=0)
    gen.__name__ = name
    return gen

def end_generator(interpolator,name):
    gen = lambda: interpolator(step=1)
    gen.__name__ = name
    return gen

def stat_tests(sample1,sample2):
    from scipy.stats import mannwhitneyu
    try:
        return mannwhitneyu(sample1,sample2)[1]
    except:
        return None

endpoints = [
    start_generator(interpolate_ER_PPM,'ER'),
    end_generator(interpolate_ER_PPM,'PPM'),
    end_generator(interpolate_ER_GCG,'GCG'),
    end_generator(interpolate_ER_inhomogeneous,'inhomogeneous'),
    end_generator(interpolate_ER_triangular,'triangular'),
    end_generator(interpolate_ER_GRG_torus,'torus'),
    end_generator(interpolate_GRG_torus_circle,'circle'),
]


# 50 packs of 100 or 50 graphs (try 100)
# npacks should be a multiple of 3
experiment = Experiment(endpoints,npacks=150,sample_size=50) # Sample size 50 for now
if compute_mmds:
    if load:
        experiment.load_graphs(g_file)
    else:
        experiment.generate_graphs(g_file)
    experiment.apply_kernels(experiment.iterator_compare_generators(),save_name=vals_file,save_mmds_name=mmds_file)

mmds = None
with open(mmds_file) as f:
    mmds = json.load(f)
#mmds = json.loads(mmds_file)
ps = defaultdict(dict)
for e1,e2 in it.combinations(endpoints,2):
    for k in mmds:
        label = e1.__name__+'_vs_'+e2.__name__
        label_null = e1.__name__+'_vs_'+e1.__name__
        ps[k][label] = stat_tests(mmds[k][label],mmds[k][label_null])
with open('separation_ps.json','w') as f:
    json.dump(ps,f)

for k in ps:
    print('{}: {}'.format(k,', '.join([
        gen_pair 
        for gen_pair in ps[k] 
        if ps[k][gen_pair]>0.05
    ])))

