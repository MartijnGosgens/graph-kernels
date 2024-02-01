from generate_graphs import (interpolate_ER_density,
                             interpolate_ER_PPM,
                             interpolate_ER_GRG_torus,
                             interpolate_ER_inhomogeneous,
                             interpolate_GRG_torus_circle,
                             interpolate_ER_GCG)
from experiment import Experiment,selected_kernels,tuple2str
from scipy.stats import spearmanr
from collections import defaultdict
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
from other_kernels import NetLSD,Gin,GraphletSampling4,DegreeHistogram
import numpy as np

fast_kernels = (
     PyramidMatch,
     NeighborhoodHash,
     ShortestPath,
     WeisfeilerLehman,
     Propagation,
     OddSth,
     WeisfeilerLehmanOptimalAssignment,
    )
just_netlsd_gin = (Gin,NetLSD)
just_graphlet4 = [GraphletSampling4]
interpolators = [interpolate_ER_density,interpolate_ER_size,interpolate_ER_GCG, interpolate_ER_inhomogeneous, interpolate_ER_triangular,interpolate_ER_GRG_torus, interpolate_GRG_torus_circle,interpolate_ER_PPM]

nsteps = 11
steps = np.linspace(0,1,nsteps)
parameters = [{'step': s} for s in steps]
nsamples = 100
npacks = 30

def calc_step_mmd_spearman(mmds):
    steps_dict = defaultdict(list)
    mmds_dict = defaultdict(list)
    for (locator1,locator2),vals in mmds.items():
        for v in vals:
            steps_dict[locator1[0]].append(float(locator1[-1]))
            mmds_dict[locator1[0]].append(v)
    return {
        m: spearmanr(steps_dict[m],mmds_dict[m]).statistic
        for m in steps_dict.keys()
    }

def load_mmds(fn):
    import json
    with open(fn) as f:
        mmds = json.load(f)
    str2locator = lambda s: (s.split(', ')[0],float(s.split(', ')[1]))
    return {
        interpolator: {
            tuple(map(str2locator,key.split('_vs_'))): vals
            for key,vals in ivals.items()
        }
        for interpolator,ivals in mmds.items()
    }

def scatter_mmds(mmds,transition_name,ax=None):
    steps_dict = defaultdict(list)
    mmds_dict = defaultdict(list)
    step_star = defaultdict(float)
    for (locator1,locator2),vals in mmds.items():
        step_star[locator1[0]] = int(locator2[-1])
        for v in vals:
            steps_dict[locator1[0]].append(abs(float(locator1[-1])-float(locator2[-1])))
            mmds_dict[locator1[0]].append(v)
    
    if ax is None:
        import matplotlib.pyplot as plt
        _,ax = plt.subplots()
    for m in steps_dict.keys():
        ax.set_xlabel(r'$\theta$' if step_star[m]==0 else r'$1-\theta$')
        ax.set_ylabel(r'MMD')
        ax.scatter(steps_dict[m],mmds_dict[m])
        ax.set_title(f'$r_{step_star[m]}={spearmanr(steps_dict[m],mmds_dict[m]).statistic:.03f}$')




load = False
perform_start = True
perform_end = True
for interpolator in interpolators:
    g_name = interpolator.__name__
    print('start',g_name)
    g_file = f'{g_name}_graphs.json'
    start_vals_file = f'{g_name}_start_vals.tsv'
    start_mmds_file = f'{g_name}_start_mmds.json'
    start_spearmans_file = f'{g_name}_start_spearmans.tsv'
    end_vals_file = f'{g_name}_end_vals.tsv'
    end_mmds_file = f'{g_name}_end_mmds.json'
    end_spearmans_file = f'{g_name}_end_spearmans.tsv'
    experiment = Experiment([interpolator],parameters,npacks,sample_size=nsamples)
    
    if load:
        experiment.load_graphs(g_file)
        print('loaded graphs')
    else:
        npacks_dict = {}
        for interpolator in interpolators:
            npacks_dict[interpolator.__name__] = {}
            npacks_dict[interpolator.__name__][tuple2str(parameters[0].values())] = npacks*(nsteps+1)
            npacks_dict[interpolator.__name__][tuple2str(parameters[-1].values())] = npacks*(nsteps+1)
        experiment.generate_graphs(g_file,npacks_dict=npacks_dict)
        print('Generated graphs',flush=True)
    
    if perform_start:
        start_mmds=experiment.apply_kernels(experiment.iterator_transitions_startcomparison(),save_name=start_vals_file,save_mmds_name=start_mmds_file,kernels=selected_kernels)
        #start_mmds=experiment.load_mmds('interpolation_start_mmds.json')
        with open(start_spearmans_file,'w') as save_file:
            for k,kmmds in start_mmds.items():
                for m,s in calc_step_mmd_spearman(kmmds).items():
                    print("\t".join(['#',k,m,str(s)]),flush=True,file=save_file)

    if perform_end:
        end_mmds=experiment.apply_kernels(experiment.iterator_transitions_endcomparison(),save_name=end_vals_file,save_mmds_name=end_mmds_file,kernels=selected_kernels)
        with open(end_spearmans_file,'w') as save_file:
            for k,kmmds in end_mmds.items():
                for m,s in calc_step_mmd_spearman(kmmds).items():
                    print("\t".join(['#',k,m,str(s)]),flush=True,file=save_file)

