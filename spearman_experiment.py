from generate_graphs import (interpolate_ER_triangular,
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
from other_kernels import NetLSD,Gin
import numpy as np

fast_kernels = (
     #RandomWalk, # ERRRs,
     #GraphletSampling,
     PyramidMatch,
     NeighborhoodHash,
     ShortestPath,
     WeisfeilerLehman,
     Propagation,
     OddSth,
     WeisfeilerLehmanOptimalAssignment,
     #NeighborhoodSubgraphPairwiseDistance,
     # SvmTheta, # ERRR
    )
just_netlsd_gin = (Gin,NetLSD)

interpolators = [interpolate_ER_GCG,interpolate_ER_PPM,interpolate_ER_GRG_torus,interpolate_ER_inhomogeneous, interpolate_GRG_torus_circle,interpolate_ER_triangular]

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

load = False
for interpolator in interpolators:
    g_name = interpolator.__name__
    print('start',g_name)
    g_file = f'{g_name}_graphs.json'
    start_vals_file = f'{g_name}_start_vals.tsv'
    start_mmds_file = f'{g_name}_start_mmds.json'
    start_spearmans_file = f'{g_name}_end_spearmans.tsv'
    end_vals_file = f'{g_name}_end_vals.tsv'
    end_mmds_file = f'{g_name}_end_mmds.json'
    end_spearmans_file = f'{g_name}_end_spearmans.tsv'
    experiment = Experiment([interpolator],parameters,npacks,sample_size=nsamples)
    npacks_dict = npacks_dict={
        interpolator.__name__: {
            tuple2str(parameters[0].values()): npacks*(nsteps+1),
            tuple2str(parameters[-1].values()): npacks*(nsteps+1),
        }
        for interpolator in interpolators 
    }
    if load:
        experiment.load_graphs(g_file)
        print('loaded graphs')
    else:
        experiment.generate_graphs(g_file,npacks_dict=npacks_dict)
        print('Generated graphs')
    
    start_mmds=experiment.apply_kernels(experiment.iterator_transitions_startcomparison(),kernels=selected_kernels,save_name=start_vals_file,save_mmds_name=start_mmds_file)
    #start_mmds=experiment.load_mmds('interpolation_start_mmds.json')
    with open(start_spearmans_file,'w') as save_file:
        for k,kmmds in start_mmds.items():
            for m,s in calc_step_mmd_spearman(kmmds).items():
                print("\t".join(['#',k,m,str(s)]),flush=True,file=save_file)

    print('end',g_name)
    end_mmds=experiment.apply_kernels(experiment.iterator_transitions_endcomparison(),kernels=selected_kernels,save_name=end_vals_file,save_mmds_name=end_mmds_file)
    with open(end_spearmans_file,'w') as save_file:
        for k,kmmds in end_mmds.items():
            for m,s in calc_step_mmd_spearman(kmmds).items():
                print("\t".join(['#',k,m,str(s)]),flush=True,file=save_file)

'''print('start experiment')
experiment = Experiment(interpolators_just_one,parameters,npacks,sample_size=nsamples)
npacks_dict = npacks_dict={
    interpolator.__name__: {
        tuple2str(parameters[0].values()): npacks*(nsteps+1),
        tuple2str(parameters[-1].values()): npacks*(nsteps+1),
    }
    for interpolator in interpolators 
}
print(npacks_dict)
experiment.generate_graphs('interpolation_graphs_just_one.json',npacks_dict=npacks_dict)
print('Generated graphs')
#experiment.load_graphs('interpolation_graphs_just_one.json')
#print('loaded graphs')
start_mmds=experiment.apply_kernels(experiment.iterator_transitions_startcomparison(),kernels=selected_kernels,save_name='interpolation_start_vals.tsv',save_mmds_name='interpolation_start_mmds.json')
#start_mmds=experiment.load_mmds('interpolation_start_mmds.json')
with open('spearmans_start.tsv','w') as save_file:
    for k,kmmds in start_mmds.items():
        for m,s in calc_step_mmd_spearman(kmmds).items():
            print("\t".join(['#',k,m,str(s)]),flush=True,file=save_file)

print('end experiment')
end_mmds=experiment.apply_kernels(experiment.iterator_transitions_endcomparison(),kernels=selected_kernels,save_name='interpolation_end_vals.tsv',save_mmds_name='interpolation_end_mmds.json')
with open('spearmans_end.tsv','w') as save_file:
    for k,kmmds in end_mmds.items():
        for m,s in calc_step_mmd_spearman(kmmds).items():
            print("\t".join(['#',k,m,str(s)]),flush=True,file=save_file)
'''

