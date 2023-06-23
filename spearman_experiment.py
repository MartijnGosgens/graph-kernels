from generate_graphs import interpolate_ER_PPM,interpolate_ER_GRG_torus,interpolate_ER_inhomogeneous,interpolate_GRG_torus_circle
from experiment import Experiment,selected_kernels
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

interpolators = [interpolate_ER_PPM,interpolate_ER_GRG_torus,interpolate_ER_inhomogeneous, interpolate_GRG_torus_circle]
nsteps = 11
steps = np.linspace(0,1,nsteps)
parameters = [{'step': s} for s in steps]
nsamples = 30
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

print('start experiment')
experiment = Experiment(interpolators,parameters,npacks,sample_size=nsamples)
experiment.generate_graphs('interpolation_graphs2.json')
print('Generated graphs')
#experiment.load_graphs('interpolation_graphs2.json')
#print('loaded graphs')
start_mmds=experiment.apply_kernels(experiment.iterator_transitions_startcomparison(),save_name='interpolation_start_vals.tsv',save_mmds_name='interpolation_start_mmds.json')
#start_mmds=experiment.load_mmds('interpolation_start_mmds.json')
with open('spearmans_start.tsv','w') as save_file:
    for k,kmmds in start_mmds.items():
        for m,s in calc_step_mmd_spearman(kmmds).items():
            print("\t".join(['#',k,m,str(s)]),flush=True,file=save_file)

print('end experiment')
end_mmds=experiment.apply_kernels(experiment.iterator_transitions_endcomparison(),save_name='interpolation_end_vals.tsv',save_mmds_name='interpolation_end_mmds.json')
with open('spearmans_end.tsv','w') as save_file:
    for k,kmmds in end_mmds.items():
        for m,s in calc_step_mmd_spearman(kmmds).items():
            print("\t".join(['#',k,m,str(s)]),flush=True,file=save_file)
