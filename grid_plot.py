from generate_graphs import interpolate_GRG_torus_circle
from other_kernels import NetLSDWave,NetLSD
from experiment import Experiment
import numpy as np
nsteps = 11
steps = np.linspace(0,1,nsteps)
parameters = [{'step': s} for s in steps]
nsamples = 100
npacks = 10
kernels = [NetLSD,NetLSDWave]


interpolator = interpolate_GRG_torus_circle
g_name = interpolator.__name__
g_file = f'smallgraphs/{g_name}_graphs.json'
vals_file = f'{g_name}_grid_vals.tsv'
mmds_file = f'{g_name}_grid_mmds.json'

import json
with open(f'{g_name}_grid_mmds.json') as f:
    all_mmds = json.load(f)


import matplotlib.pyplot as plt
import numpy as np

import matplotlib
import matplotlib as mpl

steps = [0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]

for k in kernels:
    k_name = k.__name__
    mmds = {
        tuple([float(v.split(' ')[-1][:3]) for v in k.split('_vs_')]): vs
        for k,vs in all_mmds[k_name].items()
    }
    data = np.array(
        [
            [
                abs(np.array(list(map(float,mmds[s1,s2]))).mean()) for s2 in steps
            ]
            for s1 in steps
        ]
    )
    fig, ax = plt.subplots()
    im = ax.imshow(data)

    # Create colorbar
    cbar = ax.figure.colorbar(im, ax=ax, cmap="reds")
    cbar.ax.set_ylabel("MMD", rotation=90, va="bottom")

    # Show all ticks and label them with the respective list entries
    ax.set_xticks(np.arange(len(steps)), labels=map(str,steps))
    ax.set_yticks(np.arange(len(steps)), labels=map(str,steps))
    ax.set_xlabel('$\\theta_1$')
    ax.set_ylabel('$\\theta_2$')

    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
            rotation_mode="anchor")


    #ax.set_title("NetLSD Wave MMD values on dimensionality interpolation")
    fig.tight_layout()
    plt.savefig(f'dimensionality_{k_name}_heatmap.svg')
    plt.show()