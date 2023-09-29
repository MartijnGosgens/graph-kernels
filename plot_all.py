from spearman_experiment import *
import matplotlib.pyplot as plt

from generate_graphs import (interpolate_ER_triangular,
                             interpolate_ER_PPM,
                             interpolate_ER_GRG_torus,
                             interpolate_ER_inhomogeneous,
                             interpolate_GRG_torus_circle,
                             interpolate_ER_GCG, edges2nx)
import json
import networkx as nx
import numpy as np
from grakel.kernels import GraphletSampling
from collections import defaultdict
import pandas as pd

scatter_mmds(load_mmds('interpolate_ER_inhomogeneous_start_mmds.json')['ShortestPath'],transition_name=r'ER$\leftrightarrow$CL')
plt.savefig('scatter_ER_inhomogeneous_start.jpg',bbox_inches='tight')

scatter_mmds(load_mmds('interpolate_ER_inhomogeneous_end_mmds.json')['ShortestPath'],transition_name=r'ER$\leftrightarrow$CL')
plt.savefig('scatter_ER_inhomogeneous_end.jpg',bbox_inches='tight')


#### DEGREE PLOT
interpolators = [interpolate_ER_triangular, interpolate_ER_PPM, interpolate_ER_GRG_torus, interpolate_ER_inhomogeneous, interpolate_GRG_torus_circle, interpolate_ER_GCG]


transition2name = {
    interpolate_ER_inhomogeneous: r'ER$\leftrightarrow$CL',
    interpolate_ER_triangular: r'ER$\leftrightarrow$Triadic',
    interpolate_ER_PPM: r'ER$\leftrightarrow$PP',
    interpolate_ER_GRG_torus: r'ER$\leftrightarrow$Torus',
    interpolate_GRG_torus_circle: r'Torus$\leftrightarrow$Circle',
    interpolate_ER_GCG: r'ER$\leftrightarrow$SC',
}
packvar = lambda pack: np.array(sum([
        [d for _,d in G.degree()]
        for G in pack
    ],[])).var()
packtransitivity = lambda pack: sum(map(nx.transitivity,pack))/len(pack)

fig,ax = plt.subplots()
for interpolator in interpolators:
    file_name = f'{interpolator.__name__}_graphs.json'
    with open(file_name) as f:
        graphs = json.load(f)
    graphs = graphs[interpolator.__name__]
    #print(list(graphs.keys()))
    firstpacks = {
        float(step): list(map(edges2nx,packs[0]))
        for step,packs in graphs.items()
    }
    ax.plot(list(firstpacks.keys()),[
        packvar(pack)
        for pack in firstpacks.values() 
    ], label=transition2name[interpolator])
    ax.set_xlabel(r'$\theta$')
    ax.set_ylabel('Degree variance')
plt.legend()
plt.savefig('Degree variance.jpg',bbox_inches='tight')

#### GRAPHLET COUNTS


load=True
graphlet_size = 3
normalize = False

interpolators = [interpolate_ER_triangular, interpolate_ER_PPM, interpolate_ER_GRG_torus, interpolate_ER_inhomogeneous, interpolate_GRG_torus_circle, interpolate_ER_GCG]

transition2name = {
    interpolate_ER_inhomogeneous: r'ER$\leftrightarrow$CL',
    interpolate_ER_triangular: r'ER$\leftrightarrow$Triadic',
    interpolate_ER_PPM: r'ER$\leftrightarrow$PP',
    interpolate_ER_GRG_torus: r'ER$\leftrightarrow$Torus',
    interpolate_GRG_torus_circle: r'Torus$\leftrightarrow$Circle',
    interpolate_ER_GCG: r'ER$\leftrightarrow$SC',
}

model2location = {
    'ER': (interpolate_ER_inhomogeneous,'0.0'),
    'CL': (interpolate_ER_inhomogeneous,'1.0'),
    'Triadic': (interpolate_ER_triangular,'1.0'),
    'PP': (interpolate_ER_PPM,'1.0'),
    'Torus': (interpolate_ER_GRG_torus,'1.0'),
    'Circle': (interpolate_GRG_torus_circle,'1.0'),
    'SC': (interpolate_ER_GCG,'1.0'),
}

def compute_graphlet_counts(model2location=model2location,save_name='pack1_graphlets.json'):
    graphlet2model2counts = defaultdict(lambda: defaultdict(list))
    selected_graphs = []
    generator = []
    for model,(interpolator,step) in model2location.items():
        print(f'Loading {interpolator.__name__} {step}')
        file_name = f'{interpolator.__name__}_graphs.json'
        with open(file_name) as f:
            graphs = json.load(f)
        graphs = graphs[interpolator.__name__]
        #print(list(graphs.keys()))
        selected_graphs += list(map(edges2grakel,graphs[step][0]))
        generator += [model]*len(graphs[step][0]) 
    print('Counting graphlets')
    k = GraphletSampling(normalize=True,k=graphlet_size)
    k.fit_transform(selected_graphs)
    bin_descs = [
        str(bin)
        for bin in k._graph_bins.values()
    ]
    #bin_descs = ['Wedges','Triangles']
    for v,gen in zip(k._phi_X,generator):
        if normalize:
            v = np.array(v)/(sum(np.array(v)**2)**0.5)
        for i,x in enumerate(v):
            graphlet2model2counts[bin_descs[i]][gen].append(x)

    with open(save_name, 'w') as save_file:
        json.dump(graphlet2model2counts, save_file)
    return graphlet2model2counts

rename_dict = {
    "<(0 0 [1, 2])(1 0 [0])(2 0 [0])>": "Wedges",
    "<(0 0 [1, 2])(1 0 [0, 2])(2 0 [0, 1])>": "Triangles",
}

def load_graphlet_counts(file_name=f'pack1_graphlets_k{graphlet_size}.json'):
    with open(file_name) as f:
        graphlet2model2counts = json.load(f)
    return graphlet2model2counts


def rename_graphlet_counts(graphlet2model2counts,rename_dict=rename_dict):
    return {
        rename_dict[graphlet] if graphlet in rename_dict else graphlet: model2counts
        for graphlet,model2counts in graphlet2model2counts.items()
    }

graphlet2model2counts = {}
if load:
    graphlet2model2counts = load_graphlet_counts(file_name=f'pack1_graphlets_k{graphlet_size}_{"un" if not normalize else ""}normalized.json')
else:
    graphlet2model2counts = compute_graphlet_counts(save_name=f'pack1_graphlets_k{graphlet_size}_{"un" if not normalize else ""}normalized.json')
graphlet2model2counts = rename_graphlet_counts(graphlet2model2counts)

data = [
    [graphlet]+[
        np.array(counts).mean()
        for counts in model2counts.values()
    ]
    for graphlet,model2counts in graphlet2model2counts.items()
]
print(data)
df = pd.DataFrame(data, columns = ['Graphlet']+list(model2location.keys()))
ax=df.plot(x='Graphlet', kind='bar', stacked=False, title=f'Graphlet counts for $k={graphlet_size}$',rot=0)
if normalize:
    ax.set_ylabel('Normalized graphlet count')
else:
    ax.set_ylabel('Graphlet count')

plt.savefig(f'graphlet_histograph_k{graphlet_size}_{"un" if not normalize else ""}normalized.jpg',bbox_inches='tight')
    
