from typing import Dict, List, Callable, Any
from generate_graphs import *

import plotly.graph_objects as go
import networkx as nx
import numpy as np
import netlsd


NX_Graph = nx.classes.graph.Graph # define networkx Graph type
Parameters = Dict[str, Any] # define parameters dict

# define descriptor function
DESCRIPTOR: Callable[[NX_Graph], np.ndarray] = netlsd.wave

DESCRIPTOR_OPTIONS: Dict[str, Any] = {
        "timescales": np.linspace(1, 10, 100),#np.logspace(-2, 2, 250),
    } # some descriptor optins e.g. dimensionality of a vector for netLSD

# number of graphs, number of nodes:
NUM_OF_GRAPHS_FOR_EACH_MODEL = 100

# define aliases for GNP and PA: (typo-safe)
ER = "ER"
CL = "CL"
PP = "PP"
torus = "torus"
circle = "circle"
SC = "SC"

model_name_to_generator: Dict[str, Callable[..., NX_Graph]] = {
    ER: lambda: ig2nx(interpolate_ER_inhomogeneous(step=0.0)),
    CL: lambda: ig2nx(interpolate_ER_inhomogeneous(step=1.0)),
    PP: lambda: ig2nx(interpolate_ER_PPM(step=1.0)),
    torus: lambda: ig2nx(interpolate_ER_GRG_torus(step=1.0)),
    circle: lambda: ig2nx(interpolate_GRG_torus_circle(step=1.0)),
    SC: lambda: ig2nx(interpolate_ER_GCG(step=1.0)),
} # finite machine for graphs generation

interpolators = {
    "Density": interpolate_ER_density,
    "Heterogeneity": interpolate_ER_inhomogeneous,
    "Communities": interpolate_ER_PPM,
    "Geometry": interpolate_ER_GRG_torus,
    "Dimensionality": interpolate_GRG_torus_circle,
    "Complementarity": interpolate_ER_GCG
}

def calculate_descriptors_mean(graphs: List[NX_Graph],
                               desciptor_func: Callable[[NX_Graph], np.ndarray],
                               descriptor_params: Dict[str, Any]={},
                               ):
    
    descriptors = np.array([desciptor_func(g, **descriptor_params) for g in graphs])
    return descriptors.mean(axis=0)

interpolator_name="Dimensionality"
interpolator=interpolate_GRG_torus_circle
fig = go.Figure()

for step in [0.0,0.6,1.0]:
    label = f"$\\theta = {step:.01f}$"
    graphs: List[NX_Graph] = [
        ig2nx(interpolator(step=step)) for _ in range(NUM_OF_GRAPHS_FOR_EACH_MODEL)
    ]
        
    averaged_descriptor = calculate_descriptors_mean(graphs=graphs,
                                                    desciptor_func=DESCRIPTOR,
                                                    descriptor_params=DESCRIPTOR_OPTIONS,
                                                    )
            
    
    fig.add_trace(
        go.Scatter(
            x=DESCRIPTOR_OPTIONS["timescales"], 
            y=averaged_descriptor, 
            name=label,
        ),
        
    )
fig.show()
fig.write_image('test.pdf')
# Do it again to prevent the "Loading [MathJax]/extensions/MathMenu.js" popup from appearing
fig = go.Figure()

for step in [0.0,0.6,1.0]:
    label = f"$\\theta = {step:.01f}$"
    graphs: List[NX_Graph] = [
        ig2nx(interpolator(step=step)) for _ in range(NUM_OF_GRAPHS_FOR_EACH_MODEL)
    ]
        
    averaged_descriptor = calculate_descriptors_mean(graphs=graphs,
                                                    desciptor_func=DESCRIPTOR,
                                                    descriptor_params=DESCRIPTOR_OPTIONS,
                                                    )
            
    
    fig.add_trace(
        go.Scatter(
            x=DESCRIPTOR_OPTIONS["timescales"], 
            y=averaged_descriptor, 
            name=label,
        ),
        
    )
fig.show()
        
#fig.update_xaxes(type="log")
fig.update_layout(xaxis_title="Timescales", width=500, height=300, margin=dict(l=20, r=20, t=20, b=20, pad=0))

fig.write_image("NetLSDWave_dimensionality_just3.pdf".format(interpolator_name))


fig = go.Figure()

for model_name,model_generator in model_name_to_generator.items():
        
    # define label from model name and unpacked parameters values:
    label = model_name

    graphs: List[NX_Graph] = [
        model_generator() for _ in range(NUM_OF_GRAPHS_FOR_EACH_MODEL)
    ]
    
    averaged_descriptor = calculate_descriptors_mean(graphs=graphs,
                                                        desciptor_func=DESCRIPTOR,
                                                        descriptor_params=DESCRIPTOR_OPTIONS,
                                                        )
            
    
    fig.add_trace(
        go.Scatter(
            x=DESCRIPTOR_OPTIONS["timescales"], 
            y=averaged_descriptor, 
            name=label,
        ),
        
    )
        
        
#fig.update_xaxes(type="log")
fig.update_layout(title="NetLSD Wave traces of interpolation endpoints", xaxis_title="Timescales")
fig.write_image("NetLSDWave_endpoints.svg")

for interpolator_name,interpolator in interpolators.items():
    fig = go.Figure()

    for step in np.linspace(0,1,11):
        label = f"$\\theta = {step:.01f}$"
        model_generator = model_name_to_generator[model_name]
        graphs: List[NX_Graph] = [
            ig2nx(interpolator(step=step)) for _ in range(NUM_OF_GRAPHS_FOR_EACH_MODEL)
        ]
            
        averaged_descriptor = calculate_descriptors_mean(graphs=graphs,
                                                        desciptor_func=DESCRIPTOR,
                                                        descriptor_params=DESCRIPTOR_OPTIONS,
                                                        )
                
        
        fig.add_trace(
            go.Scatter(
                x=DESCRIPTOR_OPTIONS["timescales"], 
                y=averaged_descriptor, 
                name=label,
            ),
            
        )
            
            
    #fig.update_xaxes(type="log")
    fig.update_layout(title="NetLSD Wave traces of {}".format(interpolator_name), xaxis_title="Timescales")
    fig.write_image("NetLSDWave_{}.svg".format(interpolator_name))
    
