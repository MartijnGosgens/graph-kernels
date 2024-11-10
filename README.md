# README

The experiments require matplotlib, [netlsd](https://pypi.org/project/NetLSD/), [grakel](https://github.com/ysig/GraKeL), [networkx](https://networkx.org/) and [igraph](https://python.igraph.org/en/stable/).

To create a virtual environment with these packages, run

````bash
conda create --name graphkernels python=3.12 grakel pandas networkx matplotlib pip plotly conda dglteam::dgl --channel conda-forge
pip install python-igraph netlsd kaleido torch


````



To generate graphs and perform the experiment, simply run the file `spearman_experiment.py` with python:

```bash
python spearman_experiment.py
```

To generate a latex table, run `latex_tables.py`, and to generate the figures, run `plot_all.py`.