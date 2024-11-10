# README

To create a virtual environment with all required packages for the experiments, run

````bash
conda create --name graphkernels python=3.7 grakel pandas matplotlib pip plotly --channel conda-forge
conda activate graphkernels
pip install python-igraph netlsd kaleido
git clone https://github.com/uoguelph-mlrg/GGM-metrics GgmMetrics
cd GgmMetrics/
pip install -r requirements.txt
pip install dgl==0.6.1
cd ../
ln -s GgmMetrics/evaluation


````



To generate graphs and perform the experiment, simply run the file `spearman_experiment.py` with python:

```bash
python spearman_experiment.py
```

To generate a latex table, run `latex_tables.py`, and to generate the figures, run `plot_all.py`.