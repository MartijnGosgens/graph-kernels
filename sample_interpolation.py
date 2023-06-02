import json
import sys
import numpy as np
from generate_graphs import *
''' Example usage
    python sample_interpolation.py 11 50 interpolation_graphs.json
'''

interpolators = [interpolate_ER_PPM,interpolate_ER_GRG_torus,interpolate_ER_inhomogeneous]
nsteps = int(sys.argv[1]) if len(sys.argv)>1 else 11
steps = np.linspace(0,1,11)
nsamples = int(sys.argv[2]) if len(sys.argv)>2 else 50
savename = sys.argv[3] if len(sys.argv)>3 else 'interpolation_graphs.json'
res = dict()
res['interpolators'] = dict()
res['stats'] = dict()

for interpolator in interpolators:
    name = interpolator.__name__
    print('Sampling',name)
    res['interpolators'][name] = dict()
    res['stats'][name] = dict()
    for step in steps:
        print('Sampling step',step)
        res['interpolators'][name][step] = []
        res['stats'][name][step] = []
        for _ in range(nsamples):
            g = interpolator(step)
            es = ig2edges(g)
            res['interpolators'][name][step].append( es[:] )
            res['stats'][name][step].append( nx2stats(edges2nx(es)) )
    print('done')
print(json.dumps(res), file=open(savename, 'w'))

