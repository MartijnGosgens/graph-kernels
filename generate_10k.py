import json
from generate_graphs import *


res = dict()
res['graphs'] = dict()
res['stats'] = dict()

for ggn, gg in zip(['ER', 'PPM', 'GRGt', 'GRGs'], [generate_ER, generate_PPM, generate_GRG, generate_GRG_square]):
    res['graphs'][ggn] = []
    res['stats'][ggn] = []
    for _ in range(10000):
        g = gg()
        es = ig2edges(g)
        res['graphs'][ggn].append( es[:] )
        res['stats'][ggn].append( nx2stats(edges2nx(es)) )
        
print(json.dumps(res), file=open('10k_samples.json', 'w'))

