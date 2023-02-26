import json
import sys
from generate_graphs import *


res = dict()
res['graphs'] = dict()
res['stats'] = dict()

for ggn, gg in zip(['ER', 'PPM', 'GRG', 'GRGs', 'GRGh', 'GRGt', 'PA'], [generate_ER, generate_PPM, generate_GRG, generate_GRG_square, generate_GRG_hypersphere, generate_GRG_torus, generate_PA]):
    print(ggn)
    res['graphs'][ggn] = []
    res['stats'][ggn] = []
    for _ in range(int(sys.argv[1])):
        g = gg()
        es = ig2edges(g)
        res['graphs'][ggn].append( es[:] )
        res['stats'][ggn].append( nx2stats(edges2nx(es)) )
    print('done')
print(json.dumps(res), file=open(sys.argv[2], 'w'))

