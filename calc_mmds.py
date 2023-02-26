import json
import numpy as np
from collections import defaultdict
from scipy.stats import wilcoxon
import matplotlib.pyplot as plt

_INF = 1e10
SAMPLES = 30
PACKS = 30
GENS = ['ER', 'PPM', 'GRG', 'GRGs', 'GRGh', 'GRGt', 'PA']

fn = 'kernels_v2.vals.json'
st = json.loads(open(fn).read())

print('loaded')
st2 = dict()

for _as in st:
    a,s = _as.split('\t')
    print(_as)
    for i1, m1 in enumerate(GENS):
        for i2, m2 in enumerate(GENS):
            kk = f'{a}\t{m1}\t{m2}'
            if kk not in st2:
                st2[kk] = list()

            if i1!=i2:
                vl = [v for q,p,v in st[_as][f'{m1}\t{m1}'] if int(q)<SAMPLES and int(p)<SAMPLES and not np.isnan(v) and not np.isinf(v) and abs(v)<_INF]
                val = np.mean(vl)
                vl = [v for q,p,v in st[_as][f'{m2}\t{m2}'] if int(q)<SAMPLES and int(p)<SAMPLES and not np.isnan(v) and not np.isinf(v) and abs(v)<_INF]
                val += np.mean(vl)
                vl = [v for q,p,v in st[_as][f'{m1}\t{m2}'] if int(q)<SAMPLES and int(p)<SAMPLES and not np.isnan(v) and not np.isinf(v) and abs(v)<_INF]
                val -= 2*np.mean(vl)
                st2[kk].append(val)
            else:
                vl = [v for q,p,v in st[_as][f'{m1}\t{m1}'] if int(q)<SAMPLES and int(p)<SAMPLES and not np.isnan(v) and not np.isinf(v) and abs(v)<_INF]
                val = np.mean(vl)
                vl = [v for q,p,v in st[_as][f'{m2}\t{m2}'] if SAMPLES<=int(q)<2*SAMPLES and SAMPLES<=int(p)<2*SAMPLES and not np.isnan(v) and not np.isinf(v) and abs(v)<_INF]
                val += np.mean(vl)
                vl = [v for q,p,v in st[_as][f'{m1}\t{m2}'] if int(q)<SAMPLES and SAMPLES<=int(p)<2*SAMPLES and not np.isnan(v) and not np.isinf(v) and abs(v)<_INF]
                val -= 2*np.mean(vl)
                st2[kk].append(val)

print('saving')
with open('mmd_matrix.json', 'w') as ofh:
    print(json.dumps(st2), file=ofh)
print('done')
