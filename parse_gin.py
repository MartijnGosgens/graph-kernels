import json
import numpy as np
from collections import defaultdict
from scipy.stats import wilcoxon
import matplotlib.pyplot as plt

_INF = 1e10
fn = 'gin_raw.tsv'
SAMPLES = 30
PACKS = 30
# GENS = ['ER', 'PPM', 'GRG', 'GRGs', 'GRGh', 'GRGt', 'PA']

# st = defaultdict(lambda:defaultdict(list))
st = dict()

for idx, line in enumerate(open(fn)):
    if '\t' not in line: continue
    chunks = line.strip().split('\t')
    key = str(chunks[1])+"\t"+str(chunks[2])
    if key not in st:
        st[key] = dict()
    skey = str(chunks[3])+"\t"+str(chunks[5])
    if skey not in st[key]:
        st[key][skey] = []
    try:
        val = float(chunks[-1])
    except:
        continue
    if abs(val)>_INF:
        continue
    st[key][skey].append([chunks[4], chunks[6], val])
    if not idx%1000000:
        print(idx)

with open('gin.vals.json', 'w') as ofh:
    print(json.dumps(st), file=ofh)

