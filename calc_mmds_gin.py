import json
import numpy as np
from collections import defaultdict
from scipy.stats import wilcoxon
import matplotlib.pyplot as plt

_INF = 1e10
SAMPLES = 30
PACKS = 30
GENS = ['ER', 'PPM', 'GRG', 'GRGs', 'GRGh', 'GRGt', 'PA']


fn = 'gin.vals.json'

st = json.loads(open(fn).read())

print('loaded')
st2 = dict() # defaultdict(list)

for _as in st:
    a,s = _as.split('\t')
    print(_as)
    for i1, m1 in enumerate(GENS):
        for i2, m2 in enumerate(GENS):
            kk = f'{a}\t{m1}\t{m2}'
            if kk not in st2:
                st2[kk] = list()
            # if i1>i2:
            #     continue
            if i1!=i2:
                vl = [v for q,p,v in st[_as][f'{m1}\t{m1}'] if int(q)<SAMPLES and int(p)<SAMPLES and not np.isnan(v) and not np.isinf(v) and abs(v)<_INF]
                # print(vl)
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

                # val = np.mean(st[(a,s)][(m1,m1)])
                # val += np.mean(st[(a,s)][(m2,m2)])
                # val -= 2*np.mean(st[(a,s)][(m1,m2)])
                # st2[(a,m1,m2)].append(val)
#    exit()
print('saving')
with open('gin_mmd_matrix.json', 'w') as ofh:
    print(json.dumps(st2), file=ofh)
print('done')
# for k in st2:
#     try:
#         wx = wilcoxon(st2[k]).pvalue
#     except:
#         wx = 1.
#     print(k,'\tmean:\t', np.mean(st2[k]),'\tstd:\t', np.std(st2[k]), '\twilcoxon:\t',wx)

#     plt.clf()
#     # for g in pools:
#     plt.hist(st2[k])# , label=g, alpha=.5)
#     # plt.legend()
#     suf = ''
#     if wx>.1:
#         suf = '_diff'
#     plt.title(f'{k} mean {np.mean(st2[k]):.04f} std {np.std(st2[k]):.04f} pval {wx:.02%}')
#     plt.savefig(f'mmds{suf}_{k[0]}_{k[1]}_vs_{k[2]}.png')
