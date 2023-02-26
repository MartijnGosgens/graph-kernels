import json
import numpy as np
from collections import defaultdict
from scipy.stats import wilcoxon, mannwhitneyu,ttest_ind
import matplotlib.pyplot as plt


def stat_tests(sample1,sample2):
    try:
        return mannwhitneyu(sample1,sample2)[1]
    except:
        return None

st2 = json.loads(open('mmd_matrix.json').read())
for k in st2:
    knl, g1, g2 = k.split('\t')
    k2 = f"{knl}\t{g1}\t{g1}"
    t1 = f'{g1}_vs_{g2}'
    t2 = f'{g1}_vs_{g1}'
    try:
        pv = stat_tests(st2[k],st2[k2])
    except:
        pv = 1.
    _k = k.replace('\t',' ')
    _kk = k.split('\t')
    print(f'{knl}\t{t1}\t{np.mean(st2[k])}\t{t2}\t{np.mean(st2[k2])}\tpval\t{pv}')

    try:
        plt.clf()
        plt.hist(st2[k],label=t1, alpha=.5)
        plt.hist(st2[k2],label=t2, alpha=.5)

        suf = ''
        if pv<.1:
            suf = '_diff'
        plt.title(f'{t1} {np.mean(st2[k]):.04f} {t2} {np.mean(st2[k2]):.04f} pval {pv:.02%}')
        plt.savefig(f'mmds_png/mmds{suf}_{_kk[0]}_{_kk[1]}_vs_{_kk[2]}.png')
    except:
        pass