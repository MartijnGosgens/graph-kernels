import json
import numpy as np
from scipy.stats import mannwhitneyu,ttest_ind
import matplotlib.pyplot as plt
import networkx as nx
from collections import defaultdict

fn = '1k_samples.json'

def stat_tests(sample1,sample2):
    try:
        return mannwhitneyu(sample1,sample2)[1]
    except:
        return None


d = json.loads(open(fn).read())
gens = list(d['stats'])

pools = defaultdict(list)
for g in gens:
    # print(g)
    for gg in d['graphs'][g]:
        it = nx.from_edgelist(gg)
        for _,dsts in nx.all_pairs_shortest_path_length(it):
            pools[g].extend(list(dsts.values()))

par = 'all_shortestpath_len'
# print(f'\n\n{par}')
for i1, g1 in enumerate(gens):
    for i2, g2 in enumerate(gens):
        if i1<=i2:
            print(f'{par}\t{g1}\t{np.mean(pools[g1])}\t{g2}\t{np.mean(pools[g2])}\tpval\t{stat_tests(pools[g1],pools[g2])}')

plt.clf()
plt.title(par)
for g in pools:
    plt.hist(pools[g], label=g, alpha=.5)
plt.legend()
plt.savefig(f'base_stats_{par}.png')
