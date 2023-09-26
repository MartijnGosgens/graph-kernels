from glob import glob
import json
from matplotlib import pyplot as plt
import numpy as np

def draw_distr(vals, title, ofn):
    plt.clf()
    n1,v1 = zip(*sorted(vals.items()))
    print(title)
    print(n1)
    print(v1)
    # width = .3
    plt.bar(np.array(n1), np.array(v1)/sum(v1))
    plt.suptitle(title)
    # plt.legend()
    plt.tight_layout()
    print(ofn)
    plt.savefig(ofn)
    # plt.show() 
    
# _ = draw_distr(degrees[0], 'Degrees', '_test.png')


seen = set()
def keystoint(x, remove=None):
    r = {int(k): v for k, v in x.items()}
    if remove is not None:
        r[remove] = 0
    return r

for ff in glob("graphstats/*_undir_degrees.json"):
    print(ff)
    q = ff.replace('\\','_').replace('/','_').replace('_torus', '-torus').split('_')
    data = json.loads(open(ff).read())
    now = q[2]
    if now not in seen:
        draw_distr(keystoint(data[0]), f'Degree distribution / {now}', f'undir_degrees_{now}.png')
        seen.add(now)
    now = q[3]
    if now not in seen:
        draw_distr(keystoint(data[-1]), f'Degree distribution / {now}', f'undir_degrees_{now}.png')
        seen.add(now)

seen = set()
for ff in glob("graphstats/*_undir_spaths.json"):
    print(ff)
    q = ff.replace('\\','_').replace('/','_').replace('_torus', '-torus').split('_')
    data = json.loads(open(ff).read())
    now = q[2]
    if now not in seen:
        draw_distr(keystoint(data[0], 0), f'Shortest path distribution / {now}', f'undir_spaths_{now}.png')
        seen.add(now)
    now = q[3]
    if now not in seen:
        draw_distr(keystoint(data[-1], 0), f'Shortest path distribution / {now}', f'undir_spaths_{now}.png')
        seen.add(now)
    
