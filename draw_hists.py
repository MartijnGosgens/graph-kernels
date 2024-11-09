from glob import glob
import json
from matplotlib import pyplot as plt
from matplotlib.ticker import MaxNLocator
import numpy as np

def draw_distr(vals, title, ofn):
    plt.clf()
    fig, ax = plt.subplots(figsize=(5,4))
    vals = {
        int(n): v
        for n,v in vals.items()
    }
    if 'path' in title:
        vals[0]=0
    n1,v1 = zip(*sorted(vals.items()))
    # width = .3
    ax.bar(np.array(n1), np.array(v1)/sum(v1))
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    plt.suptitle(title)
    # plt.legend()
    plt.tight_layout()
    print(ofn)
    plt.savefig(ofn)
    # plt.show() 
    
# _ = draw_distr(degrees[0], 'Degrees', '_test.png')


seen = set()
def keystoint(x):
    return {int(k): v for k, v in x.items()}

for ff in glob("graphstats/*_degrees.json"):
    print(ff)
    q = ff.replace('\\','_').replace('/','_').replace('_torus', '-torus').split('_')
    data = json.loads(open(ff).read())
    now = q[2]
    if now not in seen:
        draw_distr(keystoint(data[0]), f'Degree distribution / {now}', f'degrees_{now}.svg')
        seen.add(now)
    now = q[3]
    if now not in seen:
        draw_distr(keystoint(data[-1]), f'Degree distribution / {now}', f'degrees_{now}.svg')
        seen.add(now)

seen = set()
for ff in glob("graphstats/*_spaths.json"):
    print(ff)
    q = ff.replace('\\','_').replace('/','_').replace('_torus', '-torus').split('_')
    data = json.loads(open(ff).read())
    now = q[2]
    if now not in seen:
        draw_distr(keystoint(data[0]), f'Shortest path distribution / {now}', f'spaths_{now}.svg')
        seen.add(now)
    now = q[3]
    if now not in seen:
        draw_distr(keystoint(data[-1]), f'Shortest path distribution / {now}', f'spaths_{now}.svg')
        seen.add(now)
    
