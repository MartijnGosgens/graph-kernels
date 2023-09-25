from glob import glob
import json
from matplotlib import pyplot as plt
import numpy as np

def draw_distr(vals, title, ofn):
    plt.clf()
    n1,v1 = zip(*sorted(vals.items()))
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
def keystoint(x):
    return {int(k): v for k, v in x.items()}

for ff in glob("graphstats/*_degrees.json"):
    print(ff)
    q = ff.replace('\\','_').replace('/','_').replace('_torus', '-torus').split('_')
    data = json.loads(open(ff).read())
    now = q[2]
    if now not in seen:
        draw_distr(keystoint(data[0]), f'Degree distribution / {now}', f'degrees_{now}.png')
        seen.add(now)
    now = q[3]
    if now not in seen:
        draw_distr(keystoint(data[-1]), f'Degree distribution / {now}', f'degrees_{now}.png')
        seen.add(now)

seen = set()
for ff in glob("graphstats/*_spaths.json"):
    print(ff)
    q = ff.replace('_torus', '-torus').split('_')
    data = json.loads(open(ff).read())
    now = q[1]
    if now not in seen:
        draw_distr(keystoint(data[0]), f'Shortest path distribution / {q[1]}', f'spaths_{now}.png')
        seen.add(now)
    now = q[2]
    if now not in seen:
        draw_distr(keystoint(data[-1]), f'Shortest path distribution / {q[1]}', f'spaths_{now}.png')
        seen.add(now)
    
