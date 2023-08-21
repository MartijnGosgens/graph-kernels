from generate_graphs import (interpolate_ER_triangular, interpolate_ER_PPM, interpolate_ER_GRG_torus, interpolate_ER_inhomogeneous, interpolate_GRG_torus_circle, interpolate_ER_GCG)
from collections import defaultdict

interpolators = ['interpolate_ER_inhomogeneous','interpolate_ER_triangular', 'interpolate_ER_PPM', 'interpolate_ER_GRG_torus', 'interpolate_GRG_torus_circle', 'interpolate_ER_GCG']

def tsv2dict(tsv_file):
    output = {}
    bests = -1
    for line in tsv_file.readlines():
        output[line.split('\t')[1]] = float(line.split('\t')[-1].strip())
    return output,max(output,key=output.get)

def interpolator2label(interpolator):
    endpoints = interpolator.split('_')[1:]
    return endpoints[0]+r'$\leftrightarrow$'+endpoints[1]

results = defaultdict(dict)
best_start = {}
best_end = {}
for interpolator in interpolators:
    with open(f'{interpolator}_start_spearmans.tsv') as startfile, open(f'{interpolator}_end_spearmans.tsv') as endfile:
        start_d,best_start[interpolator]=tsv2dict(startfile)
        end_d,best_end[interpolator]=tsv2dict(endfile)
        for k in start_d.keys():
            results[k][interpolator] = (start_d[k],end_d[k])
with open('output.txt','w') as outfile:
    print(r'\begin{tabular}{c|'+'c'*len(interpolators)+r'}',file=outfile)
    print('\t&\t'.join(['']+[
        interpolator2label(interpolator)
        for interpolator in interpolators
    ])+r'\\',file=outfile)
    for k,result in results.items():
        print('\t&\t'.join([k]+[
            ('$\\textbf{'+'{:.3f}'.format(result[interpolator][0])+'}$' if best_start[interpolator] == k else '${:.3f}$'.format(result[interpolator][0]))+
            ' / '+
            ('$\\textbf{'+'{:.3f}'.format(result[interpolator][1])+'}$' if best_end[interpolator] == k else '${:.3f}$'.format(result[interpolator][1]))
            for interpolator in interpolators
        ])+r'\\', file=outfile)
    print(r'\end{tabular}',file=outfile)

