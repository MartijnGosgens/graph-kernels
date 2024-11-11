from generate_graphs import (interpolate_ER_triangular, interpolate_ER_PPM, interpolate_ER_GRG_torus, interpolate_ER_inhomogeneous, interpolate_GRG_torus_circle, interpolate_ER_GCG)
from collections import defaultdict
from os import listdir,path

folder = ''
interpolators = ['interpolate_ER_inhomogeneous','interpolate_ER_PPM',  'interpolate_ER_GRG_torus', 'interpolate_GRG_torus_circle', 'interpolate_ER_GCG']

k2display={
    "DegreeHistogram": "Degree",
    "ShortestPath": "SP",
    "WeisfeilerLehman": "WL",
    "WeisfeilerLehmanOptimalAssignment": "WL-OA",
    "GraphletSampling": "Graphlet-3",
    "GraphletSampling4": "Graphlet-4",
    "NeighborhoodSubgraphPairwiseDistance": "NSPDK",
    "PyramidMatch": "PM",
    "NetLSD": "NetLSD-heat",
    "NetLSDWave": "NetLSD-wave",
    "Gin": "RandGIN",
}

interpolation2display = {
    #"interpolate_ER_density": "Density",
    "interpolate_ER_inhomogeneous": "Heterogeneity",
    "interpolate_ER_PPM": "Communities",
    "interpolate_ER_GRG_torus": "Geometry",
    "interpolate_GRG_torus_circle": "Dimensionality",
    "interpolate_ER_GCG": "Complementarity"
}

k2interpolator2result = defaultdict(dict)
interpolator2k2result = defaultdict(dict)

def tsv2dict(tsv_file,invert=False):
    output = {}
    bests = -1
    for line in tsv_file.readlines():
        output[line.split('\t')[1]] = float(line.split('\t')[-1].strip()) * (-1 if invert else 1)
    return output,max(output,key=output.get)

# Table of running times
import json
k2inter2time = defaultdict(dict)
for interpolation in interpolators:
    with open(path.join(folder,interpolation+'_start_times.json')) as start_file,open(path.join(folder,interpolation+'_end_times.json')) as end_file:
        start_times = json.load(start_file)
        end_times = json.load(end_file)
        for k in k2display.keys():
            ktimes = sum(start_times[k].values(),[])+sum(end_times[k].values(),[])
            k2inter2time[k][interpolation] = sum(ktimes)/len(ktimes)
with open('times.txt','w') as outfile:
    print(r'\begin{tabular}{c|'+'c'*len(interpolators)+r'}',file=outfile)
    print(r'\hline',file=outfile)
    print('\t&\t'.join(['']+[
        interpolation2display[interpolator]
        for interpolator in interpolators
    ])+r'\\',file=outfile)
    for k,display in k2display.items():
        print('\t&\t'.join([k2display[k]]+[
            '${:.2f}$'.format(k2inter2time[k][interpolator])
            for interpolator in interpolators
        ])+r'\\', file=outfile)
    print(r'\end{tabular}',file=outfile)

            



for fn in listdir(folder):
    with open(path.join(folder,fn)) as file:
        end = ('_end_' in fn)
        for line in file.readlines():
            k2interpolator2result[line.split('\t')[1]][(line.split('\t')[2],end)] = float(line.split('\t')[-1].strip()) * (-1 if end else 1)
            interpolator2k2result[(line.split('\t')[2],end)][line.split('\t')[1]] = float(line.split('\t')[-1].strip()) * (-1 if end else 1)

bests = {
    interpolator: max(k2result,key=k2result.get)
    for interpolator,k2result in interpolator2k2result.items()
}
kernels = list(k2display.keys())
interpolators = list(interpolation2display.keys())
interpolator2k2average = {
    interpolator: {
        k: (interpolator2k2result[(interpolator,False)][k]+interpolator2k2result[(interpolator,True)][k])/2
        for k in kernels
    }
    for interpolator in interpolators
}
best_avgs = {
    interpolator: max(k2result,key=k2result.get)
    for interpolator,k2result in interpolator2k2average.items()
}

results = defaultdict(dict)
average=True
best_start = {}
best_end = {}
hist = []
for interpolator in interpolation2display:
    with open(f'{interpolator}_start_spearmans.tsv') as startfile, open(f'{interpolator}_end_spearmans.tsv') as endfile:
        start_d,best_start[interpolator]=tsv2dict(startfile)
        end_d,best_end[interpolator]=tsv2dict(endfile,invert=True)
        for k in start_d.keys():
            results[k][interpolator] = (start_d[k],end_d[k])
            hist.append(start_d[k]-end_d[k])
with open('output.txt','w') as outfile,open('output2.txt','w') as outfile2:
    print(r'\begin{tabular}{c|'+'c'*len(interpolators)+r'}',file=outfile)
    print(r'\hline',file=outfile)
    print('\t&\t'.join(['']+[
        interpolation2display[interpolator]
        for interpolator in interpolators
    ])+r'\\',file=outfile)
    
    print('\t'.join(['']+[
        interpolation2display[interpolator]
        for interpolator in interpolators
    ]),file=outfile2)
    for k in kernels:
        if average:
            print('\t&\t'.join([k2display[k]]+[
                ('$\\textbf{'+'{:.3f}'.format(
                    interpolator2k2average[interpolator][k]
                )+'}$' if best_avgs[interpolator] == k else '${:.3f}$'.format(
                    interpolator2k2average[interpolator][k]))
                for interpolator in interpolators
            ])+r'\\', file=outfile)
            print('\t'.join([k2display[k]]+[
                '{:.3f}'.format(
                    interpolator2k2average[interpolator][k]
                )
                for interpolator in interpolators
            ]),file=outfile2)
            print('\t'.join(['']+[
                interpolation2display[interpolator]
                for interpolator in interpolators
            ]))
        else:
            print('\t&\t'.join([k2display[k]]+[
                ('$\\textbf{'+'{:.3f}'.format(
                    interpolator2k2result[(interpolator,False)][k]
                )+'}$' if bests[(interpolator,False)] == k else '${:.3f}$'.format(
                    interpolator2k2result[(interpolator,False)][k]
                ))+
                ' / '+
                ('$\\textbf{'+'{:.3f}'.format(
                    interpolator2k2result[(interpolator,True)][k]
                )+'}$' if bests[(interpolator,True)] == k else '${:.3f}$'.format(
                    interpolator2k2result[(interpolator,True)][k]
                ))
                for interpolator in interpolators
            ])+r'\\', file=outfile)
    print(r'\end{tabular}',file=outfile)

