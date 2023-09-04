from generate_graphs import ig2edges,edges2grakel
from grakel.kernels import (RandomWalk,
                            GraphletSampling,
                            PyramidMatch,
                            NeighborhoodHash,
                            ShortestPath,
                            WeisfeilerLehman,
                            Propagation,
                            OddSth,
                            WeisfeilerLehmanOptimalAssignment,
                            NeighborhoodSubgraphPairwiseDistance)
from other_kernels import NetLSD,Gin
'''
    The kernels need not be grakel kernels, but they need to follow the same interface. That is, it should be a class
    (so that it has a .__name__), it's constructor needs to take the boolean parameter normalize and it needs to
    have a function fit_transform. That is, we should be able to call
        kernel(normalize=True,**kernel_params[kernel]).fit_transform(pack), 
    and this should return a 2D array of dimension (len(pack),len(pack)).
'''
import json
import itertools as it
from collections import defaultdict
from time import time

selected_kernels = (
     #RandomWalk, # ERRRs,
     #GraphletSampling,
     NetLSD,
     Gin,
     PyramidMatch,
     NeighborhoodHash,
     ShortestPath,
     WeisfeilerLehman,
     Propagation,
     OddSth,
     WeisfeilerLehmanOptimalAssignment,
     NeighborhoodSubgraphPairwiseDistance,
     # SvmTheta, # ERRR
    )

kernel_params = defaultdict(dict)
kernel_params[GraphletSampling] = {
    'sampling': {'n_samples': 500}
}
kernel_params[RandomWalk] = {
    'lamda': 0.1, # not 'lambda' (typo in grakel?)
    'p': 5
}

def tuple2str(t,sep=', '):
    return sep.join(map(str,t))


def nested_map(iterator,func,target_type=None,depth=-1):
    if target_type is not None and isinstance(iterator,target_type) or depth==0:
        return func(iterator)
    if hasattr(iterator,'items'):
        return {
            k: nested_map(val,func,target_type=target_type,depth=depth-1)
            for k,val in iterator.items()
        }
    elif hasattr(iterator, '__iter__'):
        return [
            nested_map(val,func,target_type=target_type,depth=depth-1)
            for val in iterator
        ]
    return func(iterator)

def locate(search_dict,locator_tuple):
    if len(locator_tuple)==1:
        return search_dict[locator_tuple[0]]
    return locate(search_dict[locator_tuple[0]],locator_tuple[1:])

def calc_mmd(vals,pack1size=None):
    if pack1size is None:
        pack1size = int(len(vals)/2)
    pack2size = len(vals)-pack1size
    idxs1 = range(pack1size)
    idxs2 = range(pack1size,pack1size+pack2size)
    return (
        (2/(pack1size*(pack1size-1)))*sum([vals[i,j] for i,j in it.combinations(idxs1,2)])
        +(2/(pack2size*(pack2size-1)))*sum([vals[i,j] for i,j in it.combinations(idxs2,2)])
        -(2/(pack1size*pack2size))*sum(vals[i,j] for i,j in it.product(idxs1,idxs2))
    )

class Experiment:
    '''
        generators: list of functions that return an ig graph
        parameters (optional): list of keyword-dictionaries for the parameters. Each parameter-tuple should be suitable for each of the generators

    '''
    def __init__(self,generators,parameters=None,npacks=None,sample_size=None):
        self.generators = generators
        self.generator_names = list(map(lambda g: g.__name__, generators))
        self.parameters = parameters
        if parameters is not None:
            self.parameter_names = list(map(lambda param: tuple2str(param.values()),parameters))
        self.npacks = npacks
        self.sample_size = sample_size

    def generate_pack(self,generator,sample_size,parameter={}):
        return [
            ig2edges(generator(**parameter))
            for _ in range(sample_size)
        ]

    '''
        Generates the graphs using the generators and parameters.
        Saves the json to save_name.
        The npacks_dict allows to specify a different number of packs for certain graph-parameter combinations. It should be a nested dict like
            npacks_dict[generator.__name__][parameter tuple string]
            If there are no parameters, the second index should be ''.
        The result will be a nested dictionary.
            If parameters is None:
                graphs[generator:self.generators][pack:int][isample:int]
            If parameters is not None:
                graphs[generator:self.generators][parameter:self.parameters][pack:int][isample:int]
    '''
    def generate_graphs(self,save_name,npacks=None,sample_size=None,npacks_dict=None):
        if npacks is None:
            npacks = self.npacks
        npacks_d = defaultdict(lambda: defaultdict(lambda: npacks))
        if npacks_dict is not None:
            for g,g_d in npacks_dict.items():
                for p,v in g_d.items():
                    npacks_d[g][p] = v
        if sample_size is None:
            sample_size = self.sample_size
        graphs = dict()
        for generator in self.generators:
            if self.parameters is None:
                graphs[generator.__name__] = [
                    self.generate_pack(generator,sample_size)
                    for _ in range(npacks_d[generator.__name__][''])
                ]
            else:
                graphs[generator.__name__] = dict()
                for parameter in self.parameters:
                    npacks_to_create = npacks_d[generator.__name__][tuple2str(parameter.values())]
                    print('Generating',npacks_to_create,'packs of',generator.__name__,'parameter',tuple2str(parameter.values()))
                    graphs[generator.__name__][tuple2str(parameter.values())] = [
                        self.generate_pack(generator,sample_size,parameter=parameter)
                        for _ in range(npacks_to_create)
                    ]
        with open(save_name, 'w') as save_file:
            json.dump(graphs, save_file)
        self.graphs = self.convert_graphs2grakel(graphs)
        

    def load_graphs(self,file_name):
        with open(file_name) as f:
            graphs = json.load(f)
        self.graphs = nested_map(graphs,edges2grakel,depth=(3 if self.parameters is None else 4))
    

    def convert_graphs2grakel(self,graphs):
        # Convert to grakel graph objects
        return nested_map(graphs,edges2grakel,depth=(3 if self.parameters is None else 4))


    def apply_kernels(self,comparison_iterator,save_name,kernels=selected_kernels,save_mmds_name=None,return_mmds=True):
        mmds = {k.__name__: defaultdict(list) for k in kernels}
        with open(save_name, 'w') as save_file, open(save_mmds_name+'.tsv','w') as mmds_tsv:
            for locator1,locator2 in comparison_iterator:
                pack1 = locate(self.graphs,locator1)
                pack2 = locate(self.graphs,locator2)
                for k in kernels:
                    #print(k.__name__)
                    start_time = time()
                    vals = k(normalize=True,**kernel_params[k]).fit_transform(pack1+pack2)
                    it1 = range(len(pack1))
                    it2 = range(len(pack1),len(pack1)+len(pack2))
                    for i,j in it.product(it1,it2):
                        print("\t".join(map(str, ['#', k.__name__]+list(locator1)+list(locator2)+[i,j-len(pack1),vals[i,j]])), flush=True,file=save_file)
                    for i,j in it.combinations(it1,2):
                        print("\t".join(map(str, ['#', k.__name__]+list(locator1)+list(locator1)+[i,j,vals[i,j]])), flush=True,file=save_file)
                    for i,j in it.combinations(it2,2):
                        print("\t".join(map(str, ['#', k.__name__]+list(locator2)+list(locator2)+[i-len(pack1),j-len(pack1),vals[i,j]])), flush=True,file=save_file)
                    mmd=calc_mmd(vals,self.sample_size)
                    mmds[k.__name__][locator1[:-1],locator2[:-1]].append(mmd)
                    print("\t".join(['#',k.__name__,tuple2str(locator1),tuple2str(locator2),str(mmd)]),flush=True,file=mmds_tsv)
        if save_mmds_name is not None:
            print(list(mmds.keys()))
            with open(save_mmds_name,'w') as mmd_file:
                json.dump({k: {
                    tuple2str(locator1)+'_vs_'+tuple2str(locator2): vals
                    for (locator1,locator2),vals in kmmds.items()
                } for k,kmmds in mmds.items()}, mmd_file)
        if return_mmds:
            return mmds
                    
                    
    def iterator_compare_generators(self):
        # For each two models, we compare the two packs of the same index
        for m1,m2 in it.combinations(self.generator_names,2):
            for p_idx in range(int(self.npacks/3)):
                yield (m1,p_idx),(m2,p_idx)
        # For each model, we also compare each pack to the next pack
        for m in self.generator_names:
            for p_idx in range(int(self.npacks/3)):
                yield (m,int(self.npacks/3)+p_idx),(m,int(self.npacks*2/3)+p_idx)
    
    def iterator_transitions(self):
        start_param = self.parameter_names[0]
        end_param = self.parameter_names[-1]
        for m in self.generator_names:
            for param in self.parameter_names:
                for p_idx in range(self.npacks):
                    # Comparison to start
                    if param!=start_param:
                        yield (m,param,p_idx),(m,start_param,p_idx) 
                    # Comparison to end
                    if param!=end_param:
                        yield (m,param,p_idx),(m,end_param,p_idx) 
                    # Comparison to the next pack of the same parameter
                    yield (m,param,p_idx),(m,param,(p_idx+1) % self.npacks) 

    def iterator_transitions_startcomparison(self):
        start_param = self.parameter_names[0]
        for m in self.generator_names:
            for i,param in enumerate(self.parameter_names):
                for p_idx in range(self.npacks):
                    # Comparison to start
                    if param!=start_param:
                        print('Compare',p_idx,'to',i*self.npacks+p_idx)
                        yield (m,param,p_idx),(m,start_param,i*self.npacks+p_idx)
                    else:
                        # Compare the first self.npacks of startparam to the last self.npacks of startparam
                        # Note that they haven't been used in the comparison to intermediate params
                        print('Compare',p_idx,'to',p_idx+len(self.parameter_names)*self.npacks)
                        yield (m,param,p_idx),(m,param,p_idx+len(self.parameter_names)*self.npacks) 

    def iterator_transitions_endcomparison(self):
        end_param = self.parameter_names[-1]
        for m in self.generator_names:
            for i,param in enumerate(self.parameter_names):
                for p_idx in range(self.npacks):
                    # Comparison to start
                    if param!=end_param:
                        yield (m,param,p_idx),(m,end_param,i*self.npacks+p_idx) 
                    else:
                        # Compare the first self.npacks of endparam to the last self.npacks of endparam
                        # Note that they haven't been used in the comparison to intermediate params
                        yield (m,param,p_idx),(m,param,p_idx+len(self.parameter_names)*self.npacks)
