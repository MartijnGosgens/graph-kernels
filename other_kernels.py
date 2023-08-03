
import numpy as np
from generate_graphs import grakel2nx

class Kernel:
    def __init__(self,normalize=True):
        self.normalize=normalize

    '''
        A function that maps a networkx graph to a numpy vector
    '''
    def embed(self,nxg):
        pass

    '''
        Returns the dotproduct if normalize==False and the cosine similarity otherwise
    '''
    def sim(self,v1,v2):
        dot = np.dot(v1,v2)
        if not self.normalize:
            return dot
        return dot / (np.linalg.norm(v1)*np.linalg.norm(v2))

    def fit_transform(self,pack):
        representations = None
        if self.embed_dataset is not None:
            representations = self.embed_dataset(map(grakel2nx,pack))
        else:
            representations = [
                self.embed(grakel2nx(grakelg))
                for grakelg in pack
            ]
        return np.array([
            [
                self.sim(v1,v2)
                for v2 in representations
            ]
            for v1 in representations
        ]) 


class NetLSD(Kernel):
    def __init__(self, normalize=True):
        super().__init__(normalize)
        import netlsd
        self.embed = netlsd.heat
        self.embed_dataset = None
        
class Gin(Kernel):
    def __init__(self, normalize=True):
        super().__init__(normalize)
        from scipy import spatial
        import torch
        import dgl
        from GgmMetrics.evaluation import gin_evaluation
        device =  torch.device('cpu')
        model = gin_evaluation.load_feature_extractor(device)
        ev = gin_evaluation.MMDEvaluation(model=model, kernel='rbf', sigma='range', multiplier='mean')
        self.embed = lambda x: ev._GINMetric__get_activations_single_dataset([dgl.DGLGraph(x)])[0]
        self.embed_dataset = lambda x: ev._GINMetric__get_activations_single_dataset(list(map(dgl.DGLGraph,x)))

    