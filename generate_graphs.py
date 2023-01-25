import numpy as np
import itertools as it
import igraph as ig
import networkx as nx

# Set parameters
n = 50
mean_deg = 8
p = mean_deg/(n-1)
p_out = p/2
p_in = p*((3/2)*n-1)/(n-1)
r = (p/np.pi)**0.5


# Graph generators
def generate_ER(n=n,p=p):
    return ig.Graph.Erdos_Renyi(n, p, directed=False, loops=False)


def generate_PPM(n=n,p_in=p_in,p_out=p_out,k=2):
    rem = n % k
    sizes = [int(n/k)+1]*rem+[int(n/k)]*(k-rem)
    ps = [
        [
            p_in if i==j else p_out
            for j in range(k)
        ]
        for i in range(k)
    ]
    return ig.Graph.SBM(n, ps, sizes, directed=False, loops=False)


def generate_GRG(n=n, r=r):
    return ig.Graph.GRG(n, r, torus=True)


def angle(x,y):
    return np.arccos((x*y).sum()/((x*x).sum() * (y*y).sum())**0.5)
def generate_GRG_hypersphere(n, d, thres_angle, return_igraph=True):
    """
        This generates a GRG on a d-sphere (that is, d=1 is a circle while d=2 is a sphere).
        The link between thres_angle and the edge-density is tricky for higher dimensions. The edge-density is given by
        Incomplete Beta function:
        p=(1/2) * Beta(sin^2(thres_angle),d/2,1/2)
        For d=1, this is simply thres_angle/pi. For d=2, it is sin^2(thres_angle/2) and for d=3, we get
        (2*thres_angle-sin(2*thres_angle))/pi and it will only get worse from there.
    """
    coords = dict(zip(range(n),np.random.normal(size=(n,d+1))))
    edges = [
        (i,j) for i,j in it.combinations(range(n),2)
        if angle(coords[i],coords[j])<thres_angle
    ]
    if return_igraph:
        g = ig.Graph()
        g.add_vertices(n)
        g.add_edges(edges)
        return g
    return edges

def generate_GRG_square(n=n, r=r):
    return ig.Graph.GRG(n, r, torus=False)

def ig2edges(g):
    return [(e.source, e.target) for e in g.es]

def edges2nx(e):
    return nx.from_edgelist(e)

def ig2nx(g):
    return edges2nx(ig2edges(g))

def nx2stats(ng):
    return (
        nx.average_clustering(ng), # avg_clustering
        len(ng.edges()), # edges
        nx.transitivity(ng), # transitivity
        np.sum(list(nx.triangles(ng).values()))/3, # triangles
        np.sum([y for x,y in ng.degree()])/len(ng) # avg_degree
    )
