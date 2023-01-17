import numpy as np
import igraph as ig
import networkx as nx

# Set parameters
n = 50
mean_deg = 7
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
