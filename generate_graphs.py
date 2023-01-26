import numpy as np
import itertools as it
import igraph as ig
import networkx as nx
import math

# Set parameters
n = 50
m = 4
# PA will have average degree 7.6
mean_deg = 7.6
p = mean_deg/(n-1)
in_out_ratio = 2
p_out = 2*mean_deg / (n+in_out_ratio * (n-2))
p_in = p_out*in_out_ratio
r = (p/np.pi)**0.5
thres_angle = 2*np.arcsin(p**0.5)


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
def torusdist(x,y):
    return sum([
        min((vx-vy)**2, (vx+1-vy)**2, (vy+1-vx)**2)
        for vx,vy in zip(x, y)
    ])**0.5
def p2torus_r(d, p):
    # Source: https://en.wikipedia.org/wiki/Volume_of_an_n-ball
    if d%2==0:
        return (math.factorial(int(d/2))*p)**(1/d) / np.pi**0.5
    # d is odd:
    k=int((d-1)/2)
    return (
        math.factorial(d)*p
        / (2*math.factorial(k)*(4*np.pi)**k)
    )**(1/d)

def p2thres_angle(d, p):
    """
        Computes the thres_angle such that a GRG on a d-sphere has density p.
        The surface area of a hyperspherical cap with angle alpha<pi/2 is a fraction
            (1/2) betainc(d/2,1/2.np.sin(alpha)**2)
        of the total surface of the d-sphere. Inverting this, gives
            np.arcsin(betaincinv(d/2, 1/2, 2*p)**0.5).
        If alpha>pi/2, then we can compute pi-alpha as the spherical cap with fraction 1-p.
    """
    from scipy.special import betaincinv
    if p == 0.5:
        return np.pi/2
    if p > 0.5:
        return np.pi - p2thres_angle(1-p, d=d)
    return np.arcsin(
        betaincinv(d/2, 1/2, 2*p) ** 0.5
    )


def generate_GRG_hypersphere(n=n, d=2, thres_angle=thres_angle, p=None, return_igraph=True):
    """
        This generates a GRG on a d-sphere (that is, d=1 is a circle while d=2 is a sphere).
        The link between thres_angle and the edge-density is tricky for higher dimensions. The edge-density is given by
        Incomplete Beta function:
        p=(1/2) * Beta(sin^2(thres_angle),d/2,1/2)
        For d=1, this is simply thres_angle/pi. For d=2, it is sin^2(thres_angle/2) and for d=3, we get
        (2*thres_angle-sin(2*thres_angle))/pi and it will only get worse from there.
    """
    if p is not None:
        thres_angle = p2thres_angle(d=d, p=p)
    coords = dict(zip(range(n), np.random.normal(size=(n, d+1))))
    edges = [
        (i, j) for i, j in it.combinations(range(n), 2)
        if angle(coords[i], coords[j]) < thres_angle
    ]
    if return_igraph:
        return edges2ig(n, edges)
    return edges


def generate_GRG_torus(n=n, r=r, p=None, d=2, return_igraph=True):
    if p is not None:
        r = p2torus_r(d=d, p=p)
    coords = dict(zip(range(n), np.random.rand(n, d)))
    edges = [
        (i, j) for i, j in it.combinations(range(n), 2)
        if torusdist(coords[i], coords[j]) < r
    ]
    if return_igraph:
        return edges2ig(n, edges)
    return edges


def generate_GRG_square(n=n, r=r):
    return ig.Graph.GRG(n, r, torus=False)


def generate_PA(n=n, m=m):
    """
        Note that the average degree will be lower than 2*m since initially, when m is larger than the number of
        vertices, it will not be possible to connect to m others. Also, this implementation seems to avoid multi-edges,
        so that the total number of edges is fixed and equal to m*(n-(m+1)/2).
    """
    return ig.Graph.Barabasi(n, m)


def edges2ig(n, edges):
    g = ig.Graph()
    g.add_vertices(n)
    g.add_edges(edges)
    return g


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
