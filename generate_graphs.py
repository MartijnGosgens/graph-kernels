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

def generate_inhomogeneous(n=n,p=p,pl_exp=3):
    gamma = 1/(pl_exp-1)
    fitness=np.exp(np.random.exponential(gamma,n))
    # The ig generator requires a fixed number of edges. We sample a Binomial number of edges, so that it resembles ER/PPM/GRG.
    n_edges=np.random.binomial(n*(n-1)/2,p)
    print('Generating inhomogeneous with exponent',pl_exp,'and',n_edges,'edges')
    return ig.Graph.Static_Fitness(n_edges,fitness)

def interpolate_ER_inhomogeneous(step,n=n,p=p):
    pl_exp = 1+1/step if step>0 else float('Inf')
    return generate_inhomogeneous(n=n,p=p,pl_exp=pl_exp)


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

def closure_graph(n,p,p1):
    G1=nx.erdos_renyi_graph(n,p1)
    p2 = (p-p1)/((1-p1)*(1-(1-p1**2)**(n-2)))
    G2=nx.erdos_renyi_graph(n,p2)
    wedges = [
        (i,j) 
        for (i,j) in it.combinations(G1.nodes,2) 
        if (i not in G1[j]) and (set(G1[i])&set(G1[j]))!=set()
    ]
    G1.add_edges_from([e for e in wedges if np.random.rand()<p2])
    return nx2ig(G1)

# Step needs to be in the interval [0,1], so that p_in=(1+step)*p_out
def interpolate_ER_PPM(step,p=p,n=n,k=2):
    in_out_ratio = 1+5.6*step
    p_out = 2*mean_deg / (n+in_out_ratio * (n-2))
    p_in = p_out*in_out_ratio
    return generate_PPM(n=n,p_in=p_in,p_out=p_out,k=k)

def interpolate_ER_triangular(step,p=p,n=n):
    p1 = p*(1-step/2)
    return closure_graph(n,p,p1)


def generate_GRG(n=n, r=r):
    return ig.Graph.GRG(n, r, torus=True)


def angle(x,y):
    return np.arccos((x*y).sum()/((x*x).sum() * (y*y).sum())**0.5)
def torusdist(x,y,sizes=None):
    if sizes is None:
        sizes = (1,)*len(x)
    return sum([
        min((px-py)**2, (px+s-py)**2, (py+s-px)**2)
        for px,py,s in zip(x, y,sizes)
    ])**0.5
def p2torus_r(d, p, h=1):
    if d==2 and h!=1:
        if h>4*p/np.pi:
            return (p*h/np.pi)**0.5
        elif h<=0:
            return p/2
        else:
            # We perform Newton-Raphson on f to find the intersection f(z)=2*p/h.
            # Note that target>pi/2 because of the previous if-statement
            target = 2*p/h
            f = lambda z: (1/z**2-1)**0.5+np.arcsin(z)/z**2
            f_prime = lambda z: -2*np.arcsin(z)/z**3
            # We start at z=1, though we could basically start anywhere in (0,1]
            z = 1
            val = f(z)
            tolerance = 0.01
            # Usually this takes less than 3-5 iterations
            while abs(val-target)>tolerance:
                step = (val-target)/f_prime(z)
                # Ensure z won't become negative
                if step>z:
                    step = 0.9*z
                z -= step
                val = f(z)
            return h/(2*z)

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


def generate_GCG_hypersphere(n=n, d=2, p=p, return_igraph=True):
    '''
        GCG stands for Geometric Complementarity Graph. Complementarity is modeled geometrically by connecting nodes that are close to being
        maximally distant (i.e., at opposite poles).
    '''
    thres_angle = np.pi-p2thres_angle(d=d, p=p)
    coords = dict(zip(range(n), np.random.normal(size=(n, d+1))))
    edges = [
        (i, j) for i, j in it.combinations(range(n), 2)
        if angle(coords[i], coords[j]) > thres_angle
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

def interpolate_ER_GRG_torus(step,p=p,n=n,d=2, return_igraph=True):
    if p is not None:
        r = p2torus_r(d=d, p=p)
    coords = dict(zip(range(n), np.random.rand(n, d)))
    edges = [
        (i, j) for i, j in it.combinations(range(n), 2)
        if (torusdist(coords[i], coords[j]) < r and np.random.rand()<p+step*(1-p)) or (torusdist(coords[i], coords[j]) > r and np.random.rand()<p*(1-step))
    ]
    if return_igraph:
        return edges2ig(n, edges)
    return edges

def interpolate_ER_GCG(step,p=p,n=n, return_igraph=True):
    thres_angle = np.pi-p2thres_angle(d=2, p=p)
    coords = dict(zip(range(n), np.random.normal(size=(n, 3))))
    edges = [
        (i, j) for i, j in it.combinations(range(n), 2)
        if (angle(coords[i], coords[j]) > thres_angle) and (np.random.rand()<p+step*(1-p)) or (angle(coords[i], coords[j]) < thres_angle) and (np.random.rand()<p-step*p)
    ]
    if return_igraph:
        return edges2ig(n, edges)
    return edges

def interpolate_GRG_torus_circle(step, p=p, n=n, return_igraph=True):
    h = 1-step
    r = p2torus_r(d=2, p=p, h=h)
    rands = np.random.rand(n, 2)
    coords = dict(zip(range(n), zip(rands[:,0],h*rands[:,1])))
    edges = [
        (i, j) for i, j in it.combinations(range(n), 2)
        if torusdist(coords[i], coords[j],sizes=(1,h)) < r
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
    return ig.Graph.Barabasi(n, m, outpref=True)


def edges2ig(n, edges):
    g = ig.Graph()
    g.add_vertices(n)
    g.add_edges(edges)
    return g

def nx2ig(g):
    return edges2ig(len(g),g.edges)


def ig2edges(g):
    return [(e.source, e.target) for e in g.es]


def edges2nx(e):
    G = nx.Graph()
    G.add_nodes_from(range(n))
    G.add_edges_from(e)
    return G

def edges2grakel(g,N=n):
    from grakel import Graph
    # Grakel assumes by default that graphs are directed
    edges = list(map(tuple,g))+list(map(tuple,map(reversed,g)))
    return Graph(edges, node_labels={i: 'A' for i in range(N)}, edge_labels={e: 'B' for e in edges})

def grakel2nx(g):
    return nx.from_edgelist(g.get_edges())

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
