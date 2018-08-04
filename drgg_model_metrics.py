import networkx as nx
from ParetoRGG import Simulator
from pajek_data import get_metrics


def simulated_graph(n, d, alpha):
    """ Create a new directed geometric graph with the given params. """
    simulated_graph = Simulator(n, d, alpha)
    radii = simulated_graph.generate_points(simulated_graph.n)
    vertices = simulated_graph.generate_radii(simulated_graph.n)
    in_list, out_list = simulated_graph.adjancency_lists(vertices, radii)
    matrix = simulated_graph.adjacency_matrix(out_list)
    g = nx.from_numpy_matrix(matrix, create_using=nx.DiGraph())
    return g


# use the number of nodes, d, and alpha
# values used to fit the PairsP data
iters = 100
acs = []   # average clustering
ds = []    # diameter of giant component
asps = []  # avg shortest paths
# sample from a number of graphs
for i in range(iters):
    g = simulated_graph(n=5019, d=3, alpha=8)
    avg_clustering, diameter, avg_shortest_path = get_metrics(g)
    acs.append(avg_clustering)
    ds.append(diameter)
    asps.append(avg_shortest_path)
print(sum(acs) / iters, sum(ds) / iters, sum(asps) / iters)
print(acs)
print(ds)
print(asps)
