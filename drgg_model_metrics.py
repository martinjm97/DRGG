import networkx as nx
from drgg import Simulator
from real_world_data_metrics import get_metrics
from tqdm import tqdm


def simulated_graph(n, d, alpha, verbose=False):
    """ Create a new directed geometric graph with the given params. """
    simulated_graph = Simulator(n, d, alpha)
    vertices = simulated_graph.generate_points(simulated_graph.n)
    radii = simulated_graph.generate_radii(simulated_graph.n)
    in_list, out_list = simulated_graph.adjacency_lists(vertices, radii, verbose)
    matrix = simulated_graph.adjacency_matrix(out_list, verbose)
    g = nx.from_numpy_matrix(matrix, create_using=nx.DiGraph())
    return g

if __name__ == "__main__":
    # use the number of nodes, d, and alpha
    # values used to fit the PairsP data
    iters = 100
    acs = []   # average clustering
    ds = []    # diameter of giant component
    asps = []  # avg shortest paths

    print("Performing DRGG simulations: ")
    # sample from a number of graphs
    for i in tqdm(range(iters)):
        g = simulated_graph(n=5019, d=3, alpha=8, verbose=False)
        avg_clustering, diameter, avg_shortest_path = get_metrics(g)
        acs.append(avg_clustering)
        ds.append(diameter)
        asps.append(avg_shortest_path)
    print("Simulated averages for {} trials.".format(iters))
    print("Average clustering is ", sum(acs) / iters)
    print("Diameter of the giant component is ", sum(ds) / iters)
    print("Average shortest path length is ", sum(asps) / iters)
