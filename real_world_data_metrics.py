import networkx as nx
import matplotlib.pyplot as plt
from collections import Counter

def pipeline(path):
    """ Run the whole pipline. """
    graph = nx.read_pajek(path) 

    # trimmed in and outdegree
    indegree, outdegree = prune_zero_outdegree(graph)

    # clustering coefficient, diameter, expected shortest path lengths
    clustering, diameter, small_world = get_metrics(graph)
    print(clustering, diameter, small_world)
    # values for PairsP 0.119 7 4.004

    # plot the indegree and outdegree histograms
    plot(indegree, outdegree)

def prune_zero_outdegree(graph):
    """ Remove the nodes with zero outdegree. """
    indegree = list(graph.in_degree().values())
    outdegree = list(graph.out_degree().values())
    outdegree = [k for k in outdegree if k>0]
    return indegree, outdegree

def get_metrics(graph):
    """ Compute important metrics of the data. """
    # we begin with a directed multigraph
    # we make undirected and remove multiedges
    udg = nx.Graph(graph.to_undirected())
    clustering = nx.average_clustering(udg)
    # use the giant component for diameter  
    # and average shortest path
    giant = max(nx.connected_component_subgraphs(udg), key=len)
    diameter = nx.diameter(giant)
    small_world = nx.average_shortest_path_length(giant)
    return clustering, diameter, small_world

def plot(indegree, outdegree):
    """ Display the results. """
    plt.figure(1)
    x,bins,count =plt.hist(indegree,  color = 'blue', alpha = 0.5, bins=33)

    plt.title("Indegree")
    plt.figure(2)
    plt.hist(outdegree, color = 'red', alpha = 0.5, bins=33)
    plt.title("Outdegree")
    plt.show()

if __name__ == "__main__":
    # load in the appropriate dataset to test.
    f = 'PairsP.net'
    pipeline(f)
