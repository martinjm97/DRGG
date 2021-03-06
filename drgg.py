import numpy as np
from scipy import special
import matplotlib.pyplot as plt
import networkx as nx
from tqdm import tqdm

feat = ['n', 'alpha', 'd', 'edges/nlogn', 'dist/logn',
        'diameter/logn', 'triangles', 'clustering']


class Simulator:
    def __init__(self, n=100, d=4, alpha=6):
        self.n = n
        self.d = d
        self.alpha = alpha
        self.V = lambda d0: (np.pi)**(d0/2.0)/special.gamma(1 + d0/2.0)
        self.r0 = (np.log(n)/(n*self.V(d)))**(1.0/d)
        self.eta = (alpha - 1)*(self.r0)**(alpha - 1) / \
            (1 - (2*self.r0)**(alpha - 1))

    def generate_points(self, n0):
        """ Sample vertices. """
        return np.random.rand(n0, self.d)

    def generate_radii(self, n0):
        """ Obtain the radii of each of the vertices sampled. """

        def f(x): return (self.r0**(1-self.alpha) -
                          ((self.alpha - 1)*x/(self.eta)))**(1.0/(1-self.alpha))
        vf = np.vectorize(f)
        return vf(np.random.rand(n0))

    def distance(self, x, y):
        """ Compute the pairwise distances. """
        a = abs(x - y)

        def temp(a): return min(a, 1-a)
        vtemp = np.vectorize(temp)
        b = vtemp(a)
        return sum(b*b)**0.5

    def adjacency_lists(self, ptList, radList, verbose=False):
        """ Generate the adjacency lists for indegree and outdegree. """
        result1 = [set() for i in range(len(ptList))]
        result2 = [set() for i in range(len(ptList))]
        if verbose:
            print("Creating adjacency lists. ")
        for inVertex in tqdm(range(len(ptList)), disable=not verbose):
            for outVertex in range(len(ptList)):
                if inVertex == outVertex:
                    continue
                if self.distance(ptList[inVertex, :], ptList[outVertex, :]) < radList[inVertex]:
                    result1[inVertex].add(outVertex)
                    result2[outVertex].add(inVertex)
        return (result1, result2)

    def indegrees(self, inadj):
        """ From the indegree adjacency list produce the indegree distribution"""
        return [len(inadj[i]) for i in range(len(inadj))]

    def outdegrees(self, outadj):
        """ From the outdegree adjacency list produce the indegree distribution"""
        return [len(outadj[i]) for i in range(len(outadj))]

    def empirical_in_out(self):
        """ Find the indegree and outdegree distributions sampling as 
            appropriate.    
        """
        points = self.generate_points(self.n)
        radii = self.generate_radii(self.n)
        inList, outList = self.adjacency_lists(points, radii)
        indeg = self.indegrees(inList)
        outdeg = self.outdegrees(outList)
        return (indeg, outdeg)

    def adjacency_matrix(self, outList, verbose=False):
        """ Produce the adjacency matrix. """
        matrix = np.zeros((self.n, self.n))
        if verbose:
            print("Creating adjacency matrix. ")
        for i in tqdm(range(self.n), disable=not verbose):
            for j in outList[i]:
                matrix[i][j] = 1
        return matrix

    def information(self, verbose=True):
        """ Provide overall statistics of the graph produced. """
        points = self.generate_points(self.n)
        radii = self.generate_radii(self.n)
        inList, outList = self.adjacency_lists(points, radii, verbose)
        matrix = self.adjacency_matrix(outList)
        G = nx.from_numpy_matrix(matrix, create_using=nx.DiGraph())
        H = G.to_undirected()
        if nx.is_strongly_connected(G):
            return True, np.array([self.n,
                                   self.alpha,
                                   self.d,
                                   matrix.sum()/(self.n*np.log(self.n)),
                                   self.expected_distance(G)/np.log(self.n),
                                   nx.diameter(G)/np.log(self.n),
                                   self.sum_values(nx.triangles(H))/3,
                                   nx.average_clustering(H)
                                   ])
        else:
            return False, np.zeros(len(feat))

    def sum_values(self, d):
        try:
            return sum(b for a, b in d.items())
        except:
            return sum(b for a, b in d)

    def expected_distance(self, G):
        """ Find the expected shortest path length between vertices 
            i.e. how small world is the graph.
        """
        d = dict(nx.all_pairs_shortest_path_length(G))
        s, c = 0, 0
        for i in d:
            for j in d[i]:
                if d[i][j] > 0:
                    s += d[i][j]
                    c += 1
        return s / c

    def hist(self, lst):
        """ Produce an example histogram. """
        hs = np.histogram(lst, bins=50, range=(0, 50), density=True)
        return hs

    def plot(self, x, first, second=None):
        """ Generate a plot from the model. """
        if second == None:
            plt.plot(x, first)
            plt.show()
        else:
            plt.plot(x, first, 'r', second, 'b')
            plt.show()

    def concat_empirical(self):
        resultIn = []
        resultOut = []
        for k in range(repl):
            intemp, outtemp = empirical_in_out()
            resultIn.extend(intemp)
            resultOut.extend(outtemp)
        return (resultIn, resultOut)

    def theoretical_out(self):
        """ The expected outdegree distribution with the theoretical model. """
        z = eta*V(d)/(d - alpha + 1) * \
            (1/2.0**(d-alpha + 1) - r0**(d-alpha + 1))

        def f(k): return special.binom(self.n-1, k)*(z**k)*(1-z)**(self.n-1-k)
        return [f(k) for k in range(0, 51)]


if __name__ == "__main__":
    T = 10
    alpha = 4
    dim = 2
    lin = [100, 200, 300, 400, 500, 1000]
    data_mean = np.zeros((len(lin), 8))
    data_vars = np.zeros((len(lin), 8))
    for i in range(len(lin)):
        data = np.zeros((T, 8))
        for t in range(T):
            worked = False
            while not worked:
                s = Simulator(n=lin[i], alpha=alpha, d=dim)
                worked, data[t] = s.information(
                    printer='Testing n = '+str(lin[i])+'\t|  '+'('+str(t+1)+'/'+str(T)+')')
            assert data[t].max() > 0
        data_mean[i] = data.mean(axis=0)
        data_vars[i] = data.var(axis=0)
        print('Completed n = '+str(lin[i]), data_mean[i])

    np.savetxt("numerical_means.csv", data_mean, delimiter=",")
    np.savetxt("numerical_vars.csv", data_vars, delimiter=",")
