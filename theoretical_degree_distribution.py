from scipy.integrate import quad
from scipy.special import binom, beta, gamma
import numpy as np
import matplotlib.pyplot as plt
from decimal import Decimal
from utils import wrap_decimal


@wrap_decimal
def integrand(r, n, k, d, vd, c, alpha):
    """ Computes the integrand using the Decimal object, which allows for
    arbitrary precision, useful for n > 300, but slows down computation
    significantly overall.  """
    nCk = Decimal(binom(float(n), float(k)))
    return float((c/r**alpha) * (vd*r**d)**k * (1-vd*r**d)**(n-k) * nCk)


def compute_m(n, vd, d):
    return (np.log(n) / (vd*n))**(1.0/d)


def compute_c(alpha, m):
    return (alpha-1) / (m**(1-alpha) - 1/2.0**(float(1-alpha)))


def compute_vd(d):
    return np.pi**(d/2.0) / gamma(1+d/2.0)


def exact_solution(n, k, d, vd, c, m, alpha):
    """ integrate the exact function from m (the min circle size) to 1/2 """
    return quad(lambda r: integrand(r, n, k, d, vd, c, alpha), m, 1/2.0)[0]


@wrap_decimal
def ramis_approx_solution(c, vd, n, k, alpha, d):
    beta = (alpha-1)/d + 1
    nCk = Decimal(binom(float(n), float(k)))
    numerator = c * nCk * vd**((alpha-1)/d)/d * Decimal(np.sqrt(2*np.pi)) * (
        k-beta)**(k-beta+Decimal(0.5)) * (n-k)**(n-k+Decimal(0.5))
    denomenator = (n-beta)**(n-beta+Decimal(1.5))
    return numerator / denomenator


def approximation_quality(exact_solutions, approx_solutions):
    print("exact sum: " + str(sum(exact_solutions)))
    print("appx sum: " + str(sum(approx_solutions)))
    exact = np.array(exact_solutions)
    apx = np.array(approx_solutions)
    mse = sum((exact-apx)**2)*1.0/len(exact)
    print("mean squared error: " + str(mse))


def compute_z(d, vd, c, m, alpha):
    alpha = float(alpha)
    return c * vd / (d-alpha+1) * (1.0/(2**(d-alpha+1)) - m**(d-alpha+1))


def out_deg_distro(n, k, d, vd, c, m, alpha):
    z = compute_z(d, vd, c, m, alpha)
    return binom(n-1, k) * z**k * (1 - z)**(n-k-1)


def compute_edge_count(alpha, n, d, vd, c, m):
    """ Expected total number of edges. """
    z = compute_z(d, vd, c, m, alpha)
    return (n - 1) * n * z


def create_plot(exact, apx, n):
    start = 5
    index = range(start, len(exact[:n]))
    axis_font = {'size': '16'}
    a, = plt.plot(np.log(index), np.log(exact[start:n]), '.', color='r')
    b, = plt.plot(np.log(index), np.log([float(i) for i in apx[start:n]]), '-', color='black')
    plt.rcParams.update({'font.size': 16})
    plt.xlabel('Log indegree', axis_font)
    plt.ylabel('Log probability', axis_font)
    plt.title('Indegree distribution for {} Nodes'.format(n))
    plt.legend([a, b], ['Exact', 'Approximation'])
    # plt.savefig('{}nodeindeg.eps'.format(n), format='eps', dpi=1000)
    plt.show()


def run_pipeline():
    n, d, alpha = (500, 2, 3)
    vd = compute_vd(d)
    m = compute_m(n, vd, d)
    c = compute_c(alpha, m)

    exact_solutions = []
    approx_solutions = []
    for k in range(2, n):
        exact_solutions.append(exact_solution(n, k, d, vd, c, m, alpha))
        approx_solutions.append(ramis_approx_solution(c, vd, n, k, alpha, d))
    # approximation_quality(exact_solutions, approx_solutions)
    create_plot(exact_solutions, approx_solutions, n)


if __name__ == "__main__":
    run_pipeline()
