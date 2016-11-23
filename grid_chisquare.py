#!/usr/bin/env python

# Grid Chi-square
# Module to perform Chi-square analysis for a grid of values.
from __future__ import division, print_function
import os
from joblib import Memory
from joblib import Parallel, delayed
import numpy as np
import matplotlib.pyplot as plt

path = "/home/jneal/Phd/Codes/equanimous-octo-tribble/"  # save path
cachedir = os.path.join(path, "cache")  # save path
memory = Memory(cachedir=cachedir, verbose=0)


# @memory.cache
def chi_squared(observed, expected, error=None):
    """Calculate chi squared.
    Same result as as scipy.stats.chisquare
    """
    if np.any(error):
        chisqr = np.sum((observed - expected) ** 2 / (error ** 2))
    else:
        # chisqr = np.sum((observed-expected)**2)
        chisqr = np.sum((observed - expected)**2 / expected)
        # When divided by exted the result is identical to scipy
    return chisqr


# @memory.cache
def grid_chisqr(t, obs, model, As, Bs, Cs):
    n_jobs = 4
    grid = Parallel(n_jobs=n_jobs)(delayed(chi_squared)(obs,
                                     model(t, a, b, c))
                                     for a in As for b in Bs for c in Cs)
    return grid


# @memory.cache
def model_func(t, a, b, c):
    return a * np.exp(-(t-b)**2 / (2*c**2))


if __name__ == "__main__":
    a = 5
    b = 5
    c = 2
    t = np.linspace(0, 5 * np.pi, 1000)
    iterable_a = np.linspace(0.1, 20, 50)
    iterable_b = np.linspace(0.1, 20, 30)
    iterable_c = np.arange(1, 4)
    simulation = model_func(t, a, b, c)

    chisqr_grid = parallel_chisqr_3D(t, simulation, model_func, iterable_a, iterable_b, iterable_c)

    print(chisqr_grid)
    X, Y, Z = np.meshgrid(iterable_a, iterable_b, iterable_c, indexing="ij")

    min_loc = np.argmin(chisqr_grid)
    print("min location", min_loc)

    a_sol = X.ravel()[min_loc]
    b_sol = Y.ravel()[min_loc]
    c_sol = Z.ravel()[min_loc]

    plt.plot(t, simulation, "x", label="sim")
    plt.plot(t, model_func(t, a_sol, b_sol, c_sol), label="solution")
    plt.legend()
    plt.show()
    print("Found Solution", a_sol, b_sol, c_sol)
    # main()

    plt.contourf(X[:, :, 1], Y[:, :, 1],
                 np.array(np.log10(chisqr_grid)).reshape(X.shape)[:, :, 1], 20)
    plt.show()
    # print("Time to run = {} seconds".format(time.time()-start))
