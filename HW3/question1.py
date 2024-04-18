# This is a sample Python script.
import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import solve
# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.

#the u can satisfy the boundry condition
u_exact = lambda x: np.sin(np.pi * x)
#the forcing function f
f = lambda x: (np.pi ** 2) * np.sin(np.pi * x)

#compute the Chebyshev differentiation matrix

def chebyshev_points(N):
    #Compute the Chebyshev points of the second kind.
    return np.cos(np.pi * np.arange(N+1) / N)
def cheb_matrix(x):
    #x is the cheb points
    N = len(x) - 1
    c = np.hstack([2, np.ones(N-1), 2]) * (-1)**np.arange(N+1)
    X = np.tile(x, (N+1, 1)).T
    dX = X - X.T
    D = np.outer(c, 1/c) / (dX + np.eye(N+1))
    return D - np.diag(np.sum(D, axis=1))


def solve_poisson(N):
    #Solve the Poisson equation using the Chebyshev spectral method.
    x = chebyshev_points(N)
    D2 = cheb_matrix(x) @ cheb_matrix(x)  # Second derivative matrix
    # Apply boundary conditions to matrix
    D2[0, :] = 0; D2[N, :] = 0
    D2[0, 0] = 1; D2[N, N] = 1

    f = -2 * np.ones_like(x)  # Forcing term
    # Apply boundary conditions to f
    f[0] = f[N] = 0

    u = solve(D2, f)
    return x, u



# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    N_values = np.geomspace(10, 1000, 10).astype(int)
    errors_l2 = []
    errors_linf = []

    for N in N_values:
        x, u_numerical = solve_poisson(N)
        u_exact = x ** 2 - 1  # Exact solution
        error_l2 = np.sqrt(np.sum((u_numerical - u_exact) ** 2) / np.sum(u_exact ** 2))
        error_linf = np.max(np.abs(u_numerical - u_exact)) / np.max(np.abs(u_exact))
        errors_l2.append(error_l2)
        errors_linf.append(error_linf)

    print(errors_l2)

    # Plotting
    plt.figure(figsize=(10, 6))
    plt.loglog(1 / N_values, errors_l2, 'o-', label='$L^2$ Error')
    plt.loglog(1 / N_values, errors_linf, 's-', label='$L^\\infty$ Error')
    plt.xlabel('$h$ (Approximate)')
    plt.ylabel('Relative Error')
    plt.legend()
    plt.grid(True)
    plt.title('Spatial Convergence Study of the Chebyshev Spectral Method')
    plt.show()

# See PyCharm help at https://www.jetbrains.com/help/pycharm/
