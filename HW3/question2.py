import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import diags
from scipy.sparse.linalg import cg
from scipy.linalg import norm

# Constants
nu = 1  # Diffusion coefficient
pi = np.pi


# Manufactured solution and its derivatives
def u_manufactured(x, y, t):
    return np.sin(pi * x) * np.cos(pi * y) * np.exp(-pi * t)


def f_manufactured(x, y, t):
    return 2 * pi * np.sin(pi * x) * np.cos(pi * y) * np.exp(-pi * t) - (pi ** 2) * u_manufactured(x, y, t) * (2)


# Function to solve the heat equation
def solve_heat_equation(N, dt, T_end):
    # Grid setup
    x = np.linspace(0, 1, N + 2)
    y = np.linspace(0, 1, N + 2)
    h = 1 / (N + 1)
    r = nu * dt / (2 * h ** 2)

    # Time setup
    t_steps = int(T_end / dt) + 1

    # Initial condition
    U = np.zeros((N + 2, N + 2))
    for i in range(N + 2):
        for j in range(N + 2):
            U[i, j] = u_manufactured(x[i], y[j], 0)

    # Construct the sparse matrix for Crank-Nicolson scheme
    main_diag = (1 + 4 * r) * np.ones(N)
    off_diag = -r * np.ones(N - 1)
    A = diags([main_diag, off_diag, off_diag], [0, -1, 1])

    for _ in range(1, t_steps):
        # Construct the right-hand side
        F = np.zeros((N, N))
        for i in range(1, N + 1):
            for j in range(1, N + 1):
                F[i - 1, j - 1] = U[i, j] + dt / 2 * f_manufactured(x[i], y[j], (_ - 0.5) * dt)


        U_interior = np.zeros((N, N))
        for i in range(N):
            b = F[:, i] + r * (U[1:-1, i] + U[1:-1, i + 2])
            U_interior[:, i], _ = cg(A, b, tol=1e-8)

        # Update solution with boundary conditions
        U[1:-1, 1:-1] = U_interior
        for i in range(N + 2):
            U[i, 0] = u_manufactured(x[i], y[0], dt)
            U[i, -1] = u_manufactured(x[i], y[-1], dt)
            U[0, i] = u_manufactured(x[0], y[i], dt)
            U[-1, i] = u_manufactured(x[-1], y[i], dt)

    return x, y, U


if __name__ == '__main__':
    dt_values = np.array([0.1, 0.05, 0.01, 0.005, 0.001])
    errors_l2 = np.zeros(len(dt_values))
    errors_linf = np.zeros(len(dt_values))
    N_spatial = 500  # Fixed spatial discretization
    T_end = 0.1  # Fixed end time for the study


    for i, dt in enumerate(dt_values):
        x, y, U_numerical = solve_heat_equation(N_spatial, dt, T_end)
        U_exact = np.zeros_like(U_numerical)
        for j in range(len(x)):
            for k in range(len(y)):
                U_exact[j, k] = u_manufactured(x[j], y[k], T_end)
        error_l2 = norm(U_numerical - U_exact) / norm(U_exact)
        error_linf = np.max(np.abs(U_numerical - U_exact)) / np.max(np.abs(U_exact))
        errors_l2[i], errors_linf[i] = error_l2, error_linf

    plt.plot(dt_values, errors_l2, '-o', label='L2 Error')
    plt.xlabel('Δt')
    plt.ylabel('Relative Error')
    plt.title('Temporal Convergence Study (L2 Error)')
    plt.xscale('log')
    plt.yscale('log')
    plt.grid(True)
    plt.legend()
    plt.show()

    plt.plot(dt_values, errors_linf, '-o', label='Linf Error')
    plt.xlabel('Δt')
    plt.ylabel('Relative Error')
    plt.title('Temporal Convergence Study (Linf Error)')
    plt.xscale('log')
    plt.yscale('log')
    plt.grid(True)
    plt.legend()
    plt.show()

    dt_values = 0.01
    h_values = np.array([0.1, 0.05, 0.01, 0.005, 0.001])  # Spatial step sizes
    N_spatial_list = (1 / h_values).astype(int) - 1
    errors_l2 = np.zeros(len(N_spatial_list))
    errors_linf = np.zeros(len(N_spatial_list))
    for i, N in enumerate(N_spatial_list):
        x, y, U_numerical = solve_heat_equation(N, dt_values, T_end)
        U_exact = np.zeros_like(U_numerical)
        for j in range(len(x)):
            for k in range(len(y)):
                U_exact[j, k] = u_manufactured(x[j], y[k], T_end)
        error_l2 = norm(U_numerical - U_exact) / norm(U_exact)
        error_linf = np.max(np.abs(U_numerical - U_exact)) / np.max(np.abs(U_exact))
        errors_l2[-i], errors_linf[-i] = error_l2, error_linf

    plt.plot(h_values, errors_l2, '-o', label='L2 Error')
    plt.xlabel('Δh')
    plt.ylabel('Relative Error')
    plt.title('Temporal Convergence Study (L2 Error)')
    plt.xscale('log')
    plt.yscale('log')
    plt.grid(True)
    plt.legend()
    plt.show()

    plt.plot(h_values, errors_linf, '-o', label='Linf Error')
    plt.xlabel('Δh')
    plt.ylabel('Relative Error')
    plt.title('Temporal Convergence Study (Linf Error)')
    plt.xscale('log')
    plt.yscale('log')
    plt.grid(True)
    plt.legend()
    plt.show()