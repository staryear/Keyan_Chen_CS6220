# This is a sample Python script.
import numpy as np
import time
import matplotlib.pyplot as plt
# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.


def runge_function(x):
    res = 1/(1+25*(x**2))
    return res

def equispaced_points(n):
    nodes = np.linspace(-1, 1, n)
    return nodes

def chebyshev_extrema(n):
    nodes = np.cos((2*np.arange(n) + 1) / (2*n) * np.pi)
    return nodes[::-1]

def invert_vandermonde_matrix(x):
    #construct inverted vandermonde matrix
    n = len(x)
    V = np.zeros((n, n))
    for i in range(n):
        V[:, i] = x**i

    invert_v = np.linalg.inv(V)
    return invert_v

def polynomial_interpolation(x):
    #interpolate runge function in to inverted V
    f = runge_function(x)
    V_inv = invert_vandermonde_matrix(x)
    res = V_inv.dot(f)
    return res

def condition_number(M):
    # calculate condition number
    M_inv = np.linalg.inv(M)
    con_number = max(np.sum(np.abs(M), axis=1)) * max(np.sum(np.abs(M_inv), axis=1))
    return con_number


def evaluate_polynomial_using_coef(coefficients, x):
    # evaluate_polynomial
    n = len(coefficients)
    y = np.zeros_like(x)
    for i in range(n):
        y += coefficients[i] * x**i
    return y

def l2_relative_error(y_true, y_pred):
    # calculate l2 error
    return np.linalg.norm(y_pred - y_true) / np.linalg.norm(y_true)


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    print(equispaced_points(8))
    print(chebyshev_extrema(8))
    print(invert_vandermonde_matrix(equispaced_points(8)))

    # max_nodes = 50
    # nodes_range = range(2, max_nodes + 1)
    # condition_numbers_equispaced = []
    # condition_numbers_chebyshev = []
    # # problem a
    # for n in nodes_range:
    #     # For equispaced nodes
    #     equispaced = equispaced_points(n)
    #     V_equispaced = invert_vandermonde_matrix(equispaced)
    #     condition_numbers_equispaced.append(condition_number(V_equispaced))
    #
    #     # For chebyshev nodes
    #     chebyshev = chebyshev_extrema(n)
    #     V_chebyshev = invert_vandermonde_matrix(chebyshev)
    #     condition_numbers_chebyshev.append(condition_number(V_chebyshev))
    #
    # plt.plot(nodes_range, condition_numbers_equispaced, label='Equispaced', marker='o', markersize=4)
    # plt.plot(nodes_range, condition_numbers_chebyshev, label='Chebyshev', marker='x', markersize=4)
    # plt.yscale('log')
    # plt.xlabel('Number of Nodes')
    # plt.ylabel('Condition Number')
    # plt.title('Condition Number vs Number of Nodes (Max Node: 50)')
    # plt.legend()
    # plt.grid(True)
    # plt.show()
    #
    #
    # #problem b
    #
    # plot_grid = np.linspace(-1, 1, 1000)
    # # when N = 9
    # equispaced = equispaced_points(9)
    # coeff_equispaced = polynomial_interpolation(equispaced)
    # y_equispaced = evaluate_polynomial_using_coef(coeff_equispaced, plot_grid)
    #
    # chebyshev = chebyshev_extrema(9)
    # coeff_chebyshev = polynomial_interpolation(chebyshev)
    # y_chebyshev = evaluate_polynomial_using_coef(coeff_chebyshev, plot_grid)
    #
    # plt.plot(plot_grid, y_equispaced, '--', label=f'Equispaced Points N=9', linewidth=1)
    # plt.plot(plot_grid, y_chebyshev, '-.', label=f'Chebyshev Points N=9', linewidth=1)
    # plt.title(f'Polynomial Interpolants for N=9')
    # plt.legend()
    # plt.xlabel('x')
    # plt.ylabel('y')
    # plt.show()
    #
    # #when N = 50
    # equispaced = equispaced_points(50)
    # coeff_equispaced = polynomial_interpolation(equispaced)
    # y_equispaced = evaluate_polynomial_using_coef(coeff_equispaced, plot_grid)
    #
    # chebyshev = chebyshev_extrema(50)
    # coeff_chebyshev = polynomial_interpolation(chebyshev)
    # y_chebyshev = evaluate_polynomial_using_coef(coeff_chebyshev, plot_grid)
    #
    # plt.plot(plot_grid, y_equispaced, '--', label=f'Equispaced Points N=50', linewidth=1)
    # plt.plot(plot_grid, y_chebyshev, '-.', label=f'Chebyshev Points N=50', linewidth=1)
    # plt.title(f'Polynomial Interpolants for N=50')
    # plt.legend()
    # plt.xlabel('x')
    # plt.ylabel('y')
    # plt.show()

    #problem c

    equispaced_nodes = np.linspace(-1, 1, 10000)
    y_true = runge_function(equispaced_nodes)
    errors = []
    times = []
    numbers = [i for i in range(3, 51)]
    for N in numbers:
        start_time = time.time()
        chebyshev = chebyshev_extrema(N)
        coeff_chebyshev = polynomial_interpolation(chebyshev)
        y_chebyshev = evaluate_polynomial_using_coef(coeff_chebyshev, equispaced_nodes)

        error = l2_relative_error(y_chebyshev, y_true)

        errors.append(error)
        end_time = time.time()
        total_time = end_time - start_time
        times.append(total_time)

    plt.plot(times, errors, marker='o')
    plt.xlabel('Total Time (seconds)')
    plt.ylabel('l2 errors')
    plt.title('Accuracy vs time')
    plt.show()

    plt.plot(numbers, errors, marker='o')
    plt.xlabel('nodes numbers')
    plt.ylabel('l2 errors')
    plt.title('Accuracy vs node')
    plt.show()

    plt.plot(numbers, times, marker='o')
    plt.xlabel('nodes numbers')
    plt.ylabel('Total Time (seconds)')
    plt.title('time vs nodes')
    plt.show()










# See PyCharm help at https://www.jetbrains.com/help/pycharm/
