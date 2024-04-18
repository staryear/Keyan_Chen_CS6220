import numpy as np
import matplotlib.pyplot as plt

nu_val = 1.0
dt_val = 0.01
h_val = 0.1

# Define the function show stability using the computed expression
def G_modulus(k, l):
    sigma_val = (nu_val * dt_val / h_val**2)
    sin_k = np.sin(h_val * k / 2)**2
    sin_l = np.sin(h_val * l / 2)**2
    numerator = 4 * sigma_val * sin_k + 4 * sigma_val * sin_l - 1
    denominator = 4 * sigma_val * sin_k + 4 * sigma_val * sin_l + 1
    return np.abs(numerator / denominator)

# Create a grid of wave number values
if __name__ == '__main__':
    k_values = np.linspace(-np.pi/h_val, np.pi/h_val, 40)
    l_values = np.linspace(-np.pi/h_val, np.pi/h_val, 40)
    K, L = np.meshgrid(k_values, l_values)

    # Compute xi modulus over the grid
    G_mod = G_modulus(K, L)

    contour = plt.contourf(K, L, G_mod, levels=np.linspace(0, 1, 50), cmap='viridis')
    cbar = plt.colorbar(contour)
    cbar.set_label('Factor G')
    plt.title('Factor G for Crank-Nicolson Scheme')
    plt.xlabel('k')
    plt.ylabel('p')
    plt.grid(True)
    plt.show()