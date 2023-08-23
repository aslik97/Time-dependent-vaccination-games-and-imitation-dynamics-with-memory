import numpy as np
import matplotlib.pyplot as plt
from scipy.special import factorial

# Define parameters
n_weak = 1
n_strong = 2
sigma_weak = 2
sigma_strong = 2
tau = 2
x = np.arange(0, 5, 0.01)  # Range of x values

# Calculate gamma distribution kernels
weak_kernel = (x**(n_weak - 1) * sigma_weak**n_weak * np.exp(-sigma_weak * x)) / factorial(n_weak - 1)
strong_kernel = (x**(n_strong - 1) * sigma_strong**n_strong * np.exp(-sigma_strong * x)) / factorial(n_strong - 1)

# Define parameters for AF distribution kernel
t = tau  # The parameter t
t1_values = np.array([0.1])  # Different t1 values

# Calculate AF distribution kernels
AF_kernels = np.zeros((len(t1_values), len(x)))
for i, t1 in enumerate(t1_values):
    t2 = t - t1
    AF_kernels[i, :] = (np.exp(-x / t1) - np.exp(-x / t2)) / (t1 - t2)

# Create the plot
plt.figure()
plt.plot(x, weak_kernel, 'b', linewidth=2, label='Weak Gamma (n=1)')
plt.plot(x, strong_kernel, 'r', linewidth=2, label='Strong Gamma (n=2)')
for i, t1 in enumerate(t1_values):
    plt.plot(x, AF_kernels[i, :], color=[0.2, 0.7, 0.2], label=f'AF (t1 = {t1})')
plt.xlabel('x')
plt.ylabel('Kernel Value')
plt.title('Weak, Strong, and AF Distribution Kernels')
plt.legend()

# Set axis limits
plt.xlim(0, 5)
plt.ylim(0, 2)

# Show grid
plt.grid(True)

# Show the plot
plt.show()
