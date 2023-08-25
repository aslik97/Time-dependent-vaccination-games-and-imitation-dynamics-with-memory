"""Stability regions for E3 with weak kernel distribution."""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# Load the parameters
mu = 3.9 * np.exp(-5)
v = 1 / 7
beta = 10 * (mu + v)
alpha = 0.002

# Define the endemic steady state
E3 = [
    (mu + v) / beta,
    (mu * alpha * (beta - (mu + v))) / (beta * (mu + alpha * (mu + v))),
    (mu * (beta - (mu + v))) / (beta * (mu + alpha * (mu + v))),
]

# Define values of tau and k for analysis
tau_values = np.arange(0.1, 10.1, 0.1)  # Mean time delay values
k_values = np.arange(0.1, 10.1, 0.1)    # Imitation rate values

# Initialize stability matrix (1: stable, 0: unstable)
stability_matrix = []

# Define the dynamics function
def dynamics(t, Y, tau, k):
    S, I, p = Y
    n_weak = 1
    sigma = 2
    g = lambda s: (s**(n_weak - 1) * sigma**n_weak * np.exp(-sigma * s)) / np.math.factorial(n_weak - 1)
    
    # Evaluate the gamma distribution kernel function
    g_value = g(t - tau)
    
    dydt = [
        mu * (1 - p) - beta * S * I - mu * S,
        beta * S * I - (mu + v) * I,
        p * k * p * (1 - p) * (1 - alpha * np.trapz(p * g_value, dx=0.1))
    ]
    
    # Flatten the dydt array
    return np.array(dydt).flatten()

# Loop through different tau and k values
for i, tau in enumerate(tau_values):
    stability_matrix_row = []  # Initialize a row for the stability matrix
    for j, k in enumerate(k_values):
        
        # Define the initial conditions for S, I, and p
        initial_conditions = [0.1, 0.1, 0.1]
        
        # Define the time span
        t_span = [0, 100]  # Time span from 0 to 100
        
        # Solve the DDE system using solve_ivp without specifying t_eval
        sol = solve_ivp(
            lambda t, Y: dynamics(t, Y, tau, k),
            t_span,
            initial_conditions
        )
        
        # Extract the final state
        final_state = sol.y[:, -1]
        
        # Calculate Jacobian matrix at the endemic steady state
        J = np.array([
            [-beta * final_state[2] - mu, -beta * final_state[0], 0],
            [beta * final_state[2], -mu - v, 0],
            [0, 0, k * (1 - 2 * final_state[2]) * (1 - 3 * final_state[2])]
        ])
        
        # Calculate eigenvalues of the Jacobian
        eig_values = np.linalg.eigvals(J)
        
        # Determine stability and append to the row
        if all(np.real(eig_values) < 0):
            stability_matrix_row.append(1)  # Stable
        else:
            stability_matrix_row.append(0)  # Unstable
            
    # Append the row to the stability matrix
    stability_matrix.append(stability_matrix_row)

# Convert the list of lists to a numpy array for easier visualization
stability_matrix = np.array(stability_matrix)

# Plot the regions of stability and instability
plt.figure()
plt.imshow(stability_matrix, extent=[0.1, 10.0, 0.1, 10.0], cmap='RdYlGn', origin='lower', aspect='auto')
plt.colorbar(label='Stability')
plt.xlabel('τ (Mean Time Delay)')
plt.ylabel('k (Imitation Rate)')
plt.title('Regions of Stability (Green) and Instability (Red)')
plt.show()
