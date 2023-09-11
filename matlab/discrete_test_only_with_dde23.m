close all
clc

% Define the range of k and tau values
k_values = 0.1:0.1:50;
tau_values = 0.1:0.1:50;

Y = [  (mu + v) / beta;
     (mu * alpha * (beta - (mu + v))) / (beta * (mu + alpha * (mu + v)));
     (mu * (beta - (mu + v))) / (beta * (mu + alpha * (mu + v)))];

% Initialize stability matrix
stability_matrix = zeros(length(k_values), length(tau_values));

for i = 1:length(k_values)
    k = k_values(i);

    for j = 1:length(tau_values)
        tau = tau_values(j);

        % Define the initial conditions for S, I, and p
        initial_conditions = [0.1, 0.1, 0.1];
        
        % Define a function for delayed values
        Z_function = @(t) interp1(sol.x, sol.y.', t - tau, 'spline').';

        % Solve the DDE system using dde23
        lags = tau; % Use tau as the lag
        sol = dde23(@(t, Y, Z) dynamics_discrete(t, Y, Z_function, tau, k), lags, initial_conditions, [0, 50]);
        
        % Extract the final state
        final_state = sol.y(:, end);
        disp(final_state);
        
        % Calculate Jacobian matrix at the endemic steady state
        J = calculate_jacobian(final_state, tau, k);

        % Calculate eigenvalues of the Jacobian
        eig_values = eig(J);

        % Check stability criteria based on eigenvalues
        if all(real(eig_values) < 0)
            stability_matrix(i, j) = 1;  % Stable
        else
            stability_matrix(i, j) = 0;  % Unstable
        end
    end
end

% Create a stability plot with colors
figure;
imagesc(tau_values, k_values, stability_matrix');
colormap([1 0 0; 0 1 0]); % Red for unstable, green for stable
colorbar('Ticks', [0, 1], 'TickLabels', {'Unstable', 'Stable'});
xlabel('\tau');
ylabel('k');
title('Stability Analysis for E3 Steady State');

function dydt = dynamics_discrete(t, Y, Z_function, tau, k)
    mu = 3.9 * 10^(-5);
    v = 1/7;
    beta = 10 * (mu + v);
    alpha = 0.002;

    S = Y(1);
    I = Y(2);
    p = Y(3);
  
    % Evaluate the discrete distribution kernel function
    g_value = exp(-tau) * Z_function(t);
    
    dydt = [
        mu * (1 - p) - beta * S * I - mu * S;
        beta * S * I - (mu + v) * I;
        p * k * p * (1 - p) * (1 - alpha * trapz(g_value)) % trapz to integrate second exp. wrt. Z
    ];
end

function J = calculate_jacobian(final_state, tau, k)
    mu = 3.9 * 10^(-5);
    v = 1/7;
    beta = 10 * (mu + v);
    alpha = 0.002;

    S = (mu + v) / beta;
    I = (mu * alpha * (beta - (mu + v))) / (beta * (mu + alpha * (mu + v)));
    p = (mu * (beta - (mu + v))) / (beta * (mu + alpha * (mu + v)));

    g_value = exp(-tau); % At Z=0

    J = [
        -beta * I - mu, -beta * S, -mu;
        beta * I, -mu - v + beta * S, 0;
        0, k * p * (1 - p) * (1 - alpha * p * g_value), k * (1 - 2 * p) * (I - alpha * p - k * alpha * p * (1 - p) * g_value)
    ];
end
