clear all
close all
clc

mu = 3.9 * 10^(-5);
v = 1/7;
beta = 10 * (mu + v);
alpha = 0.002;

% Define the range of k and tau values
k_values = 0.1:0.1:100;
tau_values = 0.1:0.1:100;

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
        
        % Calculate the result of the integral
        %result = calculate_integral(tau, k);

        % Solve the DDE system using dde23
        lags = tau; % Use tau as the lag

        sol = dde23(@(t, Y, Z) dynamics_discrete(t, Y, Z, tau, k), lags, initial_conditions, [0, 10]);
        disp(size(sol.x))
        % Extract the final state
        final_state = sol.y(:, end);
        
        % Calculate Jacobian matrix at the endemic steady state
        J = calculate_jacobian(final_state, tau, k);

        % Calculate eigenvalues of the Jacobian
        eig_values = eig(J);
        disp(eig_values);
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

function dydt = dynamics_discrete(t, Y, Z, tau, k)
    mu = 3.9 * 10^(-5);
    v = 1/7;
    beta = 10 * (mu + v);
    alpha = 0.002;

    S = Y(1);
    I = Y(2);
    p = Y(3);

     % Evaluate the discrete distribution kernel function
    sigma = 0.01; % Adjust the width as needed
    result = discrete_kernel(tau, 0.1, sigma);
   
    dydt = [
        mu * (1 - p) - beta * S * I - mu * S;
        beta * S * I - (mu + v) * I;
        p * k * p * (1 - p) * (1 - alpha * trapz(result)) % trapz to integrate second exp. wrt. Z
    ];
end

function J = calculate_jacobian(final_state, tau, k)
    mu = 3.9 * 10^(-5);
    v = 1/7;
    beta = 10 * (mu + v);
    alpha = 0.002;

    S = final_state(1);
    I = final_state(2);
    p = final_state(3);

    g_value = exp(-tau); % At Z=0

    J = [
        -beta * I - mu, -beta * S, -mu;
        beta * I, -mu - v + beta * S, 0;
        0, k * p * (1 - p) * (1 - alpha * p * g_value), k * (1 - 2 * p) * (I - alpha * p - k * alpha * p * (1 - p) * g_value)
    ];
end

function g_value = discrete_kernel(t, tau_0, sigma)
    % t: The time variable
    % tau_0: The center of the Dirac delta function
    % sigma: The width of the Gaussian approximation
    
    g_value = exp(-(t - tau_0).^2 / (2 * sigma^2)) / (sigma * sqrt(2 * pi));
end

