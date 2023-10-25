% Define parameters
mu = 3.9 * 10^(-5); % Correcting the value of mu
v = 1/7;
beta = 10 * (mu + v);
alpha = 0.002;

% Define the range of k and tau values
%k_values = linspace(0, 400, 100); % Adjusted the range based on your graph
%tau_values = linspace(0.1, 200, 100); 
k_values = 0.1:10:50;  % Jump by 10 instead of 0.1
tau_values = 0.1:10:50;  % Similarly, jump by 10


% Initialize stability matrix
stability_matrix = zeros(length(k_values), length(tau_values));

% Loop through k and tau values
for i = 1:length(k_values)
    for j = 1:length(tau_values)
        k = k_values(i);
        tau = tau_values(j);

        % Solve the DDE system using dde23
        lags = tau; % Use tau as the lag
        sol = dde23(@(t, Y, Z) dynamics_discrete(t, Y, Z, k), lags, [0.8, 0.1, 0.1], [0, 50]); 
        
        % Calculate Jacobian matrix at the endemic steady state
        Y = sol.y(:, end);
        J = jacobian_at(Y, beta, mu, v, k, alpha, tau);
        
        % Calculate eigenvalues of the Jacobian
        eig_values = eig(J);
        
        if all(real(eig_values) < 0)
            stability_matrix(i, j) = 1; % Stable
        end
    end
end

% Create a stability plot
figure;
imagesc(tau_values, k_values, stability_matrix);
colorbar;
xlabel('\tau');
ylabel('k');
title('Stability Analysis for E3 Steady State');

function dy = dynamics_discrete(t, Y, Z, k)
    mu = 3.9 * 10^5;
    v = 1/7;
    beta = 10 * (mu + v);
    alpha = 0.002;

    S = Y(1);
    I = Y(2);
    p = Y(3);

    % Using the delayed value directly
    p_delayed = Z(3);

    dS = mu * (1 - p) - beta * S * I - mu * S;
    dI = beta * S * I - (mu + v) * I;
    dp = k * p * (1 - p) * (1 - alpha * p_delayed);
    dy = [dS; dI; dp];
end

function J = jacobian_at(Y, beta, mu, v, k, alpha, tau)
    S = Y(1);
    I = Y(2);
    p = Y(3);

    J = zeros(3,3);
    J(1,:) = [-beta*I-mu, -beta*S, 0];
    J(2,:) = [beta*I, beta*S-mu-v, 0];
    J(3,:) = [0, k*p*(1-p), k*(1-2*p)*(1-alpha*p)]; 
    % simplified Jacobian matrix. 
end


