% Load the parameters and function from dynamics_function.m
mu = 3.9*exp(-5);
v = 1/7;
beta = 10 * (mu + v);
alpha = 0.002;

% Define the endemic steady state
E3 = [ (mu+v)/beta ;
      (mu*alpha*(beta - (mu+v)))/(beta*(mu + alpha*(mu+v))) ;
      (mu*(beta - (mu+v)))/(beta*(mu + alpha*(mu+v))) ];

% Define the gamma distribution kernel function
%n_weak = 1;
%sigma = 2;
%g = @(s) (s.^(n_weak - 1) .* sigma^n_weak .* exp(-sigma * s)) / factorial(n_weak - 1);

% Define values of tau and k for analysis
tau_values = 0.1:0.1:10; % Mean time delay values
k_values = 0.1:0.1:10;   % Imitation rate values

% Initialize stability matrix (1: stable, 0: unstable)
stability_matrix = zeros(length(tau_values), length(k_values));

% Loop through different tau and k values
for i = 1:length(tau_values)
    tau = tau_values(i);
    for j = 1:length(k_values)
        k = k_values(j);
        
        % Define the initial conditions for S, I, and p
        initial_conditions = [0.1, 0.1, 0.1];
        
        % Solve the DDE system using dde23
        lags = tau; % Use tau as the lag
        sol = dde23(@(t, Y, Z) dynamics(t, Y, Z, tau,k), lags, initial_conditions, [0, 100]);
        
        % Extract the final state
        final_state = sol.y(:, end);
        
        % Calculate Jacobian matrix at the endemic steady state
        J = [
            -beta * final_state(2) - mu, -beta * final_state(1), 0;
            beta * final_state(2), -mu - v, 0;
            0, 0, k * (1 - 2 * final_state(3)) * (1 - 3 * final_state(3))
        ];
        
        % Calculate eigenvalues of the Jacobian
        eig_values = eig(J);
        
        % Determine stability
        if all(real(eig_values) < 0)
            stability_matrix(i, j) = 1; % Stable
        end
    end
end

% Plot the regions of stability and instability
figure;
imagesc(tau_values, k_values, stability_matrix);
colormap([1 0 0; 0 1 0]);
colorbar;
xlabel('\tau (Mean Time Delay)');
ylabel('k (Imitation Rate)');
title('Regions of Stability (Green) and Instability (Red)');

