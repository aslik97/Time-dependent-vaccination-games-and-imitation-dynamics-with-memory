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
tau_values = 0.1:0.1:1; % Mean time delay values
k_values = 0.1:0.1:1;   % Imitation rate values

% Initialize stability matrix (1: stable, 0: unstable)
stability_matrix = zeros(length(tau_values), length(k_values));

 % Define the Laplace transform of the kernel
  %      Lg = @(lambda) sigma / (lambda + sigma)^p_star;


% Loop through different tau and k values
for i = 1:length(tau_values)
    tau = tau_values(i);
    for j = 1:length(k_values)
        k = k_values(j);
        
        % Define the initial conditions for S, I, and p
        initial_conditions = [0.1, 0.1, 0.1];
        
        % Solve the DDE system using dde23
        lags = tau; % Use tau as the lag
        sol = dde23(@(t, Y, Z) dynamics_weak_gamma(t, Y, Z, tau,k), lags, initial_conditions, [0, 1]);
        
        % Extract the final state
        final_state = sol.y(:, end);

         % Calculate coefficients
        A = 1; %TODO
        B = ;
        C = ;
        D = ;
        E = ;

        % Calculate roots of the cubic equation
        roots_Cubic = roots([A, 0, B, 0, C, 0, D]);
        
        
        % Calculate Jacobian matrix at the endemic steady state
       % J = [
        %    -beta * final_state(2) - mu, -beta * final_state(1), -mu;
         %   -beta * final_state(2), -mu - v+beta*final_state(1), 0;
         %   0, k*final_state(3)*(1-final_state(3)), k * (1 - 2 * final_state(3)) * (final_state(2) - alpha * final_state(3)-k*alpha*final_state(3)*(1-final_state(3))*Lg)
        %];
        
        % Calculate eigenvalues of the Jacobian
        %eig_values = eig(J);
        
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

