% Define parameters
mu = 3.9 * 10^(-5);
v = 1/7;
beta = 10 * (mu + v);
alpha = 0.002;

% Define the range of k and tau values
k_values = 0.1:0.1:50;
tau_values = 0.1:0.1:50;

% Initialize stability matrix
stability_matrix = zeros(length(k_values), length(tau_values));

% Define the symbolic variable lambda
syms lambda

for i = 1:length(k_values)
    k = k_values(i);

    for j = 1:length(tau_values)
        tau = tau_values(j);

        % steady state
       Y = [  (mu + v) / beta;
     (mu * alpha * (beta - (mu + v))) / (beta * (mu + alpha * (mu + v)));
     (mu * (beta - (mu + v))) / (beta * (mu + alpha * (mu + v)))];

        % Define a function for delayed values
        Z_function = @(t) interp1(sol.x, sol.y.', t - tau, 'spline').';
        
        % Define the initial conditions for S, I, and p
        initial_conditions = [0.1, 0.1, 0.1];
        
        Lg = @(s) exp(-s * tau);
        % Solve the DDE system using dde23
        lags = tau; % Use tau as the lag
        sol = dde23(@(t, Y, Z) dynamics_discrete(t, Y, Z, tau,k), lags, initial_conditions, [0, 1]);
        
        % Extract the final state
        final_state = sol.y(:, end);
        disp(final_state);

         % Calculate Jacobian matrix at the endemic steady state
        J = [
            -beta * final_state(2) - mu - lambda, -beta * final_state(1), -mu;
            -beta * final_state(2), -mu - v+beta*final_state(1)- lambda, 0;
            0, k*final_state(3)*(1-final_state(3)), k * (1 - 2 * final_state(3)) * (final_state(2) - alpha * final_state(3)-k*alpha*final_state(3)*(1-final_state(3))*Lg(lambda)-lambda)
        ];
        
        % Calculate eigenvalues of the Jacobian
        eig_values = eig(J);
        
         if all(real(eig_values) < 0)
            stability_matrix(i, j) = 1; % Stable
         end
        % Calculate coefficients
        %A = 1;
        %B = (beta * final_state(2) + mu)^2 - 2 * beta^2 * final_state(1) * final_state(2) - (k * alpha * final_state(3) * (1 - final_state(3)))^2;
        %C = (beta^2 * final_state(1) * final_state(2))^2 - 2 * k * mu * beta * final_state(3) * (1 - final_state(3) * final_state(2)) * (beta * final_state(2) + mu) - (k * alpha * p_star * (1 - final_state(3)))^2 * ((beta * final_state(2) + mu)^2 - 2 * beta^2 * final_state(1) * final_state(2));
        %D = (k * beta * final_state(3) * (1 - final_state(3)) * final_state(2))^2 * (mu^2 - alpha^2 * beta^2 * final_state(1)^2);
        %coefficients = [A,B,C,D];
        % Calculate roots of the cubic equation
        %roots_Cubic = roots(coefficients );
        
        % Extract real parts of roots
        %max_real_parts = max(real(roots_Cubic));
        
        %disp(max_real_parts);
        % Check stability criteria based on the real parts of roots
        %if (max_real_parts > 0)
         %  stability_matrix(i, j) = 0;  % Unstable
        %else
         %  stability_matrix(i, j) = 1;  % Stable
        %end
    end
end


%disp(stability_matrix);
% Create a stability plot
figure;
imagesc(tau_values, k_values, stability_matrix');
colormap([1 0 0; 0 1 0]); % Red for unstable, green for stable
colorbar('Ticks', [0, 1], 'TickLabels', {'Unstable', 'Stable'});
xlabel('\tau');
ylabel('k');
title('Stability Analysis for E3 Steady State');

function dydt = dynamics_discrete(t, Y, Z, tau,k)
    
    mu = 3.9 * 10^(-5);
    v = 1/7;
    beta = 10 * (mu + v);
    alpha = 0.002;
    S = Y(1);
    I = Y(2);
    p = Y(3);
  
    g = @(s) exp(-s * tau);
    % Evaluate the discrete distribution kernel function
    g_value = g(t - tau - Z);
    
    dydt = [
        mu * (1 - p) - beta * S * I - mu * S;
        beta * S * I - (mu + v) * I;
        p * k * p * (1 - p) * (1 - alpha * trapz(Z, p .* g_value)) % trapz to integrate second exp. wrt. Z
    ];
end

