

% Define parameters
mu = 3.9 * 10^(-5);
v = 1/7;
beta = 10 * (mu + v);
alpha = 0.002;

% Define the range of k and tau values
k_values = 0.1:0.1:500;
tau_values = 0.1:0.1:300;

% Initialize stability matrix
stability_matrix = zeros(length(k_values), length(tau_values));

for i = 1:length(k_values)
    k = k_values(i);

    for j = 1:length(tau_values)
        tau = tau_values(j);

        % steady state
        S_star = (mu + v) / beta;
        I_star = (mu * alpha * (beta - (mu + v))) / (beta * (mu + alpha * (mu + v)));
        p_star = (mu * (beta - (mu + v))) / (beta * (mu + alpha * (mu + v)));
        
        % Calculate coefficients
        A = 1;
        B = (beta * I_star + mu)^2 - 2 * beta^2 * S_star * I_star - (k * alpha * p_star * (1 - p_star))^2;
        C = (beta^2 * S_star * I_star)^2 - 2 * k * mu * beta * p_star * (1 - p_star) * I_star * (beta * I_star + mu) ...
            - (k * alpha * p_star * (1 - p_star))^2 * ((beta * I_star + mu)^2 - 2 * beta^2 * S_star * I_star);
        D = (k * beta * p_star * (1 - p_star) * I_star)^2 * (mu^2 - alpha^2 * beta^2 * S_star^2);

        % Calculate roots of the cubic equation
        roots_Cubic = roots([A, 0, B, 0, C, 0, D]);