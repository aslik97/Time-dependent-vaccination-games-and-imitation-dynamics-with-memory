function dydt = dynamics(t, y, mu, v, beta, k_theta, alpha_hat, gamma_hat,R_0)

    S = y(1);
    I = y(2);
    p = y(3);

    % Cp_imitation_3
    p_imitation_3 = 0.90;%(mu * p_c) / (mu + (mu + v) * alpha_hat);

    % Immigration rate (1 individual per week in a population of 5e6)
    Imm_rate = 1 / (5e6 *7);

    % Differential equations
    % Define the system of differential equations
    dSdt = mu*(1-p_imitation_3 ) - mu*S - beta*S*I; % dS/dt
    dIdt = beta*S*I - (mu+v)*I;        % dI/dt
    dpdt = k_theta*(1-p_imitation_3 ) * ((I - alpha_hat*p)*p_imitation_3  + gamma_hat);%1-(1/R_0)-(mu + v)/(mu)*I; %k_theta*(1-p) * ((I - alpha_hat*p)*p + gamma_hat); % dp/dt

    %column vector
    dydt = [dSdt; dIdt; dpdt];
end
