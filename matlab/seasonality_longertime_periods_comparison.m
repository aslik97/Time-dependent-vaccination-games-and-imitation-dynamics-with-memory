% Extended time span
tspan = [0, 365*50]; % For example, 20 years

% Initialize arrays to store total infections
total_infections_seasonal = zeros(length(epsilon_values), 1);
total_infections_nonseasonal = zeros(length(epsilon_values), 1);

for i = 1:length(epsilon_values)
    epsilon = epsilon_values(i);
  
 % Simulation for the seasonal model
    [T, Z] = ode45(@(t, y) seasonalityDynamics(t, y, mu, alpha_hat, k_theta, gamma,beta,v,epsilon), tspan, [S0 I0 p0]);

    % Simulation for the non-seasonal model
    [K, M] = ode45(@(k, m) nonseasonalityDynamics(k, m, mu, alpha_hat, k_theta, gamma,beta,v), tspan, [S0 I0 p0]);

    % Calculate total infections for seasonal model
    total_infections_seasonal(i) = trapz(T, Z(:,2));

    % Calculate total infections for non-seasonal model
    total_infections_nonseasonal(i) = trapz(K, M(:,2));

    subplot(4, 2, i * 2 - 1);
    plot(T/365, R_0*Z(:,1));
    title(['R_E(t) vs t for \epsilon = ', num2str(epsilon)]);
    xlabel('Time (years)');
    ylabel('R_E(t)');
    hold on;
    plot(K/365, R_0*M(:,1)); % R_E(t) = R_0*S(t)
    legend('with seasonality', 'without seasonality');

    % Subplot for I(t)
    subplot(4, 2, i * 2);
    plot(T/365, Z(:,2));
    title(['I(t) vs t for \epsilon = ', num2str(epsilon)]);
    xlabel('Time (years)');
    ylabel('I(t)');
    hold on;
    plot(K/365, R_0*M(:,2)); % R_E(t) = R_0*S(t)
    legend('with seasonality', 'without seasonality');
end

% Display the results
for i = 1:length(epsilon_values)
    disp(['Total infections with seasonality (epsilon = ', num2str(epsilon_values(i)), '): ', num2str(total_infections_seasonal(i))]);
    disp(['Total infections without seasonality (epsilon = ', num2str(epsilon_values(i)), '): ', num2str(total_infections_nonseasonal(i))]);
end

function dydt = seasonalityDynamics(t, y, mu, alpha_hat, k_theta, gamma_hat, beta, v, epsilon)
    S = y(1);
    I = y(2);
    p = y(3);
    p_3 = 0.90;
     T = 365;
 % The period of the function, one year
 C = 0; 

    c_t = epsilon * sin((2*pi/T)*(t -C)) + 1; % Should be 1 when epsilon = 0
    
    dSdt = mu * (1 - p_3) - c_t * beta * S * I - mu * S;
    dIdt = c_t * beta * S * I - (mu + v) * I;
    dpdt = k_theta * (1 - p_3) * ((I - alpha_hat * p_3) * p_3 + gamma_hat);
    
    dydt = [dSdt; dIdt; dpdt];
end

% Non-Seasonality Dynamics
function dmdt = nonseasonalityDynamics(k, m, mu, alpha_hat, k_theta, gamma_hat, beta, v)
    S = m(1);
    I = m(2);
    p = m(3);
    p_3=0.90;
    dSdt = mu * (1 - p_3) - beta * S * I - mu * S;
    dIdt = beta * S * I - (mu + v) * I;
    dpdt = k_theta * (1 - p_3) * ((I - alpha_hat * p_3) * p_3 + gamma_hat);
    
    dmdt = [dSdt; dIdt; dpdt];
end

