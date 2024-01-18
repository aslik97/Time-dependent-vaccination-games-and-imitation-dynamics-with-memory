close all
clear all
% Define parameters
gamma = 1.41e-4;

mu = 1/(11*365);   % Birth and death rate (1/L, where L is life expectancy)
v = 1/7;           % Rate of recovery from infection
R_0 = 15;          % Basic Reproduction Number
beta = R_0*(mu+v); % Transmission rate
p_3 = 0.90;             %critical elimination coverage
p_c =0.933; 
gamma_hat = 0;                    %for the I-model
alpha_hat = 0.091e-4;       %defined from me
k_theta = 40;                  %from the condition k*theta = 40 and taking theta(I) = theta_A*I
epsilon_values = [ 0, 0.01, 0.1, 0.5];

% Initial conditions
S0 = 1.04/R_0;
I0 = 0.82e-5;
p0 = 0.95;

% Time span
tspan = [0, 365*11];
for i = 1:length(epsilon_values)
    epsilon = epsilon_values(i);
% Simulation and Plotting
    [T, Z] = ode45(@(t, y) seasonalityDynamics(t, y, mu, alpha_hat, k_theta, gamma,beta,v,epsilon), tspan, [S0 I0 p0]);
    [K,M] = ode45(@(k, m) nonseasonalityDynamics(k, m, mu, alpha_hat, k_theta, gamma,beta,v), tspan, [S0 I0 p0]);

    if epsilon == 0
        disp(['Max difference in R_E(t) when epsilon is ', num2str(epsilon), ': ', num2str(max(abs(R_0*Z(:,1) - R_0*M(:,1))))]);
        disp(['Max difference in I(t) when epsilon is ', num2str(epsilon), ': ', num2str(max(abs(Z(:,2) - M(:,2))))]);
    end
     % Subplot for R_E(t)
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
% Seasonality Dynamics
function dydt = seasonalityDynamics(t, y, mu, alpha_hat, k_theta, gamma_hat, beta, v, epsilon)
    S = y(1);
    I = y(2);
    p = y(3);
    p_3 = 0.90;
     T = 365;
 % The period of the function, one year
 C = 0; 

    c_t = epsilon * sin((2*pi/T)*(t -C)) + 1; % Should be 1 when epsilon = 0
    if(epsilon == 0)
    disp(c_t);
    end
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

