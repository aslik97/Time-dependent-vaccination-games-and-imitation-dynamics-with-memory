
% Define parameters
gamma = 1.41e-4;
sigma = 0.9; % Different sigma values for comparison
mu = 1/(11*365);          % Birth and death rate (1/L, where L is life expectancy)
v = 1/7;           % Rate of recovery from infection
R_0 = 15;          % Basic Reproduction Number
beta = R_0*(mu+v); % Transmission rate
p_3 = 0.90;             %critical elimination coverage
p_c =0.933; 
gamma_hat = 0;                    %for the I-model
alpha_hat = 0.091e-4;       %defined from me
k_theta = 40;                  %from the condition k*theta = 40 and taking theta(I) = theta_A*I

% Initial conditions
S0 = 1.04/R_0;
I0 = 0.82e-5;
p0 = 0.95;

% Time span
tspan = [0, 365*11];

% Simulation and Plotting
    [T, Z] = ode45(@(t, y) seasonalityDynamics(t, y, mu, alpha_hat, k_theta, gamma, sigma,beta,v), tspan, [S0 I0 p0]);
    [K,M] = ode45(@(t, y) nonseasonalityDynamics(t, y, mu, alpha_hat, k_theta, gamma, sigma,beta,v), tspan, [S0 I0 p0]);
    
% Plots
figure;
subplot(1,2,1);
plot(T/365, R_0*Z(:,1)); % R_E(t) = R_0*S(t)
title('R_E(t) vs t');
xlabel('time');
ylabel('R_E(t)');
hold on;
plot(K/365, R_0*M(:,1)); % R_E(t) = R_0*S(t)
legend('with seasonality', 'without seasonality');

subplot(1,2,2);
plot(T/365, Z(:,2));
title('I(t) vs t');
xlabel('time');
ylabel('I(t)');
hold on;
plot(K/365, R_0*M(:,2)); % R_E(t) = R_0*S(t)
legend('with seasonality', 'without seasonality');
%subplot(1,3,3);
%plot(T/365, Y(:,3));
%title('p(t) vs t');
%xlabel('time');
%ylabel('p(t)');
%hold on;
%line([0 365*11],[p_c p_c],'Color','red','LineStyle','--');
%legend('p(t)', 'p_c');

% Differential equations function
function dydt = seasonalityDynamics(t, y, mu, alpha_hat, k_theta, gamma_hat, sigma, beta, v)
    S = y(1);
    I = y(2);
    p = y(3);
    p_3 = 0.90;    
    % Define c(t) with sigma
    epsilon= 0.1;
    T = 365;        % The period of the function, one year
    C = 0;          % Phase shift (in days)
    c_t= epsilon * sin(2*pi*(t - C)/T) + 1;
    % c_t = (1 - sigma) * sin(2*pi*t/365) + 1 + sigma;

    % System of equations
    dSdt = mu * (1 - p_3) - c_t * beta * S * I - mu * S;
    dIdt = c_t * beta * S * I - (mu + v) * I;
    dpdt = k_theta*(1-p_3) * ((I - alpha_hat*p_3)*p_3 + gamma_hat);
    
    dydt = [dSdt; dIdt; dpdt];
end

function dydt = nonseasonalityDynamics(t, y, mu, alpha_hat, k_theta, gamma_hat, sigma, beta, v)
    S = y(1);
    I = y(2);
    p = y(3);
    p_3 = 0.90;    
    % Define c(t) with sigma
    
    T = 365;        % The period of the function, one year
    C = 0;          % Phase shift (in days)
    
    % System of equations
    dSdt = mu * (1 - p_3) - beta * S * I - mu * S;
    dIdt = beta * S * I - (mu + v) * I;
    dpdt = k_theta*(1-p_3) * ((I - alpha_hat*p_3)*p_3 + gamma_hat);
    
    dydt = [dSdt; dIdt; dpdt];
end
