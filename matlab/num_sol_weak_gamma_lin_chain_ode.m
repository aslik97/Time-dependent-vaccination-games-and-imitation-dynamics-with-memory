clear all;
close all;

% Define vector with all delays
lags = [20, 50, 100, 150];

for j = 1:length(lags)
    tau = lags(j);
    % Initial conditions
    mu = 3.9 * 10^(-5);
    v = 1/7;
    R0 = 10;
    S0 = 1/R0;
    I0 = mu * (1 - 1/R0) / (mu + v);
    p0 = 0.05;
    M0 = p0; % Initial value of M can be set to the initial value of p
    
    y0 = [S0; I0; p0; M0];
    
    [t, y] = ode45(@(t,y) rhs_ode(t, y,tau), [0, 4000], y0);
    
    % Plot the results
    figure;
    plot(t, y(:,2), 'r-', 'Linewidth', 2);
    title(['\tau = ', num2str(tau)], 'FontSize', 16);
    xlabel('time, days', 'FontSize', 14);
    ylabel('I(t)', 'FontSize', 14);
    ylim([0 3e-4]);
end

function dydt = rhs_ode(t, y,tau)
    S = y(1);
    I = y(2);
    p = y(3);
    M = y(4);
    
    k = 400;
    mu = 3.9 * 10^(-5);
    v = 1/7;
    beta = 10 * (mu + v);
    alpha = 0.002;
    sigma = 1/tau; % passed to the function
    
    dSdt = mu * (1 - p) - beta * S * I - mu * S;
    dIdt = beta * S * I - (mu + v) * I;
    dpdt = k * p * (1 - p) * (I - alpha * M);
    dMdt = sigma * p - sigma * M;
    
    dydt = [dSdt; dIdt; dpdt; dMdt];
end
