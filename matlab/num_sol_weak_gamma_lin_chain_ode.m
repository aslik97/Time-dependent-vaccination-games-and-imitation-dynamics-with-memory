clear all;
close all;

% Define vector with all delays
lags = [20, 50, 100, 150];

 % Initial conditions
    mu = 3.9e-5;
    v = 1/7;
    R0 = 10;
    S0 = 1/R0;
    I0 = mu * (1 - 1/R0) / (mu + v);
    alpha= 0.002;
    k=400;
    p0 = 0;
    beta = R0 * (mu + v);
for j = 1:length(lags)
    tau = lags(j);
    sigma = 1 / tau;   
    integrand = @(s) sigma *tau* exp(-sigma * s);  %integrand function
    M0 = integral(integrand, 0, Inf);  % integral from 0 to infinity
    disp(M0);
   %M0 =  tau; %int_0^infty (sigma*exp(-sigma*s)*p(-s)) ds; % Initial value
   %of M can be set to the initial value of p, because it converges to p0
    [t, y] = ode45(@(t,y) rhs_ode(t, y,tau,k,mu,v,beta,alpha,sigma), [0, 4000], [S0 I0 p0 M0]);
    
    % Plot the results
    figure;
    plot(t, y(:,2), 'r-', 'Linewidth', 2);
    title(['\tau = ', num2str(tau)], 'FontSize', 16);
    xlabel('time, days', 'FontSize', 14);
    ylabel('I(t)', 'FontSize', 14);
    ylim([0 3e-4]);
end

function dydt = rhs_ode(t, y,tau,k,mu,v,beta,alpha,sigma)
    S = y(1);
    I = y(2);
    p = y(3);
    M = y(4);
    
   
    dSdt = mu * (1 - p) - beta * S * I - mu * S;
    dIdt = beta * S * I - (mu + v) * I;
    dpdt = k * p * (1 - p) * (I - alpha * M);
    dMdt = sigma * p - sigma * M;
    
    dydt = [dSdt; dIdt; dpdt; dMdt];
end
