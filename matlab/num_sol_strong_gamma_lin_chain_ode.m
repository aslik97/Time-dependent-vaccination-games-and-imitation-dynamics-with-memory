clear all
close all
clc
   mu = 3.9e-5;
    v = 1/7;
    R0 = 10;
    S0 = 1/R0;
    I0 = mu * (1 - 1/R0) / (mu + v);
    alpha= 0.002;
    k=400;
    p0 = 0.95;
    beta = R0 * (mu + v);
tau_vector=[20,50,100,150]; %constant delay 

for i = 1:length(tau_vector)
     tau = tau_vector(i);
    sigma = 1 / tau;   
    integrand = @(s) sigma * exp(-sigma * s)*tau;  %integrand function
    M10 = integral(integrand, 0, Inf);  % integral from 0 to infinity
    disp(M10);
    M20=integral(integrand, 0, Inf);  % integral from 0 to infinity
   %M0 =  tau; %int_0^infty (sigma*exp(-sigma*s)*p(-s)) ds; % Initial value
   %of M can be set to the initial value of p, because it converges to p0
   tspan = [0, 4000];
   [t, y] = ode45(@(t,y) rhs_ode(t, y,tau,k,mu,v,beta,alpha,sigma), tspan, [S0 I0 p0 M10 M20]);
    
    figure;
    plot(t, y(:,2), 'r-', 'Linewidth', 2);
    title(['\tau = ', num2str(tau)], 'FontSize', 16);
    xlabel('time, days', 'FontSize', 14);
    ylabel('I(t)', 'FontSize', 14);
    ylim([0 3e-4]);
end
hold off



function dydt = rhs_ode(t, y,tau,k,mu,v,beta,alpha,sigma)
    S = y(1);
    I = y(2);
    p = y(3);
    M1 = y(4);
    M2 = y(5);

    dSdt = mu * (1 - p) - beta * S * I - mu * S;
    dIdt = beta * S * I - (mu + v) * I;
    dpdt = k * p * (1 - p) * (I - alpha * M2);
    dM1dt = sigma * p - sigma * M1;
    dM2dt= sigma*M1-sigma*M2;
    
    dydt = [dSdt; dIdt; dpdt; dM1dt;dM2dt];
end

