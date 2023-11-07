clear all
close all

% Parameter definitions
mu = 1/(75*365);  % Birth and death rate
v = 1/(7);        % Rate of recovery from infection
R_0 = 10;         % Basic Reproduction Number
beta = R_0 * (mu + v); % Transmission rate
p_c = 0.9;      % Critical elimination coverage
theta= 15000;
k_values=[0.0005,0.002,0.000035]; % Array of k values

% Initial conditions
S0 = 1/R_0;
I0 = (mu*(1-(1/R_0)))/(mu+v);
p0 = 0.95;        

% Time span
tspan = [0, 40 *365];
for i = 1:length(k_values)
    k = k_values(i);
% using ode45
[t, Y] = ode45(@(t, y) dynamics(t, y, mu, v, beta,theta,k), tspan, [S0 I0 p0]);

% Plots
figure;

subplot(1,2,1);
plot(t/365, Y(:,2));
title(['I(t) with k = ' num2str(k) ' days^{-1}']);
xlabel('time');
ylabel('I(t)');
ylim([0 8e-3]);

subplot(1,2,2);
plot(t/365, Y(:,3));
 title(['p(t) with k = ' num2str(k) ' days^{-1}']);
xlabel('time');
ylabel('p(t)');
hold on;
line([0 40],[p_c p_c],'Color','red','LineStyle','--');
legend('p(t)', 'p_c');
end
function dydt = dynamics(t, y, mu, v, beta,theta,k)
   
    S = y(1); % Susceptible population
    I = y(2); % Infected population
    p = y(3); % Perception, which is equivalent to p_3 but dynamic

    % derivatives
    dSdt = mu * (1 - p) - mu * S - beta * S * I;
    dIdt = beta * S * I - (mu + v) * I;
    dpdt = k * (theta*I - p) *p*(1-p);

    % Combine the derivatives into a column vector
    dydt = [dSdt; dIdt; dpdt];
end
