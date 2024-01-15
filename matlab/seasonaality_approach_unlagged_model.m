clear all
close all

% Parameter definitions
mu = 1/(75*365);  % Birth and death rate
v = 1/(7);        % Rate of recovery from infection
R_0 = 10;         % Basic Reproduction Number
beta = R_0 * (mu + v); % Transmission rate
%k_theta = 40;     % Combined constant (k*theta)
p_c = 0.9;      % Critical elimination coverage
theta= 15000;
k=0.002;

% Initial conditions
S0 = 1/R_0;
I0 = (mu*(1-(1/R_0)))/(mu+v);
p0 = 0.95;        % Starting value for perception

% Time span
tspan = [0, 11 *365];

% using ode45
[t, Y] = ode45(@(t, y) dynamics(t, y, mu, v, beta,theta,k), tspan, [S0 I0 p0]);

% Plots
figure;
subplot(1,3,1);
plot(t/365, R_0*Y(:,1)); % R_E(t) = R_0*S(t)
title('R_E(t) vs t');
xlabel('time');
ylabel('R_E(t)');

subplot(1,3,2);
plot(t/365, Y(:,2));
title('I(t) vs t');
xlabel('time');
ylabel('I(t)');
ylim([0 8e-3]);

subplot(1,3,3);
plot(t/365, Y(:,3));
title('p(t) vs t');
xlabel('time');
ylabel('p(t)');
hold on;
line([0 11],[p_c p_c],'Color','red','LineStyle','--');
legend('p(t)', 'p_c');

function dydt = dynamics(t, y, mu, v, beta,theta,k)
    % Extract the current states for readability
    S = y(1); % Susceptible population
    I = y(2); % Infected population
    p = y(3); % Perception, which is equivalent to p_3 but dynamic
    sigma=0.1;
      c_t = (1 - sigma) * sin(2*pi*t/365) + 1 + sigma;
    % Calculate the derivatives
    dSdt = mu * (1 - p) - mu * S - beta * S * I *c_t ;
    dIdt = beta * S * I*c_t - (mu + v) * I;
    dpdt = k * (theta*I - p) *p*(1-p);

    % Combine the derivatives into a column vector
    dydt = [dSdt; dIdt; dpdt];
end