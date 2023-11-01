clear all;
close all;
% Parameter definitions
mu = 1/78;          % Birth and death rate (1/L, where L is life expectancy)
v = 1/7;           % Rate of recovery from infection
R_0 = 15;          % Basic Reproduction Number
beta = R_0*(mu+v); % Transmission rate
p_2008 = 0.90;     % Constant vaccination coverage

% Differential equations
f = @(t,y) [
    mu*(1-p_2008) - mu*y(1) - beta*y(1)*y(2);  % dS/dt
    beta*y(1)*y(2) - (mu+v)*y(2);              % dI/dt
];

% Initial conditions
S0 = 1.04/R_0;
I0 = 0.82e-5;

% Time span
tspan = [0, 80*365];

% Solve the differential equations using ode45
[t, Y] = ode45(f, tspan, [S0 I0]);

% Plotting
figure;
subplot(1,3,1); % Adjusted to 2 subplots
plot(t, R_0*Y(:,1)); % R_E(t) = R_0*S(t)
title('R_E(t) vs t');
xlabel('time');
ylabel('R_E(t)');

subplot(1,3,2);
plot(t, Y(:,2));
title('I(t) vs t');
xlabel('time');
ylabel('I(t)');

subplot(1,3,3);
plot(t, ones(size(t))*p_2008);
title('p(t)');
xlabel('time (years)');
ylabel('p');
