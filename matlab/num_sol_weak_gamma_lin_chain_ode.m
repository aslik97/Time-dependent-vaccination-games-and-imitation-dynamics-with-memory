clear all;
close all;

% Define model parameters
mu = 3.9 * 10^(-5);
v = 1/7;
beta = 10 / (mu + v); 
k = 1; 
sigma = 0.1; 
beta = 10 * (mu + v);
alpha = 0.002;
R0=10;

%  weak gamma distribution kernel function
g = @(s) sigma * exp(-sigma * s);

%  initial conditions
I0 = mu * (1 - 1/R0) / (mu + v);
S0 = 1/R0;
p0 = 0.05; 

% Define the time span for integration
tspan = [0, 4000];

% Define the system of ODEs
odefun = @(t, Y) [
    mu * (1 - Y(3)) - beta * Y(1) * Y(2) - mu * Y(1);
    beta * Y(1) * Y(2) - (mu + v) * Y(2);
    k * Y(3) * (1 - Y(3)) * (Y(2) - alpha * Y(4));
    sigma * Y(3) - sigma * Y(4);
];

% Solve the system of ODEs using ode45
[t, Y] = ode45(odefun, tspan, [S0, I0, p0, 0]);

% Extract the results
S = Y(:, 1);
I = Y(:, 2);
p = Y(:, 3);
M_values = Y(:, 4);

% Plot the results
figure;
hold on;
plot(t, I, 'r-', 'LineWidth', 2);
title('I');
xlabel('Time');
ylabel('I');
legend('I');
ylim([0 3e-4]);


hold off;
