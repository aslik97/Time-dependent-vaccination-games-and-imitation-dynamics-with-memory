clear all;
close all;

% Define model parameters
mu = 0.000039; % Example value, adjust as needed
v = 1 / 7; % Example value, adjust as needed
beta = 10 * (mu + v); % Example value, adjust as needed
k = 1; % Example value, adjust as needed
alpha = 0.002; % Example value, adjust as needed
sigma = 0.1; % Example value, adjust as needed

% Define the weak gamma distribution kernel function
g = @(s) sigma * exp(-sigma * s);

% Define the initial conditions
S0 = 0.8; % Example initial value, adjust as needed
I0 = 0.1; % Example initial value, adjust as needed
p0 = 0.05; % Example initial value, adjust as needed

% Define the time span for integration
tspan = [0, 4000];

% Define the system of ODEs
odefun = @(t, Y) [
    mu * (1 - Y(3)) - beta * Y(1) * Y(2) - mu * Y(1);
    beta * Y(1) * Y(2) - (mu + v) * Y(2);
    k * Y(3) * (1 - Y(3)) * (Y(2) - alpha * (sigma * Y(3) - sigma * Y(4)));
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
subplot(2, 1, 1);
plot(t, S, 'b-', 'LineWidth', 2);
hold on;
plot(t, I, 'r-', 'LineWidth', 2);
title('S and I');
xlabel('Time');
ylabel('S and I');
legend('S', 'I');
ylim([0, 1]);

subplot(2, 1, 2);
plot(t, p, 'g-', 'LineWidth', 2);
hold on;
plot(t, M_values, 'm-', 'LineWidth', 2);
title('p and M');
xlabel('Time');
ylabel('p and M');
legend('p', 'M');

hold off;
