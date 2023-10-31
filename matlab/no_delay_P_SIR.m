clear all;
close all;
% Define the parameters
k = 40;
R_0 = 15; % Basic Reproduction Number
mu = 1/78; % per year
v = 1/7; % per day
beta = R_0*(mu+v); % Transmission rate
S0 = 1.04/R_0;
I0 = 0.82e-5;
p0 = 0.95; % starting value of p
p2008 = 0.90;
gamma = 1.253*10^(-4);        
alpha_A = gamma/0.766; 
% Define the system of differential equations
f = @(t,y) [
    mu*(1-p2008) - mu*y(1) - beta*y(1)*y(2);                                   % dS/dt
    beta*y(1)*y(2) - (mu+v)*y(2);                                             % dI/dt
     k*(1-p2008) * ((k*y(2) - alpha_A*p2008)*p2008 + gamma);            % dp/dt
];
 
% Set the initial conditions
y0 = [S0; I0; p0];

% Set the time span
tspan = [0 80*365];

% Solve the system
[t, Y] = ode45(f, tspan, y0);

% Plot the results
figure;
subplot(1,3,1);
plot(t, R_0*Y(:,1));
title('R_E(t)');
xlabel('time (years)');
ylabel('R_E');

subplot(1,3,2);
plot(t, Y(:,2));
title('I(t)');
xlabel('time (years)');
ylabel('I');

subplot(1,3,3);
plot(t, ones(size(t))*p2008);
title('p(t)');
xlabel('time (years)');
ylabel('p');
