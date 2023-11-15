clear all
close all
% Parameter definitions
mu = 1/(78*365);  % Birth and death rate (1/L, where L is life expectancy)
v = 1/(7);% Rate of recovery from infection
R_0 = 15;% Basic Reproduction Number
beta = R_0*(mu+v); % Transmission rate
p_3 = 0.90;
k_theta = 40; % Combined constant ( given k*theta = 40)          
p_c =0.933; %critical elimination coverage
gamma_hat = 1.253e-4;%for G-model
alpha_hat= 0.091e-4; 
 
% Initial conditions
S0= 1.04/R_0;
%I0=(mu*(1-(1/R_0)))/(mu+v);
%S0 = 1.04/R_0;
I0 = 0.82e-5;
p0 = 0.95;

% Time span
tspan = [0, 80*365];

% using ode45
[t, Y] = ode45(@(t,y) dynamics(t, y, mu, v, beta, k_theta, alpha_hat, gamma_hat), tspan, [S0 I0 p0]);

% Plots
figure;
subplot(1,3,1);
plot(t/365, R_0*Y(:,1)); % R_E(t) = R_0*S(t)
title('R_E(t) vs t');
xlabel('time');
ylabel('R_E(t)');
xlim([0 80]);

subplot(1,3,2);
plot(t/365, Y(:,2));
title('I(t) vs t');
xlabel('time');
ylabel('I(t)');

subplot(1,3,3);
plot(t/365, Y(:,3));
title('p(t) vs t');
xlabel('time');
ylabel('p(t)');
hold on;
line([0 80],[p_c p_c],'Color','red','LineStyle','--');
legend('p(t)', 'p_c');
