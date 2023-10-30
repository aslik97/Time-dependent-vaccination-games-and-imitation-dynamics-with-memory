% Parameter definitions
mu = 1/78;          % Birth and death rate (1/L, where L is life expectancy)
v = 1/7;           % Rate of recovery from infection
R_0 = 15;          % Basic Reproduction Number
beta = R_0*(mu+v); % Transmission rate
k = 40;            % Combined constant ( given k*theta = 40)
p_c= 0.56;
%p_c = 1 - 1/R_0;   % Critical elimination coverage
%p_c =0.933; 
disp(p_c);
gamma = 0;         % For the I-model
alpha_A = 1.638* 10^(-4);     %defined from me
theta_A = k;       % From the condition k*theta = 40 and taking theta(I) = theta_A*I

% Differential equations
f = @(t,y) [
    mu*(1-y(3)) - mu*y(1) - beta*y(1)*y(2);                                   % dS/dt
    beta*y(1)*y(2) - (mu+v)*y(2);                                             % dI/dt
    k*(1-y(3)) * ((theta_A*y(2) - alpha_A*y(3))*y(3) - gamma);                % dp/dt
];

% Initial conditions
S0 = 1.04/R_0;
I0 = 0.82e-5;
p0 = 0.95;

% Time span
tspan = [0, 80*365];

% using ode45
[t, Y] = ode45(f, tspan, [S0 I0 p0]);

% Plots
figure;
subplot(1,3,1);
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
plot(t, Y(:,3));
title('p(t) vs t');
xlabel('time');
ylabel('p(t)');
hold on;
line([0 80*365],[p_c p_c],'Color','red','LineStyle','--');
legend('p(t)', 'p_c');
