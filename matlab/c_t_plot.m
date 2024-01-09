% Define the parameters
epsilon = 0.1;  % The amplitude of the seasonal variation
T = 365;        % The period of the function, one year
C = 0;          % Phase shift (in days)

% Time vector for one year with daily resolution
t = 0:1:T-1;  % time vector from day 0 to day 364

% Define c(t)
c_t = epsilon * sin(2*pi*(t - C)/T) + 1;

figure;
plot(t, c_t);
hold on; 
plot(t, ones(size(t)), 'r--'); 
hold off; 
% Plot c(t)
figure;
plot(t, c_t);
xlabel('Time (days)');
ylabel('c(t)');
title('Seasonal Variation c(t) over One Year');
grid on;
