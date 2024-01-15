% Define the parameters
epsilon = 0.1;  % The amplitude of the seasonal variation
T = 365;        % The period of the function, one year
C = 0;          % Phase shift (in days)
sigma = 0.2;    % Define sigma (you need to specify an appropriate value)

% Time vector for one year with daily resolution
t = 0:1:T-1;  % time vector from day 0 to day 364

% Define c(t)
c_t = epsilon * sin(2*pi*(t - C)/T) + 1;
c_t_old = (1 - sigma) * sin(2*pi*t/365) + 1 + sigma;

% Plot c(t), c_t, and c_t_old in one figure
figure;
plot(t, c_t, 'b');      % Plot c_t in blue
hold on; 
plot(t, ones(size(t)), 'r--'); % Plot c(t)=1 in red dashed line

hold off; 

% Add labels and title
xlabel('Time (days)');
ylabel('c(t)');
title('Seasonal Variation c(t) over One Year');
grid on;
legend('c(t)', 'c(t)=1'); % Add a legend

