% Define parameters
n_weak = 1;
n_strong = 2;
sigma_weak = 2;
sigma_strong = 2;
tau = 2;
x = 0:0.01:5; % Range of x values

% Calculate gamma distribution kernels
weak_kernel = (x.^(n_weak - 1) .* sigma_weak^n_weak .* exp(-sigma_weak * x)) / factorial(n_weak - 1);
strong_kernel = (x.^(n_strong - 1) .* sigma_strong^n_strong .* exp(-sigma_strong * x)) / factorial(n_strong - 1);

% Define parameters for AF distribution kernel
t = tau; % The parameter t
t1_values = 0.1; % Different t1 values

% Calculate AF distribution kernels
AF_kernels = zeros(length(t1_values), length(x));
for i = 1:length(t1_values)
    t1 = t1_values(i);
    t2 = t-t1;
    AF_kernels(i, :) = (exp(-x/t1) - exp(-x/t2)) / (t1 - t2);
end

% Create the plot
figure;
hold on;
plot(x, weak_kernel, 'b', 'LineWidth', 2, 'DisplayName', 'Weak Gamma (n=1)');
plot(x, strong_kernel, 'r', 'LineWidth', 2, 'DisplayName', 'Strong Gamma (n=2)');
for i = 1:length(t1_values)
    plot(x, AF_kernels(i, :), 'Color', [0.2, 0.7, 0.2], 'DisplayName', ...
         ['AF (t1 = ' num2str(t1_values(i)) ')']);
end
hold off;

% Set labels and title
xlabel('x');
ylabel('Kernel Value');
title('Weak, Strong, and AF Distribution Kernels');
legend;

% Set axis limits
xlim([0, 5]);
ylim([0, 2]);

% Show grid
grid on;
