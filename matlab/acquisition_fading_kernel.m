% Define parameters
mu = 3.9 * exp(-5);
v = 1/7;
beta = 10 * (mu + v);
alpha = 0.002;
tau_0 = 2;

% Define the range of k and tau values
k_values = 0.1:0.1:10;
tau_values = 0.1:0.1:10;

% Initialize stability matrix
stability_matrix = zeros(length(k_values), length(tau_values));

for i = 1:length(k_values)
    k = k_values(i);

    for j = 1:length(tau_values)
        tau = tau_values(j);

        % Calculate characteristic equation coefficients
        T1 = (0.5 + 0.4) * tau;
        T2 = (0.5 - 0.4) * tau;

        % Define the Laplace transform of the kernel
        Lg = @(lambda) 1 / ((1 + lambda * T1) * (1 + lambda * T2));

        % Calculate endemic steady state values
        E3_S = (mu + v) / beta;
        E3_I = (mu * alpha * (beta - (mu + v))) / (beta * (mu + alpha * (mu + v)));
        E3_p = (mu * (beta - (mu + v))) / (beta * (mu + alpha * (mu + v)));

        % Calculate characteristic equation coefficients
        a = k * mu * beta * E3_p * (1 - E3_p) * E3_I;
        b = k * alpha * E3_p * (1 - E3_p) * Lg(1i) + 1i;
        c = [1, (mu + beta * E3_I), beta^2 + E3_S * E3_I];

        % Calculate eigenvalues of the characteristic equation
        eigenvalues = roots(c);

        % Check eigenvalues for stability
        if all(real(eigenvalues) < 0)
            stability_matrix(i, j) = 1;  % Stable
        else
            stability_matrix(i, j) = 0;  % Unstable
        end
    end
end

% Plot the regions of stability and instability
figure;
imagesc(tau_values, k_values, stability_matrix');
colormap([1 0 0; 0 1 0]); % Red for unstable, green for stable
colorbar('Ticks', [0, 1], 'TickLabels', {'Unstable', 'Stable'});
xlabel('\tau');
ylabel('k');
title('Stability Analysis for E3 Steady State');
