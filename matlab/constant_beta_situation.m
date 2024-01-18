% Define common parameters
gamma = 1.41e-4;
sigma = 0.9;
mu = 1/(11*365);
v = 1/7;
R_0 = 15;
beta = R_0*(mu+v);
p_3 = 0.90;
p_c = 0.933;
gamma_hat = 0;
alpha_hat = 0.091e-4;
k_theta = 40;

% Initial conditions
S0 = 1.04/R_0;
I0 = 0.82e-5;
p0 = 0.95;

% Time span
tspan = [0, 365*11];

% Amplitude values to test
epsilon_values = [0, 0.01, 0.1, 0.5, 1.2]; % Include 0 for constant beta

% Prepare figure for subplotting
figure;

% Loop through each amplitude value
for i = 1:length(epsilon_values)
    epsilon = epsilon_values(i);
    
    % Solve the ODEs with the current amplitude value
    [T, Z] = ode45(@(t, y) diseaseDynamics(t, y, mu, alpha_hat, k_theta, gamma, sigma, beta, v, epsilon, R_0, p_3), tspan, [S0 I0 p0]);
    
    % Subplot for R_E(t)
    subplot(length(epsilon_values), 2, i * 2 - 1);
    plot(T/365, R_0*Z(:,1));
    title(['R_E(t) vs t for \epsilon = ', num2str(epsilon)]);
    xlabel('Time (years)');
    ylabel('R_E(t)');
    
    % Subplot for I(t)
    subplot(length(epsilon_values), 2, i * 2);
    plot(T/365, Z(:,2));
    title(['I(t) vs t for \epsilon = ', num2str(epsilon)]);
    xlabel('Time (years)');
    ylabel('I(t)');
end

% Differential equations function definition
function dydt = diseaseDynamics(t, y, mu, alpha_hat, k_theta, gamma_hat, sigma, beta, v, epsilon, R_0, p_3)
    S = y(1);
    I = y(2);
    p = y(3);
    
    % Define c(t) for seasonality
    T = 365;        % The period of the function, one year
    C = 0;          % Phase shift (in days)
    c_t = epsilon * sin(2*pi*(t - C)/T) + 1;

    % Adjust beta based on epsilon
    if epsilon == 0
        beta_t = beta; % Constant beta
    else
        beta_t = c_t * beta; % Seasonal beta
    end

    % System of equations
    dSdt = mu * (1 - p_3) - beta_t * S * I - mu * S;
    dIdt = beta_t * S * I - (mu + v) * I;
    dpdt = k_theta*(1-p_3) * ((I - alpha_hat*p_3)*p_3 + gamma_hat);
    
    dydt = [dSdt; dIdt; dpdt];
end
