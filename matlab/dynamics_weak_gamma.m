function dydt = dynamics_weak_gamma(t, Y, Z, tau,k)
    S = Y(1);
    I = Y(2);
    p = Y(3);
    % Load the parameters and function from dynamics_function.m
    %mu = 3.9*exp(-5);
    %v = 1/7;
    %beta = 10 * (mu + v);
    %alpha = 0.002;

    n_weak = 1;
    sigma = 2;
    g = @(s) (s.^(n_weak - 1) .* sigma^n_weak .* exp(-sigma * s)) / factorial(n_weak - 1);
    % Evaluate the gamma distribution kernel function
    g_value = g(t - tau - Z);
    
    dydt = [
        mu * (1 - p) - beta * S * I - mu * S;
        beta * S * I - (mu + v) * I;
        p * k * p * (1 - p) * (1 - alpha * trapz(Z, p .* g_value)) % trapz to integrate second exp. wrt. Z
    ];
end
