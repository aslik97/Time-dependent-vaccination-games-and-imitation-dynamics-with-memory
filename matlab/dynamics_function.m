function dydt = dynamics_function(t, Y, Z, tau)
    S = Y(1);
    I = Y(2);
    p = Y(3);
    
    % Evaluate the gamma distribution kernel function
    g_value = g(t - tau - Z);
    
    dydt = [
        mu * (1 - p) - beta * S * I - mu * S;
        beta * S * I - (mu + v) * I;
        p * k * p * (1 - p) * (1 - alpha * trapz(Z, p .* g_value))
    ];
end
