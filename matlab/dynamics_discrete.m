function dydt = dynamics_discrete(t, Y, Z, tau)
    S = Y(1);
    I = Y(2);
    p = Y(3);
  
    g = @(s) exp(-s * tau);
    % Evaluate the discrete distribution kernel function
    g_value = g(t - tau - Z);
    
    dydt = [
        mu * (1 - p) - beta * S * I - mu * S;
        beta * S * I - (mu + v) * I;
        p * k * p * (1 - p) * (I - alpha * trapz(Z, p .* g_value)) % trapz to integrate second exp. wrt. Z
    ];
end


