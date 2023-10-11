% values for the delays

clear all:
close all;
% define vector with all delays
lags=[20,50,100,150];


for j = 1:length(lags)
        tau = lags(j);
    fprintf('lags(delay parameter): %d\n', tau);
    %options = ddeset('MaxStep', 0.01); % Adjust the value as needed
    solution = dde23(@(t,y,Z) rhs_dde23(t,y,Z,tau),tau,@(t) dde23history(t,tau), [0,4000]);
    disp(solution.y(2,:));

% plot the results
figure
plot(solution.x,solution.y(2,:),'b-','Linewidth',2)
hold on
title('Solution','FontSize',16)
xlabel('time','FontSize',14)
ylim([1.5e-4 3e-4]);
ylabel('solution','FontSize',14)
legend('x','y','z')
hold off
end

function x_h = dde23history(t, tau)
   persistent saved_I; % persistent variable to store I
    S = 0.8;

    if t <= tau
        % For t <= tau, smoothly varying initial conditions
        p = 0; % initial value for p
        % smoothly varying sine function for I
        I = 0.5 * (1 - sin(pi * t / tau)); % Adjust the frequency as needed
        saved_I = I;
    else
        % For t > tau, use the saved value of I with an oscillatory component
        p = rand() * (1 - S);
        % Add an oscillatory term to I
        oscillation_frequency = 2 * pi; % Adjust the frequency 
        I = saved_I + 0.1 * sin(oscillation_frequency * (t - tau)); % Amplitude and frequency
        I = abs(I); % Ensure non-negative values
    end

    x_h = [S I p];
end



function xdot = rhs_dde23(t,y,Z,tau)

    k = 400;
    mu = 3.9 * 10^(-5);
    v = 1/7;
    beta = 10 * (mu + v);
    alpha = 0.002;
   
    % Assuming Z corresponds to delayed values of dynamical system at time t-tau
    p_delayed = Z(3); % delayed p for the third differential equation
   
    % function to be integrated 
    f = @(tau) f_function(p_delayed,tau,t);

    % integral function for numerical integration
    result = integral(f, 0, inf);
  
    xdot = [
        mu * (1 - y(3)) - beta * y(1) * y(2) - mu * y(1) ;
        beta * y(1) * y(2) - (mu + v) * y(2);
        k * y(3) * (1 - y(3)) * (1 - alpha * result);% Here, delayed p directly
    ];

end


function f = f_function(p_delayed,tau,t)% Gaussian approximation
sigma = 0.01;  % Small standard deviation
dirac_approximation = 1 / (sigma * sqrt(2 * pi)) * exp(-(t- tau).^2 / (2 * sigma^2));
 f= p_delayed * dirac_approximation;

end
