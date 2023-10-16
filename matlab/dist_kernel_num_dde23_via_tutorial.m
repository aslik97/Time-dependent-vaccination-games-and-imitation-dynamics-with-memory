clear all:
close all;
% define vector with all delays
lags=[20,50,100,150];
%lags=[20];
for j = 1:length(lags)
        tau = lags(j);
 %options = ddeset('RelTol',1e-6,'AbsTol',1e-8);
solution = dde23(@(t,x,Z) rhs_dde23(t,x,Z,tau),tau,@dde23history, [0,4000]);
disp(size(solution.y(2,:)));
% plot the results
figure;
   plot(solution.x,solution.y(2,:),'r-','Linewidth',2)
    title(['\tau = ', num2str(tau)], 'FontSize', 16);
    xlabel('time, days', 'FontSize', 14);
    ylabel('I(t)', 'FontSize', 14);
    ylim([0 3e-4]);
hold off
end

function x_h = dde23history(t)
    mu = 3.9 * 10^(-5);
    %mu=1/(70*365);
    v = 1/7;
    R0= 10;
    % Default case for t < 0
    p = 0;
    I = mu * (1 - 1/R0) / (mu + v);
    S = 1/R0;
    %I= 0.2;
    %S=1-I;
    %disp(I);
    if t == 0
    %   a = 0.1;      % Control the steepness. Adjust as needed.
    %t0 = 0;      % The time when p = 0.5. Adjust as needed.
    %p = 1 / (1 + exp(-a * (t - t0)));
    %disp(p);
    p=0.00001;
    %I = mu * (1 - 1/R0) / (mu + v);
    %S = 1/R0;
    end
    x_h = [S I p];
end


function xdot = rhs_dde23(t,x,Z,tau)
    k = 400;
    mu = 3.9 * 10^(-5);
    %mu=1/(70*365);
    v = 1/7;
    R0=10;
    beta = R0 * (mu + v);
    alpha = 0.002;
   
    % Assuming Z corresponds to delayed values of dynamicasl system at time t-tau
    p_delayed = Z(3); % delayed p for the third differential equation
    
    % Discrete delay using Dirac delta
    if t < tau
        integral_term = 0;
    else
        integral_term = p_delayed;
    end
    xdot = [
        mu * (1 - x(3)) - beta * x(1) * x(2) - mu * x(1) ;
        beta * x(1) * x(2) - (mu + v) * x(2);
        k * x(3) * (1 - x(3)) * (x(2) - alpha * integral_term);% Here, delayed p directly
    ];
    %disp(['Time: ', num2str(t), ' Delayed Time: ', num2str(t-tau), 'p_delayed: ', num2str(p_delayed)]);
end
