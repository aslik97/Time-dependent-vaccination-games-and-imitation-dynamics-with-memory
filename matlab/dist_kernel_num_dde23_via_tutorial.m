% values for the delays

clear all:
close all;
% define vector with all delays
lags=[20,50,100,150];


for j = 1:length(lags)
        tau = lags(j);
  fprintf('lags(delay parameter): %d\n', tau);

solution = dde23(@(t,y,Z) rhs_dde23(t,y,Z,tau),tau,@(t) dde23history(t,tau), [0,2000]);
disp(size(solution.x));

% plot the results
figure
plot(solution.x,solution.y(2,:),'b-','Linewidth',2)
hold on
title('Solution','FontSize',16)
xlabel('time','FontSize',14)
ylim([0 3e-6]);
ylabel('solution','FontSize',14)
legend('x','y','z')
hold off
end
function x_h = dde23history(t,tau)
   % S+I=1
   S=0.1;
   I=1-S;
   x_h = [S I 0]; %TODO THE PARAMETERS SHOULD PROBABLY BE DIFFERENT!!!!!!
end

function xdot = rhs_dde23(t,x,Z,tau)

    k = 400;
    mu = 3.9 * 10^(-5);
    v = 1/7;
    beta = 10 * (mu + v);
    alpha = 0.002;
    t_0 = tau;
    % Assuming Z corresponds to delayed values of dynamical system at time t-tau
    p_delayed = Z(3); % delayed p for the third differential equation
    
    dirac_function = dirac(t- t_0);
    
    % large upper limit for integration 
    upper_limit = 1000;
    
    % the function to be integrated 
    f = @(tau) f_function(p_delayed,tau,t);

    % quad function for numerical integration
    result = quad(f, 0, upper_limit);

    
    xdot = [
        mu * (1 - x(3)) - beta * x(1) * x(2) - mu * x(1) ;
        beta * x(1) * x(2) - (mu + v) * x(2);
        k * x(3) * (1 - x(3)) * (1 - alpha *result);% Here, delayed p directly
    ];

end
function g_value = discrete_kernel( t,tau_0)
    % t: The time variable
    % tau_0: The center of the Dirac delta function
    % sigma: The width of the Gaussian approximation
    sigma = 0.1;
    g_value = exp(-(t - tau_0).^2 / (2 * sigma^2)) / (sigma * sqrt(2 * pi));
end

function f = f_function(p_delayed,tau,t)
 f= p_delayed * dirac(t- tau);
end
