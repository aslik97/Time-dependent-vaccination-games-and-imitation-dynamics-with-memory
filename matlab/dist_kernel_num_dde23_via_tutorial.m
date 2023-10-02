% values for the delays


% define vector with all delays
%lags=[20,50,100,150];
lags=[150];

for j = 1:length(lags)
        tau = lags(j);
  fprintf('lags(delay parameter): %d\n', tau);

solution = dde23(@rhs_dde23,tau,@dde23history, [0, 2000]);
disp(size(solution.x));

% plot the results
figure
plot(solution.x,solution.y(2,:),'b-','Linewidth',2)
%hold on
%plot(solution.x,solution.y(2,:),'r-','Linewidth',2)
%hold on
%plot(solution.x,solution.y(3,:),'g-','Linewidth',2)
hold on
title('Solution','FontSize',16)
xlabel('time','FontSize',14)
ylim([0 3e-6]);
ylabel('solution','FontSize',14)
legend('x','y','z')
hold off
end
function x_h = dde23history(t)
   % S+I=1
   S=0.1;
     x_h = [S 1-S 0.1];
end

function xdot = rhs_dde23(t,x,Z)

    k = 400;
    mu = 3.9 * 10^(-5);
    v = 1/7;
    beta = 10 * (mu + v);
    alpha = 0.002;
    t_0 = 0.01;
    dirac_function = dirac(t- t_0);
    % Assuming Z corresponds to delayed values of dynamicasl system at time t-tau
    p_delayed = Z(:,1); % delayed p for the third differential equation
    
    xdot = [
        mu * (1 - x(3)) - beta * x(1) * x(2) - mu * x(1) ;
        beta * x(1) * x(2) - (mu + v) * x(2);
        k * x(3) * (1 - x(3)) * (1 - alpha *trapz(dirac_function*p_delayed(3)));% Here, delayed p directly
    ];

end

%function g_value = discrete_kernel( t,tau_0, sigma)
    % t: The time variable
    % tau_0: The center of the Dirac delta function
    % sigma: The width of the Gaussian approximation
    
 %   g_value = exp(-(t - tau_0).^2 / (2 * sigma^2)) / (sigma * sqrt(2 * pi));
%end