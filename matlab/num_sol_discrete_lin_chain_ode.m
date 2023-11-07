clear all
close all
clc
    k = 400;
    mu = 3.9 * 10^(-5);
    v = 1/7;
    beta = 10 * (mu + v);
    alpha = 0.002;

tau_vector=[20,50,100,150]; %constant delay 

for i = 1:length(tau_vector)
    tau= tau_vector(i);
    disp(tau);
    t_hist_s=tau; 
    t_hist_e=t_0;

    t_history=linspace(t_hist_s,t_hist_e,100);
    n_history=size(t_history,1); % dimension of the history
     disp(n_history);
    for j=1:n_history
        t_j = t_history(j);    
        
        x_history(:, j) = history_fct(t_j);
       
    end

    solution = ode45(@(t,y,Z) dynamics_discrete(t,y,Z,tau), tau, @history_fct, [t_0,4000]);
    
    subplot(4, 1, i); 
    plot(solution.x, solution.y(2,:),'b-');
    title(['\tau = ', num2str(tau)], 'FontName', 'helvetica', 'FontSize', 16);
    xlabel('time','FontName','helvetica','FontSize',16)
    ylabel('I(t)','FontName','helvetica','FontSize',16)
    xlim([0 1000])
    ylim([1.5e-4 3e-4])
end
hold off

function x_h = history_fct(t)

   x_h = [0.8; 0.1; 0];
end

function dydt = dynamics_discrete(t, Y, Z, tau) 
    k = 400;
    mu = 3.9 * 10^(-5);
    v = 1/7;
    beta = 10 * (mu + v);
    alpha = 0.002;
    

    S = Y(1);
    I = Y(2);
    p = Y(3);

    % Assuming Z corresponds to delayed values of Y at time t-tau
    p_delayed = Z(3); % delayed p for the third differential equation

    dydt = [
        mu * (1 - p) - beta * S * I - mu * S; 
        beta * S * I - (mu + v) * I;
        k * p * (1 - p) * (I - alpha * p_delayed) % Here, delayed p directly
    ];
    %disp(dydt);
end

