%% Project 2 (Influenza)
clear; clc; close all;

% Parameters for influenza
R0 = 1.5;           % Basic reproduction number
gamma = 0.33;       % Recovery rate (days^-1) -> 10 days infectious window
beta = 0.50;        % Infection rate (days^-1)

% Simulation parameters
t_max = 150;        % Extended to 150 days to see the full curve with gamma = 0.10
num_stoch_runs = 5; % Number of stochastic realizations
N_sizes = [100, 10000, 1000000];         

figure('Position', [50, 100, 1500, 450]);

for p = 1:length(N_sizes)
    N = N_sizes(p);
    % Susceptible individuals is 5% of the population
    I0 = N * 0.05;
    S0 = N - I0;
    
    % Select subplot (Grid with 1 row and 3 columns)
    subplot(1, 3, p); 
    hold on;
    
    % STOCHASTIC SIMULATION (Gillespie's Direct Method)
    for run = 1:num_stoch_runs
        max_events = min(2 * N + 100, 50000); 
        
        t_stoch = zeros(max_events, 1);
        S_stoch = zeros(max_events, 1);
        I_stoch = zeros(max_events, 1);
        
        t_stoch(1) = 0; 
        S_stoch(1) = S0; 
        I_stoch(1) = I0;
        
        count = 1;
        while t_stoch(count) < t_max && I_stoch(count) > 0
            current_S = S_stoch(count);
            current_I = I_stoch(count);
            
            % Calculate Propensities
            a1 = (beta * current_S * current_I) / N;    % Infection propensity
            a2 = gamma * current_I;                     % Recovery propensity
            a0 = a1 + a2;                               % Total propensity
            
            if a0 == 0; break; end
            
            tau = -log(rand()) / a0;
            r = rand() * a0;
            count = count + 1;
            
            if count > length(t_stoch)
                t_stoch = [t_stoch; zeros(20000, 1)]; 
                S_stoch = [S_stoch; zeros(20000, 1)]; 
                I_stoch = [I_stoch; zeros(20000, 1)]; 
            end
            
            t_stoch(count) = t_stoch(count-1) + tau;
            if r < a1
                S_stoch(count) = current_S - 1;
                I_stoch(count) = current_I + 1;
            else
                S_stoch(count) = current_S;
                I_stoch(count) = current_I - 1;
            end
        end
        
        t_plot = t_stoch(1:count);
        S_plot = S_stoch(1:count);
        I_plot = I_stoch(1:count);
        
        h_stoch = plot(t_plot, I_plot / N, 'Color', [1.0, 0.4, 0.4, 0.4], 'LineWidth', 1.2);
        plot(t_plot, S_plot / N, 'Color', [0.4, 0.6, 1.0, 0.2], 'LineWidth', 1.0);
    end
    
    % DETERMINISTIC SIMULATION
    dt = 0.05; 
    time_steps = 0:dt:t_max;
    
    S_det = zeros(size(time_steps));
    I_det = zeros(size(time_steps));
    
    S_det(1) = S0 / N;
    I_det(1) = I0 / N;
    
    for step = 1:(length(time_steps) - 1)
        ds = -beta * S_det(step) * I_det(step);
        di = beta * S_det(step) * I_det(step) - gamma * I_det(step);
        
        S_det(step+1) = S_det(step) + ds * dt;
        I_det(step+1) = I_det(step) + di * dt;
    end
    
    h_ode_I = plot(time_steps, I_det, 'r-', 'LineWidth', 2.5);
    h_ode_S = plot(time_steps, S_det, 'b--', 'LineWidth', 2.0);
    
    % GRAPH FORMATTING
    title(['Influenza Simulation: N = ', num2str(N)]);
    xlabel('Time (days)');
    ylabel('Proportion of Population');
    ylim([0 1]);
    grid on;
    
    if p == 1
        legend([h_stoch, h_ode_I, h_ode_S], ...
            {'Stochastic (Infected)', 'ODE (Infected)', 'ODE (Susceptible)'}, ...
            'Location', 'best');
    end
end