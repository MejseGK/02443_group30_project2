% Project 2,part 2(a) - Extended SEIR Model for COVID-19
clear; clc; close all;

% Parameters for COVID-19 
gamma = 0.10;         % Recovery rate from table (days^-1) -> 10 days infectious
beta_I = 0.25;        % Symptomatic infection rate from table (days^-1)

% New parameters required for the Part II (a) extension
sigma = 1 / 4;        % Incubation rate (average of 4 days spent in E before symptoms)
beta_E = beta_I * 0.5; % Pre-symptomatic contagion rate is half of symptomatic rate

% Simulation parameters
t_max = 150;          % Extended timeline to see the full SEIR curve progression
num_stoch_runs = 5; 
N_sizes = [100, 10000, 1000000];           

figure('Position', [50, 100, 1500, 450]);

for p = 1:length(N_sizes)
    N = N_sizes(p);
    
    % Initial conditions: 5% of the population starts exposed/pre-symptomatic
    E0 = N * 0.05;   
    I0 = 0;              
    S0 = N - E0;
    
    subplot(1, 3, p); 
    hold on;
    
    % STOCHASTIC SIMULATION (Gillespie's Direct Method)
    for run = 1:num_stoch_runs
        max_events = min(3 * N + 100, 60000); 
        
        t_stoch = zeros(max_events, 1);
        S_stoch = zeros(max_events, 1);
        E_stoch = zeros(max_events, 1);
        I_stoch = zeros(max_events, 1);
        
        t_stoch(1) = 0; 
        S_stoch(1) = S0; 
        E_stoch(1) = E0;
        I_stoch(1) = I0;
        
        count = 1;
        while t_stoch(count) < t_max && (E_stoch(count) > 0 || I_stoch(count) > 0)
            current_S = S_stoch(count);
            current_E = E_stoch(count);
            current_I = I_stoch(count);
            
            % Calculate Propensities using both beta_E and beta_I
            a1 = (beta_E * current_S * current_E / N) + (beta_I * current_S * current_I / N); 
            a2 = sigma * current_E;                                                       
            a3 = gamma * current_I;                                                       
            a0 = a1 + a2 + a3;                                                             
            
            if a0 == 0; break; end
            
            tau = -log(rand()) / a0;
            r = rand() * a0;
            count = count + 1;
            
            if count > length(t_stoch)
                t_stoch = [t_stoch; zeros(20000, 1)]; 
                S_stoch = [S_stoch; zeros(20000, 1)]; 
                E_stoch = [E_stoch; zeros(20000, 1)]; 
                I_stoch = [I_stoch; zeros(20000, 1)]; 
            end
            
            t_stoch(count) = t_stoch(count-1) + tau;
            
            if r < a1
                % S -> E (Infection)
                S_stoch(count) = current_S - 1;
                E_stoch(count) = current_E + 1;
                I_stoch(count) = current_I;
            elseif r < (a1 + a2)
                % E -> I (Symptom Onset)
                S_stoch(count) = current_S;
                E_stoch(count) = current_E - 1;
                I_stoch(count) = current_I + 1;
            else
                % I -> R (Recovery)
                S_stoch(count) = current_S;
                E_stoch(count) = current_E;
                I_stoch(count) = current_I - 1;
            end
        end
        
        t_plot = t_stoch(1:count);
        I_plot = I_stoch(1:count);
        E_plot = E_stoch(1:count);
        
        h_stoch_I = plot(t_plot, I_plot / N, 'Color', [1.0, 0.4, 0.4, 0.3], 'LineWidth', 1.0);
        h_stoch_E = plot(t_plot, E_plot / N, 'Color', [1.0, 0.7, 0.2, 0.2], 'LineWidth', 1.0);
    end
    
    % DETERMINISTIC SIMULATION
    dt = 0.05; 
    time_steps = 0:dt:t_max;
    
    S_det = zeros(size(time_steps));
    E_det = zeros(size(time_steps));
    I_det = zeros(size(time_steps));
    
    S_det(1) = S0 / N;
    E_det(1) = E0 / N;
    I_det(1) = I0 / N;
    
    for step = 1:(length(time_steps) - 1)
        % Expanded differential equations with dual transmission parameters
        ds = - (beta_E * S_det(step) * E_det(step)) - (beta_I * S_det(step) * I_det(step));
        de = (beta_E * S_det(step) * E_det(step)) + (beta_I * S_det(step) * I_det(step)) - sigma * E_det(step);
        di = sigma * E_det(step) - gamma * I_det(step);
        
        S_det(step+1) = S_det(step) + ds * dt;
        E_det(step+1) = E_det(step) + de * dt;
        I_det(step+1) = I_det(step) + di * dt;
    end
    
    h_ode_I = plot(time_steps, I_det, 'r-', 'LineWidth', 2.5);
    h_ode_E = plot(time_steps, E_det, 'g-', 'LineWidth', 2.0);
    
 % GRAPH FORMATTING
    title(['COVID-19 SEIR Model: N = ', num2str(N)], 'FontSize', 16); % Extra large for title
    xlabel('Time (days)', 'FontSize', 14);
    ylabel('Proportion of Population', 'FontSize', 14);

    if p == 1
        legend([h_ode_E, h_ode_I, h_stoch_E, h_stoch_I], ...
            {'ODE (Pre-Symptomatic E)', 'ODE (Symptomatic I)', 'Stoch (E)', 'Stoch (I)'}, ...
            'Location', 'best');
    end
end