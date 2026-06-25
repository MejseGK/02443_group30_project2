clear; 
clc; 
close all;

N = 1000000;        % Fixed population size
t_max = 150;        % Simulation time window (days)
dt = 0.05;          % Time step for Euler integration
time_steps = 0:dt:t_max;
num_steps = length(time_steps);

% Initial conditions
i0 = 0.05;          % 5% initially infected
s0 = 1.0 - i0;      % 95% initially susceptible
r0 = 0.0;           % 0% initially recovered

% Exact pathogen configurations from your referenced papers [Beta, Gamma, Name]
pathogens = {
    0.50, 0.33, 'Seasonal Influenza';
    0.25, 0.10, 'COVID-19 (Wuhan)';
    0.15, 0.10, 'Ebola Virus'
};

num_pathogens = size(pathogens, 1);

% Preallocate arrays to store results for comparison
S_results = zeros(num_pathogens, num_steps);
I_results = zeros(num_pathogens, num_steps);
R_results = zeros(num_pathogens, num_steps);

% DETERMINISTIC SIMULATION (Euler Integration)
for p = 1:num_pathogens
    beta  = pathogens{p, 1};
    gamma = pathogens{p, 2};
    
    % Set initial states
    S_results(p, 1) = s0;
    I_results(p, 1) = i0;
    R_results(p, 1) = r0;
    
    for step = 1:(num_steps - 1)
        s_curr = S_results(p, step);
        i_curr = I_results(p, step);
        
        % Evaluate classical SIR differential equations
        ds = -beta * s_curr * i_curr;
        di = beta * s_curr * i_curr - gamma * i_curr;
        dr = gamma * i_curr;
        
        % Update states using Euler step
        S_results(p, step+1) = s_curr + ds * dt;
        I_results(p, step+1) = i_curr + di * dt;
        R_results(p, step+1) = R_results(p, step) + dr * dt;
    end
end

% Plotting
figure('Position', [100, 100, 850, 500]);
hold on;

colors = {'b-', 'r-', 'g-'}; 
h_lines = zeros(num_pathogens, 1);

for p = 1:num_pathogens
    h_lines(p) = plot(time_steps, I_results(p, :), colors{p}, 'LineWidth', 2.5);
end

xlabel('Time (days)', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('Proportion of Infected ($I/N$)', 'Interpreter', 'latex', 'FontSize', 12);
ylim([0 0.35]); 
grid on;

legend(h_lines, pathogens(:, 3), 'Interpreter', 'latex', 'Location', 'best', 'FontSize', 11);
figure('Position', [50, 150, 1500, 450]); 

for p = 1:num_pathogens
    subplot(1, 3, p); % Grid with 1 row and 3 columns
    hold on;
    
    % Plot S, I, R lines for the specific pathogen
    h_s = plot(time_steps, S_results(p, :), 'b--', 'LineWidth', 2.0);
    h_i = plot(time_steps, I_results(p, :), 'r-', 'LineWidth', 2.5);
    h_r = plot(time_steps, R_results(p, :), 'g-.', 'LineWidth', 2.0);
    
    % Formatting each individual subplot
    title(pathogens{p, 3}, 'Interpreter', 'latex', 'FontSize', 13);
    xlabel('Time (days)', 'Interpreter', 'latex', 'FontSize', 11);
    ylabel('Proportion of Population', 'Interpreter', 'latex', 'FontSize', 11);
    ylim([0 1]);
    grid on;
    
    % Place legend on each subplot for standalone clarity
    legend([h_s, h_i, h_r], ...
           {'Susceptible ($S/N$)', 'Infected ($I/N$)', 'Recovered ($R/N$)'}, ...
           'Interpreter', 'latex', 'Location', 'best', 'FontSize', 10);
end