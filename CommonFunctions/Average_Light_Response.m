function [avg_values, avg_pulse_widths] = Average_Light_Response(x3_all, t)
    % Average_Light_Response: Calculate the average value_decreases for each pulse 
    % and pulse_widths over multiple trials.
    %
    % INPUTS:
    % x3_all: Matrix where each column represents the LightResponse for a trial.
    % t: Time vector corresponding to the light_pulse signal.
    %
    % OUTPUTS:
    % avg_values: Average of the value_decreases for each pulse across trials.
    % avg_pulse_widths: Average of the pulse_widths across trials.

    n_trials = size(x3_all, 2); % Number of trials

    % Get the number of pulses from the first trial (assuming consistent across trials)
    [~, ~, sample_pulse_widths] = Light_Response(x3_all(:, 1), t);
    num_pulses = length(sample_pulse_widths);
    
    % Initialize accumulators
    total_values = zeros(1, num_pulses);
    total_pulse_widths = zeros(1, num_pulses);
    
    for i = 1:n_trials
        [value_decreases, ~, pulse_widths] = Light_Response(x3_all(:, i), t);
        
        % Update the accumulators for each pulse
%         for j = 1:num_pulses
            total_values = total_values + value_decreases';
            total_pulse_widths = total_pulse_widths + pulse_widths;
%         end
        
    end
    
    % Calculate the averages
    avg_values = total_values / n_trials;
    avg_pulse_widths = total_pulse_widths / n_trials;
end
