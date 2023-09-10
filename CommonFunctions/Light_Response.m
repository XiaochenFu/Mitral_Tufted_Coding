function [value_decreases, time_decreases, pulse_widths] = Light_Response(light_pulse, t)
    % LightResponse: Analyze dramatic increases and decreases in a signal.
    %
    % INPUTS:
    % light_pulse: Signal to analyze.
    % t: Time vector corresponding to the light_pulse signal.
    %
    % OUTPUTS:
    % value_decreases: Values of light_pulse where dramatic decreases occur.
    % time_decreases (optional): Times corresponding to the dramatic decreases.
    % pulse_widths (optional): Time differences between the dramatic increases and decreases.
    %
    % EXAMPLES:
    % [vals] = Light_Response(light_pulse, t);
    % [vals, times] = Light_Response(light_pulse, t);
    % [vals, times, widths] = Light_Response(light_pulse, t);

    % Finding dramatic decreases in the signal
    [~, locs_decrease] = findpeaks(-1*diff(light_pulse), 'MinPeakHeight', 2);

    % Extracting the times and corresponding values of the decreases
    time_decreases = t(locs_decrease); % the sample before onset
    value_decreases = light_pulse(locs_decrease-1); % the sample before going down

    % If only one output is requested, return value_decreases
    if nargout == 1
        value_decreases = value_decreases;
        return;
    end

    % Find the dramatic increase
    [~, locs_increase] = findpeaks(diff(light_pulse), 'MinPeakHeight', 1);

    % Extracting the times and corresponding values of the increases
    time_increases = t(locs_increase);

    % The pulse width would be defined as the time differences between the
    % increases and the decreases
    pulse_widths = time_decreases - time_increases;

    % If two outputs are requested, return value_decreases and time_decreases
    if nargout == 2
        value_decreases = value_decreases;
        time_decreases = time_decreases;
        return;
    end
end

%% Scripts before turning into a function 
% % % addpath('C:\Users\yycxx\Dropbox (OIST)\Fukunaga_Lab_Joined\Code\Useful_Functions')
% % light_pulse = x3_all(:,1);
% % 
% % 
% % % Finding dramatic decreases in the signal
% % [peakVals, locs_decrease] = findpeaks(-1*diff(light_pulse), 'MinPeakHeight', 1);
% % 
% % % Extracting the times and corresponding values of the decreases
% % time_decreases = t(locs_decrease);
% % value_decreases = light_pulse(locs_decrease);
% % 
% % % Plotting the original signal with marked decrease points
% % figure;
% % plot(t, light_pulse, 'b'); % Plotting the light_pulse signal
% % hold on;
% % plot(time_decreases, value_decreases, 'ro'); % Marking the decrease points with red circles
% % xlabel('Time');
% % ylabel('light_pulse Value');
% % title('Dramatic Decreases in light_pulse Signal');
% % legend('light_pulse', 'Dramatic Decrease Points');
% % grid on;
% % 
% % % Find the dramatic increase 
% % [~, locs_increase] = findpeaks(diff(light_pulse), 'MinPeakHeight', 1);
% % 
% % % Extracting the times and corresponding values of the increases
% % time_increases = t(locs_increase);
% % 
% % % The pulse width would be defined as the time differences between the
% % % increases and the decreases
% % time_decreases - time_increases
