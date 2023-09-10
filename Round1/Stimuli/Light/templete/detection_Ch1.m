% Clear variables and close any existing figures.
clear
close all
clc

% Define core, target intensities, driving currents, and pulse duration.
Core = 2;
target_intensity1 = 25;
target_intensity2 = 25;
driving_current1 = 191;
driving_current2 = 192;
pulse_duration = 5;% Duration of pulse in milliseconds.
saveplot = 0; % Set to 1 if plots should be saved, 0 otherwise.

% Add required paths for further code execution.
addpath('C:\Users\yycxx\Dropbox (OIST)\Fukunaga_Lab_Joined\Code\Basic_Settings')
addpath('C:\Users\yycxx\Dropbox (OIST)\Fukunaga_Lab_Joined\Code\Useful_Functions')

% Run the startup script and set colour configurations.
startup
Colours

c_Sp = c_M1_blue;
c_Sm = c_M2_orange;

% Read data from CSV file and load conversion from Ch1 (100mW range) to
% mW
DATA = readmatrix(['Light_Core2_19mA_192mA_Ch1_230909131059.csv']);
load('MISC_Ch1.mat')

% Define thresholds and sample rates.
spike2_LED_threshold = 0.5;
spike2_trailtype_threshold = 0.5;
min_pulse_duration_second = 0.001;
min_trail_duration_second = 0.4;

fs_all = 5000;
t_all = (1:length(DATA))/fs_all;
interval_all = 1/fs_all;
fs_spike2 = fs_all;

% Define which core data to extract based on the value of 'Core'.
switch Core
    case 1
        LED_command_spike2 = DATA(:,3); % medal core
    case 2
        LED_command_spike2 = DATA(:,4); % latel core
end

LED_time = t_all;
LED_fs = fs_all;
LED_command_interval = interval_all;
Optometer = DATA(:,1);
Optometer_time = t_all;
Optometer_interval = interval_all;
% Convert optometer data to light power in milliwatts.
Light_Power_mW = fitobject(Optometer);
window = [-0.01 0.1];
c = 'k';
global showplot
trail_type_low_samplerate_spike2 = DATA(:,2);
trail_type_low_samplerate_spike2_time_aligned = LED_time;

% Extract onset times and trail types from voltage data.
[pulse_onset_spike2, pulse_type_spike2] = Trail_Type_From_voltage(LED_command_spike2, spike2_LED_threshold, fs_spike2, min_pulse_duration_second);
[trail_onset_spike2, trail_type_spike2] = Trail_Type_From_voltage(LED_command_spike2, spike2_trailtype_threshold, fs_spike2, min_pulse_duration_second, min_trail_duration_second);

% Align trail onset and extract stimulus information.
trail_onset_spike2_aligned = trail_onset_spike2;
Stimuli_Info = struct('TrailType', []);

for FF = 1:length(trail_onset_spike2_aligned)
    trail_onset_spike2_FF = trail_onset_spike2_aligned(FF);
    [~, sniff_onset_FF_trial_type] = read_value_from_time_point(trail_onset_spike2_FF, trail_type_low_samplerate_spike2, trail_type_low_samplerate_spike2_time_aligned);
    if sniff_onset_FF_trial_type > spike2_trailtype_threshold
        Stimuli_Info(FF).TrailType = 'SPlus';
    else
        Stimuli_Info(FF).TrailType = 'SMinus';
    end
end

% Plot the average light response for both trail types.
figure
TrailType = extractfield(Stimuli_Info, 'TrailType');
eventTimes1 = trail_onset_spike2;
subplot(2, 1, 1)
[SP_all, ~] = segment_with_onset_time(Light_Power_mW, Optometer_time, eventTimes1(strcmp(TrailType, 'SPlus')), window);

[SP_all,~] = segment_with_onset_time(Light_Power_mW, Optometer_time,  eventTimes1(strcmp(TrailType,'SPlus')), window);
x3 = mean(SP_all,2);
x3_std = std(SP_all,0,2);
itv = Optometer_interval;
t = ((1:length(x3))-1)*itv+window(1);
plot_mean_std(t, x3, x3_std, c_Sp);
yLimit2 = get(gca,'YLim'); hold on
ylabel("Output (mW)")
xlabel('Time from pulse onset (s)')
tgt = yline(target_intensity1,'r');
legend(tgt,[num2str(target_intensity1) 'mW'])
title('SPlus')
[SP_vtg, SP_widths] = Average_Light_Response(SP_all, t)

% For SMinus plot
subplot(2, 1, 2)
[SM_all,~] = segment_with_onset_time(Light_Power_mW, Optometer_time,  eventTimes1(strcmp(TrailType,'SMinus')), window);
x3 = mean(SM_all,2);
x3_std = std(SM_all,0,2);
itv = Optometer_interval;
t = ((1:length(x3))-1)*itv+window(1);
plot_mean_std(t, x3, x3_std, c_Sm);
yLimit2 = get(gca,'YLim'); hold on
ylabel("Output (mW)")
xlabel('Time from pulse onset (s)')
tgt = yline(target_intensity2,'r');
legend(tgt,[num2str(target_intensity2) 'mW'])
title('SMinus')
[SM_vtg, SM_widths] = Average_Light_Response(SM_all, t)

% Save the plot if 'saveplot' is set to 1.
output_filename = sprintf('Core%d_%dmA%dmA_Ch1_SpSm_%dms.jpg', Core, driving_current1, driving_current2, pulse_duration);
if saveplot
    saveas(gcf, output_filename)
end

% Copy and rename the current script for logging.
currentScript = strcat(mfilename('fullpath'), '.m');
[path, ~, ext] = fileparts(currentScript);
outname = sprintf('Core%d_%dmA%dmA_Ch1_%dms', Core, driving_current1, driving_current2, pulse_duration);
newScript = fullfile(path, [outname, ext]);
copyfile(currentScript, newScript);
