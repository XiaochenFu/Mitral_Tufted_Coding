%%
% Ch1: Lick, high for off, low for on
% Ch2: FV, high for on, low for off
% Ch3: trial type, high for S-, low for S+
% Ch4: Water, high for on, low for off
% Ch5: Sniff, downside for inhalation
% Ch6: Light command
% Ch7: Recorded Sniff Phase, high for inhalation, low for exhalation

%%
% Necessary output: Trial number, Trial Type, Lick time in a time window,
% time of the first water reward
mydir  = pwd;
idcs   = strfind(mydir,'\');
newdir = mydir(1:idcs(end)-1);
DataPath = fullfile(newdir,'Recording');
ResultPath = fullfile(newdir,'Behaviour_Preprocess');

window_lick_check = [-1 9];
window_water_check = [0 10]; % check the first water onset within the FV opens

%%
DATA = readmatrix(fullfile(DataPath,FileName));
if contains(FileName,'Sampling5KHz')
    fs_all = 5000;
elseif contains(FileName,'Sampling2KHz')
    fs_all = 2000;
else
    fs_all = 1000;
end
t_all = (1:length(DATA))/fs_all;
interval_all = 1/fs_all;
spike2_latency = 0;
trial_type_low_samplerate_spike2 = DATA(:,3);
trial_type_spike2_fs = fs_all;
trial_type_t = t_all;
trial_type_interval = interval_all;

DO_comFV = DATA(:,2);
FV_fs = fs_all;

LED_command_spike2 = DATA(:,6);
LED_command_t = t_all;
LED_command_fs = fs_all;
LEDCom_inverval = interval_all;

sniff_spike2 = DATA(:,5);
sniff_fs = fs_all;
sniff_inverval = interval_all;


sniff_trigger = DATA(:,7);
sniff_trigger_fs = fs_all;
sniff_trigger_inverval = interval_all;


tD = interval_all;
tA = interval_all;
Water = DATA(:,4);
Water_fs = fs_all;
Lick = DATA(:,1);
Lick_fs = fs_all;

ChA = trial_type_low_samplerate_spike2;% Channel type.
% ChB= LineBCh.values;
%% Maximum Time between FV open and water delivery
% in Aliya's setup is about 3s.
% in Janine's setup, about 4.18s
t_FV_water = 4.18;

%% Use ChB or ChA as Experiment start
locs_ChA = find(ChA>2);
t_Exp_start = locs_ChA(1)*tA;
%% resample and binary the data
Water_Bi = sig_Bi(Water,3,0);
Lick_bi = sig_Bi(-Lick,0,-5);
ChA_Bi = sig_Bi(ChA,3,4);
% ChB_Bi = sig_Bi(-ChB, -1.5, -2);
DO_comFV_Bi = sig_Bi(DO_comFV,3.4,0);
%% find the onset of the odor
[~,locs_D,w,p] = findpeaks(DO_comFV_Bi,FV_fs);
oderlength = 0.29;
% ind_locs_D = locs_D(w>oderlength)*FV_fs;
odour_onsets = locs_D(w>oderlength);
%% find the onset of Water reward during the experiment
[~,locs_Water,w_Water,p_Water] = findpeaks(Water_Bi,Water_fs);
locs_Water_AfterOnset = locs_Water>t_Exp_start;
locs_Water_Odour = w_Water<oderlength;
pulse_length = 0; % there's a brief pulse at the beginning of the trial
locs_Water_Odour_not_briefpulse = w_Water>pulse_length;
% ind_locs_W = locs_Water(locs_Water_AfterOnset&locs_Water_Odour&locs_Water_Odour_not_briefpulse)*Water_fs;
water_onsets = locs_Water(locs_Water_AfterOnset&locs_Water_Odour&locs_Water_Odour_not_briefpulse);
%% find the onset of Licks
[~,locs_Lick,~,~] = findpeaks(Lick_bi,Lick_fs);
% lick_onsets = zeros(size(Lick_bi));
% lick_onsets(locs_Lick) = 1;
% t_FV_water = min(ind_locs_W(9)*tD-ind_locs_D(1)*tD)
%%
Behaviour_Info = struct('index',[],...
    'TrialOnset',[],...
    'TrialType',[],...
    't_FV_Lick_s',[],...
    't_FV_Water_s',[]...
    );

%%
Trial_Type = [];
Lick_Number_before_water = [];
Trial_Correct = [];
for FF = 1:length(odour_onsets)
    Behaviour_Info(FF).index = FF;
    FV_onset = odour_onsets(FF);
    Behaviour_Info(FF).TrialOnset = FV_onset;
    % Check water time
    Water_onset = water_onsets(water_onsets>(FV_onset+window_water_check(1))&water_onsets<(FV_onset+window_water_check(2)));
    Lick_onsets = locs_Lick(locs_Lick>(FV_onset+window_lick_check(1))&locs_Lick<(FV_onset+window_lick_check(2)));
    Behaviour_Info(FF).t_FV_Water_s = Water_onset-FV_onset;
    Behaviour_Info(FF).t_FV_Lick_s = Lick_onsets-FV_onset;
    if isempty(Water_onset)
        Behaviour_Info(FF).Lick_Number_before_water =  sum(length(Lick_onsets(Lick_onsets>FV_onset&Lick_onsets<(FV_onset+t_FV_water))));
    else
        Behaviour_Info(FF).Lick_Number_before_water =  sum(length(Lick_onsets(Lick_onsets>FV_onset&Lick_onsets<Water_onset(1))));
    end

    if ChA_Bi(int32(FV_onset*tD/tA)*FV_fs) % S+
        Behaviour_Info(FF).TrialType = 'SPlus';
    else
        Behaviour_Info(FF).TrialType = 'SMinus';
    end

end
%%
fprintf(['Water volume = ', num2str(sum(locs_Water_Odour)*10)])

% close()
Result_Title = Convert_Behaviour_FileName(FileName);
ResultName = [Result_Title,'_BehaviourInfo.mat'];
save(fullfile(ResultPath,ResultName), 'Behaviour_Info','locs_Lick')