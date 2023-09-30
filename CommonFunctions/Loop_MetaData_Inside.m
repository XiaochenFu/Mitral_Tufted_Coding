mydir  = pwd;
idcs   = strfind(mydir,'\');
newdir = mydir(1:idcs(end)-1);
DataPath = fullfile(newdir,'Behaviour_Preprocess');
ResultPath = fullfile(newdir,'Behaviour_Accuracy');
Result_Title = erase(FileName,'_BehaviourInfo.mat');
n = split(Result_Title,"_");
training_day = n{1};
pair = n{2};
DataFile = [Result_Title,'_BehaviourInfo.mat'];
load(fullfile(DataPath,DataFile))
%% plot the lick for S+ and S-, then get a threshold and a window for seperating S+ and S-
% psthAndBA(spikeTimes, eventTimes, window, psthBinSize)
eventTimes = extractfield(Behaviour_Info , 'TrialOnset');
TrialType = extractfield(Behaviour_Info , 'TrialType');
splus_index = strcmp('SPlus',TrialType);
% I want to fet the mean of the cumulative sum. Since it equals to the
% cumulative sum of the mean, I first calculate the mean using PSTH, then
% cumsum
calcWindow = [0, 9];
binSize = 0.1;
[psth_SP, bins_SP, rasterX_SP, rasterY_SP, spikeCounts_SP, ba_SP] = psthAndBA(locs_Lick, eventTimes(splus_index), calcWindow, binSize);
psth_SP_norm = psth_SP*binSize;
[psth_SM, bins_SM, rasterX_SM, rasterY_SM, spikeCounts_SM, ba_SM] = psthAndBA(locs_Lick, eventTimes(~splus_index), calcWindow, binSize);
psth_SM_norm = psth_SM*binSize;
t_FV_Water_s_all = [];
for FF = 1:length(Behaviour_Info)
    t_FV_Water_s = Behaviour_Info(FF).t_FV_Water_s;
    if isempty(t_FV_Water_s)
    else
        t_FV_Water_s_all = [t_FV_Water_s_all min(t_FV_Water_s)];
    end
end
t_FV_Water_s_avg = mean(t_FV_Water_s_all);

%% USe water time as the responding window
% %% find the time when the differences betweeen licks for S+ and S- get most

lick_difference = cumsum(psth_SP_norm)-cumsum(psth_SM_norm);
[maxidiff,maxdiffind] = max(lick_difference(bins_SM<t_FV_Water_s_avg));
if maxidiff<=1
    response_window = [0 t_FV_Water_s_avg];
    [~,maxdiffind] = max(bins_SM(bins_SM<t_FV_Water_s_avg)); % if mouse always licks more for s-, use time before FV opening as threshold
    cum_lick_SP = cumsum(psth_SP_norm);
    cum_lick_SM = cumsum(psth_SM_norm);
    lick_threshold = (cum_lick_SP(maxdiffind)+cum_lick_SM(maxdiffind))/2;
elseif max(cumsum(psth_SM_norm))<1
    response_window = [0 bins_SM(maxdiffind)];
    lick_threshold = 1;

else
    response_window = [0 bins_SM(maxdiffind)];
    cum_lick_SM = cumsum(psth_SM_norm);
    lick_threshold = cum_lick_SM(maxdiffind)+maxidiff/2;
end

%% calculate the accuracy according to a response window

for FF = 1:length(Behaviour_Info)
    lick_time = Behaviour_Info(FF).t_FV_Lick_s;
    lick_in_response_window = lick_time(lick_time>response_window(1)&lick_time<response_window(2));
    if length(lick_in_response_window)>lick_threshold
        if strcmp('SPlus',Behaviour_Info(FF).TrialType)
            Behaviour_Info(FF).Response = 'HIT';
        else
            Behaviour_Info(FF).Response = 'FA';
        end
    else
        if strcmp('SPlus',Behaviour_Info(FF).TrialType)
            Behaviour_Info(FF).Response = 'MISS';
        else
            Behaviour_Info(FF).Response = 'CR';
        end
    end
end
%% calculate the Dprime with a sliding window
Response = extractfield(Behaviour_Info , 'Response');
signal_label = 'HIT';
noise_label = 'FA';

% Extract the Behaviour_Info_cat labels from the Behaviour_Info structure
labels = {Behaviour_Info.Response};

% Calculate the hit and false alarm counts for each window
hit_idx = strcmp(labels, signal_label);
miss_idx = strcmp(labels, 'MISS');
fa_idx = strcmp(labels, noise_label);
cr_idx = strcmp(labels, 'CR');

total_length = length(labels);
window_steps = 5;
window_width = 20;

window_starts = 1:window_steps:(total_length - window_width + 1);
window_indices = bsxfun(@plus, window_starts', (0:window_width-1));

n_hit = sum(hit_idx(window_indices), 2);
n_miss = sum(miss_idx(window_indices), 2);
n_fa = sum(fa_idx(window_indices), 2);
n_cr = sum(cr_idx(window_indices), 2);

n_splus = n_hit + n_miss;
n_sminus = n_fa + n_cr;
%%
dprimes_Loglinear50 = arrayfun(@(hit, fa, splus, sminus) Dprime_Loglinear_norm50(hit, fa, splus, sminus), ...
                      n_hit, n_fa, n_splus, n_sminus);

dprimes_Loglinear = arrayfun(@(hit, fa, splus, sminus) Dprime_Loglinear(hit, fa, splus, sminus), ...
                      n_hit, n_fa, n_splus, n_sminus);

dprimes_2N = arrayfun(@(hit, fa, splus, sminus) Dprime_2N(hit, fa, splus, sminus), ...
                      n_hit, n_fa, n_splus, n_sminus);
%%
dprimes_Loglinear50_max = max(dprimes_Loglinear50);
dprimes_Loglinear_max = max(dprimes_Loglinear);
dprimes_2N_max = max(dprimes_2N);

%%
criteria.test_method = "dprime";
criteria.test_threshold = 2.5;
criteria.num_trials_for_test = 25;
criteria.num_trials_steps = 5;

criteria.drpime_function = @Dprime_Loglinear_norm50;
num_trial_Loglinear50 = Trials_to_Criteria(Behaviour_Info, criteria);

criteria.drpime_function = @Dprime_Loglinear;
num_trial_Loglinear = Trials_to_Criteria(Behaviour_Info, criteria);

criteria.drpime_function = @Dprime_2N;
num_trial_2N = Trials_to_Criteria(Behaviour_Info, criteria);


%%
HIT = strcmp('HIT',Response);
FA = strcmp('FA',Response);
MISS = strcmp('MISS',Response);
CR = strcmp('CR',Response);
% SP_accuracy = sum(HIT)/(sum(HIT)+ sum(MISS))
% SM_accuracy = sum(CR)/(sum(CR)+ sum(FA))
All_Accuracy = (sum(HIT)+ sum(CR))/(sum(HIT)+ sum(MISS)+sum(CR)+ sum(FA));
All_dprime =  Dprime_2N(sum(HIT), sum(FA), (sum(HIT)+ sum(MISS)), (sum(CR)+ sum(FA)));
fprintf('')