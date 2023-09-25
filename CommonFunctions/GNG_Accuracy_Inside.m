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
% plot(bins_SP,psth_SP)
% plot(rasterX, rasterY,'.')
plot(bins_SP,cumsum(psth_SP_norm),'r')
hold on
[psth_SM, bins_SM, rasterX_SM, rasterY_SM, spikeCounts_SM, ba_SM] = psthAndBA(locs_Lick, eventTimes(~splus_index), calcWindow, binSize);
psth_SM_norm = psth_SM*binSize;
plot(bins_SM,cumsum(psth_SM_norm),'b')
xlim([0 5])
legend({'S+','S-'})
xlabel('time from FV(s)')
ylabel('cumulated lick')
% xline(3.45)
% automatically get a threshold
% get the time when water comes
t_FV_Water_s_all = [];
for FF = 1:length(Behaviour_Info)
    t_FV_Water_s = Behaviour_Info(FF).t_FV_Water_s;
    if isempty(t_FV_Water_s)
    else
        t_FV_Water_s_all = [t_FV_Water_s_all min(t_FV_Water_s)];
    end
end
t_FV_Water_s_avg = mean(t_FV_Water_s_all);

title([training_day ' ' pair ' Mean Cumsum Lick'])
filename = sprintf('%s_%s_Mean_Cumsum_Lick.jpg',training_day,pair);
saveas(gcf,fullfile(ResultPath,filename))
%% USe water time as the responding window
% %% find the differences betweeen licks for S+ and S- get at the end of
% the responding window
figure
plot(bins_SM,cumsum(psth_SP_norm)-cumsum(psth_SM_norm))
xlim([0 5])
hold on 
xline(t_FV_Water_s_avg)
lick_difference = cumsum(psth_SP_norm)-cumsum(psth_SM_norm);
lick_difference_before_water = lick_difference(bins_SM<t_FV_Water_s_avg);
maxidiff = lick_difference_before_water(end);
maxdiffind = length(lick_difference_before_water);
if maxidiff<=1
    response_window = [0 t_FV_Water_s_avg];
    [~,maxdiffind] = max(bins_SM(bins_SM<t_FV_Water_s_avg)); % if mouse always licks more for s-, use time before FV opening as threshold
    cum_lick_SP = cumsum(psth_SP_norm);
    cum_lick_SM = cumsum(psth_SM_norm);
    lick_threshold = (cum_lick_SP(maxdiffind)+cum_lick_SM(maxdiffind))/2;
elseif max(cumsum(psth_SM_norm(bins_SM<t_FV_Water_s_avg)))<1
    response_window = [0 bins_SM(maxdiffind)];
    lick_threshold = 0.5;
else
    response_window = [0 bins_SM(maxdiffind)];
    cum_lick_SM = cumsum(psth_SM_norm);
    lick_threshold = cum_lick_SM(maxdiffind)+maxidiff/2;
end

% lick_threshold = 0.5;
xlabel('time from FV(s)')
ylabel('differences in cumulated lick')
% lick_in_response_window = lick_time(lick_time>response_window(1)&lick_time<response_window(2));
title([training_day ' ' pair ' Diff Mean Cumsum Lick'])
filename = sprintf('%s_%s_Diff_Mean_Cumsum_Lick.jpg',training_day,pair);
% saveas(gcf,fullfile(ResultPath,filename))
saveimg(gcf,ResultPath,Result_Title,filename,001)

%% Sanity check: show raster plot of licks, response window
figure
f = tiledlayout(1,2);
% Tile 1
nexttile
% lick raster for all trials 
[binnedArray, bins] = timestampsToBinned(locs_Lick, eventTimes, binSize, calcWindow);
% lick raster for only s plus trials
binnedArray_SP = binnedArray;
binnedArray_SP(~splus_index,:) = 0;
[tr,b] = find(binnedArray_SP);
[rasterX,yy] = rasterize(bins(b));
rasterY = yy+reshape([zeros(size(tr'));tr';zeros(size(tr'))],1,length(tr)*3); 
rasterX(rasterY==0) = [];rasterY(rasterY==0) = [];% remove zeros
plot(rasterX, rasterY,'r.')
hold on 
xline(t_FV_Water_s_avg)
ylim([0, length(Behaviour_Info)]);xlim([0 8]);
% Tile 1
nexttile
% same thing for sminus
binnedArray_SM = binnedArray;
binnedArray_SM(splus_index,:) = 0;
[tr,b] = find(binnedArray_SM);
[rasterX,yy] = rasterize(bins(b));
rasterY = yy+reshape([zeros(size(tr'));tr';zeros(size(tr'))],1,length(tr)*3); 
rasterX(rasterY==0) = [];rasterY(rasterY==0) = [];% remove zeros
plot(rasterX, rasterY,'b.')
hold on 
xline(t_FV_Water_s_avg)
ylim([0, length(Behaviour_Info)]);xlim([0 8]);
% add legend
L1 = plot(nan, nan, 'r.');
L2 = plot(nan, nan, 'b.');
legend([L1, L2], {'S+','S-'})
% title(f,[training_day ' ' pair ' Lick Rasters'],'FontSize',20)
% xlabel(f,'Time from FV opening','FontSize',18)
% ylabel(f,'Trial','FontSize',18)
title(f,[training_day ' ' pair ' Lick Rasters'],'FontSize',30,'FontName', 'Times New Roman','FontWeight','bold')
xlabel(f,'Time from FV opening','FontSize',28,'FontName', 'Times New Roman')
ylabel(f,'Trial','FontSize',28,'FontName', 'Times New Roman')
filename = sprintf('%s_%s_LIck_Raster.jpg',training_day,pair);
saveas(gcf,fullfile(ResultPath,filename))
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
%% calculate the Dprime wiht a sliding window every 5 trials
figure
Response = extractfield(Behaviour_Info , 'Response');
HIT = strcmp('HIT',Response);
FA = strcmp('FA',Response);

signal_label = 'HIT';
noise_label = 'FA';
labels = {Behaviour_Info.Response};

% Calculate the hit and false alarm counts for each window
hit_idx = strcmp(labels, signal_label);
miss_idx = strcmp(labels, 'MISS');
fa_idx = strcmp(labels, noise_label);
cr_idx = strcmp(labels, 'CR');

total_length = length(labels);
window_steps = 5;
window_width = 20;
n_windows = length(1:window_steps:(total_length - window_width + 1));
dprimes = zeros(1,n_windows);
Session = 1:window_steps:(total_length - window_width + 1);
k = 1;
for i = 1:window_steps:(total_length - window_width + 1)
    window_indices = i:(i + window_width - 1);
    n_hit = sum(hit_idx(window_indices));
    n_miss = sum(miss_idx(window_indices));
    n_fa = sum(fa_idx(window_indices));
    n_cr = sum(cr_idx(window_indices));
    n_splus = n_hit + n_miss;
    n_sminus = n_fa + n_cr;
    dprimes(k) = Dprime_2N(n_hit, n_fa, n_splus, n_sminus);
    k = k+1;
end
plot(Session,dprimes,'k')
ylim([-1.5 4])
xlabel('Trials')
ylabel('d prime')
title([training_day ' ' pair ' Learning Curve'])
filename = sprintf('%s_%s_Learning_Curve.jpg',training_day,pair);
saveas(gcf,fullfile(ResultPath,filename))
filename = sprintf('%s_%s_Learning_Curve.svg',training_day,pair);
saveas(gcf,fullfile(ResultPath,filename))

MISS = strcmp('MISS',Response);
CR = strcmp('CR',Response);
SP_accuracy = sum(HIT)/(sum(HIT)+ sum(MISS))
SM_accuracy = sum(CR)/(sum(CR)+ sum(FA))
All_Accuracy = (sum(HIT)+ sum(CR))/(sum(HIT)+ sum(MISS)+sum(CR)+ sum(FA))
%%
resultname = [Result_Title,'_BehaviourResult.mat'];
save(fullfile(ResultPath,resultname),'Behaviour_Info')