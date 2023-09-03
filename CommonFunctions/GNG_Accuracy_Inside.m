mydir  = pwd;
idcs   = strfind(mydir,'\');
newdir = mydir(1:idcs(end)-1);
DataPath = fullfile(newdir,'Behaviour_Preprocess');
ResultPath = fullfile(newdir,'Behaviour_Accuracy');
% Result_Title = strrep(FileName,'_BehaviourInfo.mat',[]);
Result_Title = erase(FileName,'_BehaviourInfo.mat');
n = split(Result_Title,"_");
training_day = n{1}; 
pair = n{2}; 
% Result_Title = Convert_Behaviour_FileName(FileName);

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
% %% find the time when the differences betweeen licks for S+ and S- get most
figure
plot(bins_SM,cumsum(psth_SP_norm)-cumsum(psth_SM_norm))
xlim([0 5])
hold on 
xline(t_FV_Water_s_avg)
% about 3.45
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

% lick_threshold = 0.5;
xlabel('time from FV(s)')
ylabel('differences in cumulated lick')
% lick_in_response_window = lick_time(lick_time>response_window(1)&lick_time<response_window(2));
title([training_day ' ' pair ' Diff Mean Cumsum Lick'])
filename = sprintf('%s_%s_Diff_Mean_Cumsum_Lick.jpg',training_day,pair);
% saveas(gcf,fullfile(ResultPath,filename))
saveimg(gcf,ResultPath,Result_Title,filename,101)
% %%
% figure
% % lick raster for all trials 
% [binnedArray, bins] = timestampsToBinned(locs_Lick, eventTimes, binSize, calcWindow);
% % lick raster for only s plus trials
% binnedArray_SP = binnedArray;
% binnedArray_SP(~splus_index,:) = 0;
% [tr,b] = find(binnedArray_SP);
% [rasterX,yy] = rasterize(bins(b));
% rasterY = yy+reshape([zeros(size(tr'));tr';zeros(size(tr'))],1,length(tr)*3); % note from XFu. The original code duplicate the dots to make the ticks longer, which I believe is cheating. 
% rasterX(rasterY==0) = [];rasterY(rasterY==0) = [];% remove zeros
% plot(rasterX, rasterY,'r.')
%% Sanity check: show raster plot of licks, response window
figure
subplot 121
% lick raster for all trials 
[binnedArray, bins] = timestampsToBinned(locs_Lick, eventTimes, binSize, calcWindow);
% lick raster for only s plus trials
binnedArray_SP = binnedArray;
binnedArray_SP(~splus_index,:) = 0;
[tr,b] = find(binnedArray_SP);
[rasterX,yy] = rasterize(bins(b));
rasterY = yy+reshape([zeros(size(tr'));tr';zeros(size(tr'))],1,length(tr)*3); % note from XFu. The original code duplicate the dots to make the ticks longer, which I believe is cheating. 
rasterX(rasterY==0) = [];rasterY(rasterY==0) = [];% remove zeros
plot(rasterX, rasterY,'r.')
hold on 
xline(t_FV_Water_s_avg)
yl1 = ylim;
% hold on
subplot 122
% same thing for sminus
binnedArray_SM = binnedArray;
binnedArray_SM(splus_index,:) = 0;
[tr,b] = find(binnedArray_SM);
[rasterX,yy] = rasterize(bins(b));
rasterY = yy+reshape([zeros(size(tr'));tr';zeros(size(tr'))],1,length(tr)*3); % note from XFu. The original code duplicate the dots to make the ticks longer, which I believe is cheating. 
rasterX(rasterY==0) = [];rasterY(rasterY==0) = [];% remove zeros
plot(rasterX, rasterY,'b.')
hold on 
xline(t_FV_Water_s_avg)
xlim(calcWindow)
% add legend
L1 = plot(nan, nan, 'r.');
L2 = plot(nan, nan, 'b.');
legend([L1, L2], {'S+','S-'})


yl2 = ylim;
% Find leftmost xLeft
ylow = min([yl1(1), yl2(1)]);
% Find rightmost xRight
yhigh = max([yl1(2), yl2(2)]);
subplot 121
ylim([ylow, yhigh]);
subplot 122
ylim([ylow, yhigh]);
title([training_day ' ' pair ' Lick Rasters'])
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
%% calculate the Dprime every 10 trials 
figure
Response = extractfield(Behaviour_Info , 'Response');
HIT = strcmp('HIT',Response);
FA = strcmp('FA',Response);
for session = 1:length(Behaviour_Info)/10
    n = 1:10;
    n10 = (session-1)*10+n;
    s_HIT = sum(HIT(n10));
    s_FA = sum(FA(n10));
    s_Plus = sum(splus_index(n10));
    s_Minus = 10-s_Plus;
    dP = Dprime_Loglinear_norm50(s_HIT,s_FA, s_Plus, s_Minus);
    plot(session,dP,'kx')
    hold on
end
xlabel('x10 trials')
ylabel('D prime')
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