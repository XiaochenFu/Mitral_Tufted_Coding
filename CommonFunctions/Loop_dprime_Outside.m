addpath(genpath('C:\Users\yycxx\OneDrive - OIST\Ephys_Code\Library\spikes-master'));
addpath('C:\Users\yycxx\Dropbox (OIST)\Fukunaga_Lab_Joined\Code\Basic_Settings')
addpath('C:\Users\yycxx\Dropbox (OIST)\Fukunaga_Lab_Joined\Code\Useful_Functions')
addpath('C:\Users\yycxx\Dropbox (OIST)\Fukunaga_Lab_Joined\Code\Useful_Functions\folderfile\natsortfiles')
addpath('C:\Users\yycxx\OneDrive - OIST\Thesis\Behaviour\Group_Analyses\Summary_Stimuli')
set(0,'DefaultAxesFontSize',6)
startup
mydir  = pwd;
idcs   = strfind(mydir,'\');
newdir = mydir(1:idcs(end)-1);
DataPath = fullfile(newdir,'Behaviour_Preprocess');
searchPath = [DataPath ,'\**\*.mat']; % Search in folder and subfolders for  *.csv
FileNames      = dir(searchPath); % Find all .csv files
FileNames = natsortfiles(FileNames);
figure('units','normalized','outerposition',[0 0 1 1])
f = tiledlayout(5,6);
for i = 1:length(FileNames)
    FileName = FileNames(i).name;
    %     if contains(FileName, 'Day10')
    if exist('Restriction','var')
        filenames = get_defined_file_names(searchPath,Restriction);
        if any(contains(filenames,FileName))
            Loop_dprime_Inside ; %clearvars -except i FileNames;
        end
    end
    % end
end
title(f,replace(currD,'_',' '),'FontSize',20)
xlabel(f,'x10 trials','FontSize',18)
ylabel(f,'D prime','FontSize',18)