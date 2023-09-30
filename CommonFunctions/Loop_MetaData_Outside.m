% Calculate the Meta data, d' of the mouse with Loop_dprime_Inside.m. To use, call
% this function in Excute_Defined_Animal.m
addpath('C:\Users\yycxx\Dropbox (OIST)\Fukunaga_Lab_Joined\Code\Useful_Functions\folderfile')
addpath(genpath('C:\Users\yycxx\OneDrive - OIST\Ephys_Code\Library\spikes-master'));
addpath('C:\Users\yycxx\Dropbox (OIST)\Fukunaga_Lab_Joined\Code\Basic_Settings')
addpath('C:\Users\yycxx\Dropbox (OIST)\Fukunaga_Lab_Joined\Code\Useful_Functions')
addpath('C:\Users\yycxx\Dropbox (OIST)\Fukunaga_Lab_Joined\Code\Useful_Functions\folderfile\natsortfiles')
% addpath('C:\Users\yycxx\OneDrive - OIST\Thesis\Behaviour\Group_Analyses\Summary_Stimuli')
set(0,'DefaultAxesFontSize',6)
startup

csv_filename = 'experiment_metadata.xlsx';
output_path = 'C:\Users\yycxx\OneDrive - OIST\Mitral_Tufted_Coding\Meta_Data';



mydir  = pwd;
idcs   = strfind(mydir,'\');
newdir = mydir(1:idcs(end)-1);
DataPath = fullfile(newdir,'Behaviour_Preprocess');
searchPath = [DataPath ,'\**\*.mat']; % Search in folder and subfolders for  *.csv
FileNames      = dir(searchPath); % Find all .csv files
FileNames = natsortfiles(FileNames);
outputTable = table('Size', [0 11],...
    'VariableTypes', {'string', 'string', 'string','string', 'string', 'string','string', 'string', 'string','string', 'string'},...
    'VariableNames', ...
    {'Session','Stimuli	Date','Trials',...
    'Accuracy %','Overall_dprime',...
    'highest_dprime_logniear50','trials_to_criteria_loglinear50',...
    'highest_dprime_logniear','trials_to_criteria_loglinear',...
    'highest_dprime_2N','trials_to_criteria_2N'});

for i = 1:length(FileNames)
    FileName = FileNames(i).name;
    %     if contains(FileName, 'Day10')
    if exist('Restriction','var')
        filenames = get_defined_file_names(searchPath,Restriction);
        if any(contains(filenames,FileName))
            Loop_MetaData_Inside ; %clearvars -except i FileNames;
            newRow = {training_day, pair, length(Behaviour_Info),...
                 round(All_Accuracy*100 * 10) / 10, All_dprime,...
                dprimes_Loglinear50_max, num_trial_Loglinear50,...
                dprimes_Loglinear_max, num_trial_Loglinear,...
                dprimes_2N_max, num_trial_2N,...
                };
            
            outputTable = [outputTable; newRow];
            writetable(outputTable, fullfile(output_path,csv_filename),'Sheet',currD);
        end
    end
end
