addpath(genpath('C:\Users\yycxx\OneDrive - OIST\Ephys_Code\spikes-master'));
addpath('C:\Users\yycxx\Dropbox (OIST)\Fukunaga_Lab_Joined\Code\Basic_Settings')
addpath('C:\Users\yycxx\Dropbox (OIST)\Fukunaga_Lab_Joined\Code\Useful_Functions')
startup
mydir  = pwd;
idcs   = strfind(mydir,'\');
newdir = mydir(1:idcs(end)-1);
DataPath = fullfile(newdir,'Behaviour_Preprocess');
searchPath = [DataPath ,'\**\*.mat']; % Search in folder and subfolders for  *.csv
FileNames      = dir(searchPath); % Find all .csv files
for i = 1:length(FileNames)
    FileName = FileNames(i).name;
    %     if contains(FileName, 'Day10')
    if exist('Restriction','var')
        filenames = get_defined_file_names(searchPath,Restriction);
        if any(contains(filenames,FileName))
            close all;GNG_Accuracy_Inside; %clearvars -except i FileNames;
        end
    end
    % end
end