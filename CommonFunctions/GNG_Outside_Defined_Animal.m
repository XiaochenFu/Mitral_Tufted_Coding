addpath(genpath('C:\Users\yycxx\OneDrive - OIST\Ephys_Code\Library\spikes-master'));
addpath('C:\Users\yycxx\Dropbox (OIST)\Fukunaga_Lab_Joined\Code\Basic_Settings')
addpath('C:\Users\yycxx\Dropbox (OIST)\Fukunaga_Lab_Joined\Code\Useful_Functions')
addpath('C:\Users\yycxx\Dropbox (OIST)\Fukunaga_Lab_Joined\Code\Useful_Functions\folderfile')
startup
close all
mydir  = pwd;
idcs   = strfind(mydir,'\');
newdir = mydir(1:idcs(end)-1);
DataPath = fullfile(newdir,'Recording');
ResultPath = fullfile(newdir,'Behaviour_Preprocess');
searchPath = [DataPath ,'\**\*.csv']; % Search in folder and subfolders for  *.csv
FileNames      = dir(searchPath); % Find all .csv files

processedFiles = {}; % Initialize empty cell array to keep track of processed files

for i = 1:length(FileNames)
    FileName = FileNames(i).name;

    % If the file has already been processed, skip this iteration
    if ismember(FileName, processedFiles)
        continue;
    end

    if exist('Restriction','var')
        filenames = get_defined_file_names(searchPath,Restriction);
        filteredCell = myCell(~contains(myCell, 'abort','IgnoreCase',true));

        if any(contains(filenames,FileName))
            close all;
            

            % Extract day from the filename
            dayPattern = 'Day\d+'; % Regular expression to match 'Day' followed by numbers
            dayMatch = regexp(FileName, dayPattern, 'match');

            if ~isempty(dayMatch)
                currentDay = dayMatch{1};

                % Find all files corresponding to the same day
                sameDayFiles = filenames(contains(filenames, currentDay));
                % If there are multiple files from the same day, concatenate their data
                concatenatedData = [];
                if length(sameDayFiles) > 1
                    for j = 1:length(sameDayFiles)
                        tempData = readmatrix(fullfile(DataPath, sameDayFiles{j}));
                        concatenatedData = [concatenatedData; tempData];
                        processedFiles{end+1} = sameDayFiles{j}; % Add the file to the processed list
                    end
                    DATA = concatenatedData;
                    GNG_Inside; 
                else
                    DATA = readmatrix(fullfile(DataPath,FileName));
                    GNG_Inside; 
                    processedFiles{end+1} = FileName; % Add the file to the processed list
                end
            else
                processedFiles{end+1} = FileName; % Add the file to the processed list
            end
        end
    end
end


