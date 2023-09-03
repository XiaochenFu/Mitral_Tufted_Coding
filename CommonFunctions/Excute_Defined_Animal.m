addpath('C:\Users\yycxx\OneDrive - OIST\Mitral_Tufted_Coding\CommonFunctions')
cd('C:\Users\yycxx\OneDrive - OIST\Mitral_Tufted_Coding\CommonFunctions')
Batch_IDs = {'Round1'};%,'Third','Forth','Fifth'};
% Batch_IDs = {'Second', 'Third','Forth','Fifth'};
% Mice_ID
mydir  = pwd;
idcs   = strfind(mydir,'\');
home_path = mydir(1:idcs(end)-1);
for kkk = 1:length(Batch_IDs)
    Batch_ID = Batch_IDs{kkk};
    % enter the batch folder
    batch_path = fullfile(home_path,[Batch_ID]);
    cd(batch_path )
    D = dir;
    for jjj = 1:length(D)
        currD = D(jjj).name;
        if isfolder(currD)
            if contains(currD,["TbxAi32_86"])
                cd([currD])
                cd('Code')
                %                 Extract_Sniff_Pattern_Individual_Mouse
                %                 Transfer_Learning_Curve_fig
                %                 Extract_learning_Curve_from_Fig
                %
                Restriction = 'Day1';


                GNG_Outside_Defined_Animal
                Restriction = 'Day1';
                GNG_Accuracy_Outside_Defined_Animal
                %                 clearvars -except Batch_IDs home_path kkk Batch_ID batch_path D jjj

%                                 close all
                
            end
            cd(batch_path)
        end
    end
    cd(home_path)
end