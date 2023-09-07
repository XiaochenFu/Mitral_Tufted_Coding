function Result_Title = Convert_Behaviour_FileName(FileName)
idcs   = strfind(FileName,'_');
Training_Date = FileName((idcs(1)-2):(idcs(1)-1));  % The number of training day is the part after the first unserscore
Training_Date =  erase( Training_Date , 'y' );
Stimuli0 = FileName(idcs(1)+1:idcs(2)-1);   % The stimuli is the part after the first unserscore
% for odour, we don't need to check the light command (and there should be no command if we run the dummy.vi)
if contains(Stimuli0,'Odour')
    Stimuli_Type = 'Odour';
    Stimuli = Stimuli0;
else
    Stimuli_Type = 'Light';
    switch Stimuli0
        case 'Detection3vs0'
            Stimuli = 'Detetion';
    end
    
%     if contains(FileName,'vs0') || contains(FileName,'v0')
%         Stimuli = 'Detection';% should be odour pair or light(Light detection intensity wi/wo odour, Light discrimination wi/wo odour)
%     else
%         Stimuli = 'Discrimination';
%         %         phase_name = asManyOfPattern(digitsPattern(1))+ 'ms'+'vs'+ asManyOfPattern(digitsPattern(1));
%         %         phase_pattern = namedPattern(phase_name);
%         %         extract(Stimuli0,phase_pattern)
%         if contains(FileName,'Delay')
%             Stimuli = 'Phase';
%         end
%         if contains(FileName,'Reverse')
%             Stimuli = 'Reverse';
%         end
% 
%     end
%     if contains(FileName,'mW')
%         idcs   = strfind(FileName,'mW');
%         Stimuli = [Stimuli FileName(idcs(1)-3:idcs(1)+1)];
%         Stimuli =  erase( Stimuli , '_' );
%     else
%         Stimuli =  'noLight';
%     end
%     if contains(FileName,'MV')
%         Stimuli =  strcat(Stimuli, 'MV');
%     end
%     if contains(FileName,'MT')
%         Stimuli =  strcat(Stimuli, 'MT');
%     end
end
Result_Title = sprintf('Day%s_%s',Training_Date, Stimuli)