% Function to check if burst are inside good trials or not, to keep only
% the good ones

function [Data] = CleanBurst(Data)

good_trials = Data.good_trials;
cycles = Data.Recorded_Data.cicles.data;
fS = Data.VICON.fS_EMG;
cycles = round(resample(cycles,fS,Data.Recorded_Data.fS_robot));

% save name of parameters extracted for the burst because it has not been
% saved in the previous code
name_burst = {'start','mean','max','duration','AUC'};
Data.VICON.name_burst = name_burst;
burst = Data.VICON.burst;

for i_EMG = 1:size(burst,1) % scroll all the recorded muscles
    bb = burst{i_EMG,1};
    start = bb(:,1)*fS;
    duration = bb(:,4)*fS;
    for i_bb = length(start):-1:1
        single_cy = unique(cycles(start(i_bb):start(i_bb)+duration(i_bb))); % find the cycle in which the burst is active
        if sum(ismember(single_cy,good_trials))<length(single_cy) % at least one part of the burst is outside the trial
            bb(i_bb,:)=[];
        end
    end
    burst{i_EMG,1} = bb;
end

Data.VICON.burst = burst;

end