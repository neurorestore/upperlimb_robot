%% Analysis of the peaks of load cell signal and of EMG recordings
% 1. join different recordings of the same day of the same type
% 2. find the good trials
% 3. EMG burst extraction (in all signal)
% 4. Force peaks extraction (in all the 6 channels)


clear all
close all
Path = cd;

LoadcellName = {'Fx','Fy','Fz'};
LoadcellNameM = {'Mx','My','Mz'};

%% Load
cd ('C:\R-Platform\DATA\2021_02_Robot Stim group\Robot Data')
[filename, pathname] = uigetfile({'*.mat','Matlab struct x02 (*.mat)'}, 'Choose Files:','MultiSelect', 'on'); % Select File
% Load file data synchronized
load([pathname,filename]);
filename_new = filename(1:end-7);

%% EMG extraction of bursts
cd([Path,'\SubFunctions'])
% extraction of all burst 
[EMGenv,burst] = CalculateEMGparam(Data.VICON.EMG,Data.VICON.fS_EMG);
Data.VICON.EMGfilt = EMGenv;
Data.VICON.burst = burst;

% Delete burst that are outside good trials
[Data] = CleanBurst(Data);

[EMGav_tot,EMGav_pull,EMGav_push] = CoactivationEMG (Data.VICON.EMG,Data.VICON.fS_EMG,Data.good_trials,Data.Recorded_Data.T_status.data,pathname,filename_new,Data.Recorded_Data.fS_robot);
Data.VICON.Coactivation.EMGav_tot = EMGav_tot;
Data.VICON.Coactivation.EMGav_pull = EMGav_pull;
Data.VICON.Coactivation.EMGav_push = EMGav_push;

%% Force signal analysis (find force peaks and all correlated information)

% THIS function has been moved before in x02, to do it file by file
% the force signal is opposite in the right paw, so we invert the signal to
% have the same data for everyone
% if contains(Data.info.Paw,'right') % check if I have to change also moment!
%     Data.Recorded_Data.Fz.data = -Data.Recorded_Data.Fz.data;
%     Data.Recorded_Data.Mx.data = -Data.Recorded_Data.Mx.data;
% end

for i_force = 1:size(LoadcellName,2)
    if strcmp(Data.info.Status,'active') % for task where the animal is completly passive (info.Status == active) we remove the force generate by the robot
        [Data.Recorded_Data.(LoadcellName{:,i_force}).data] = Remove_artificialForce(Data.Recorded_Data.(LoadcellName{:,i_force}).data,Data.Recorded_Data.T_status.data, Data.good_trials,Data.Recorded_Data.fS_robot);
    end
    [peaks,peaks_contra,names] = CalculatePEAKS (Data.Recorded_Data.(LoadcellName{:,i_force}).data,Data.Recorded_Data.fS_robot,Data.Recorded_Data.T_status.data);
    new_fields{1} = [LoadcellName{:,i_force},'peaks'];
    new_fields{2} = [LoadcellName{:,i_force},'peaks_contra'];
    Data.Recorded_Data.Analysis.(new_fields{:,1})= peaks;
    Data.Recorded_Data.Analysis.(new_fields{:,2})= peaks_contra;
    [peaks,peaks_contra,names] = CalculatePEAKS_moments (Data.Recorded_Data.(LoadcellNameM{:,i_force}).data,Data.Recorded_Data.fS_robot);
    new_fields{1} = [LoadcellNameM{:,i_force},'peaks'];
    new_fields{2} = [LoadcellNameM{:,i_force},'peaks_contra'];
    Data.Recorded_Data.Analysis.(new_fields{:,1})= peaks;
    Data.Recorded_Data.Analysis.(new_fields{:,2})= peaks_contra;
end
new_fields{3} = ('name_peaks_matrix');
Data.Recorded_Data.Analysis.(new_fields{:,3})= names;

%% Eliminate peaks in not good trials and add information about number of
% cicle and status of the task
cd([Path,'\SubFunctions'])
[Analysis] = InsideGoodTrials (Data.Recorded_Data.Analysis, Data.Recorded_Data.cicles.data, Data.good_trials, Data.Recorded_Data.T_status.data);
Data.Recorded_Data.Analysis = Analysis;

%% SIMI extraction of movements (find kinematic parameters)
% check if kinematic info have been upload
if isfield(Data.SIMI,'x')
    Data = CalculateKINparam(Data);
else
    fprintf('Kinematic Data are not present in the selected struct:\n Parameters can not be extracted\n')
    Data.Recorded_Data.Analysis.Fzpeaks(:,16) = NaN;
end

%% Save index start movements
% data from the kinematic analysis for Pull-active task
% data from movement of spindle-drive motor for Passive task

if strcmp(Data.info.Status,'Push_active') || strcmp(Data.info.Status,'Pull_active')
    if isfield(Data.SIMI,'pks')
        Data.SIMI.start = Data.SIMI.pks{1,3}(Data.SIMI.pks{1,6}==1); %only movements pks associated with a force peaks
    else
        Data.SIMI.start = NaN;
    end
elseif strcmp(Data.info.Status,'active')
    status = Data.Recorded_Data.T_status.data;
    inin = find(status==2);
    index = inin(find (diff(inin)>1.5)+1);
    index = [inin(1),index];
    Data.SIMI.start = index/Data.Recorded_Data.fS_robot;    
end

%% Save data
save([pathname,filename_new,'_AnalysisPeaks.mat'],'Data');



