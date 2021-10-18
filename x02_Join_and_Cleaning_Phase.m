%% Join of multiple data and cleaning of bad trials
% 1. join different recordings of the same day of the same type
% 2. find the good trials


clear all
close all
Path = cd;

LoadcellName = {'Fx','Fy','Fz'};
LoadcellNameM = {'Mx','My','Mz'};

%% Load
cd ('C:\R-Platform\DATA\2021_02_Robot Stim group\Robot Data')
[filename, pathname] = uigetfile({'*.mat','Matlab extracted data (*.mat)'}, 'Choose Files:','MultiSelect', 'on'); % Select File
% Load file data synchronized
load([pathname,filename]);

% Load other files if more recording of the same type in one day
prompt = 'Do you want to upload other files? (y/n)';
selection = input(prompt,'s');
 
if contains(selection,'y')
    cd([Path,'\SubFunctions'])
    Data = Join_multiple_Recordings(Data,selection,pathname); 
    cd(Path)
else
    %timeEMG = size(Data.VICON.EMG,1)/Data.VICON.fS_EMG;
    timeRobot = length(Data.Recorded_Data.t.data)/Data.Recorded_Data.fS_robot;
    timeKIN = Data.SIMI.duration;
    Data.Props.SingleTime = min([timeEMG,timeRobot,timeKIN]);
    % change direction force for right forelimb
    if contains(Data.info.Paw,'right') % check if I have to change also moment!
        Data.Recorded_Data.Fz.data = -Data.Recorded_Data.Fz.data;
        Data.Recorded_Data.Mx.data = -Data.Recorded_Data.Mx.data;
    end
    % remove offset to force signal
    for i_force = 1:size(LoadcellName,2)
        fx = Data.Recorded_Data.(LoadcellName{i_force}).data;
        mfx  = median(fx(Data.Recorded_Data.T_status.data==0));
        Data.Recorded_Data.(LoadcellName{i_force}).data = fx-mfx;
    end
end

%% Clean Data (eliminate last part for longer force), select only good trials
% save number of good trials
cd([Path,'\SubFunctions'])
Data = Cleaning_Trials(Data);

%% New filename
% new name to save the data including all files
if selection == 'y'
    num_files = size (Data.Props.name,1);
    filename_new = Data.Props.name(1,:);
    for i=1:num_files-1
        add_file = Data.Props.name(i+1,end-5:end);
        filename_new =[filename_new,'_',add_file];
    end    
else
    filename_new = filename(1:end-4);
end

%% Save data

save([pathname,filename_new,'x02.mat'],'Data');
cd (Path)