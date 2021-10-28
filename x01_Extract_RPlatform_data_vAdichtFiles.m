%% Function to obtain raw data recorded during task with the R-Platform and synchronized EMG with robot data %%
% Extraction of the data from .tdms files (files saved by cRIO-9030)
% Extraction of the data from .c3d files (files saved by VICON nexus)
% Extraction of the LED from .avi files (files saved by SIMI)
% Synchronization with the trigger channel and save of the LED on

clear all
close all
Path = cd;

fS_robot = 100; %sampling frequency
fS_VICON = 2000; % sampling frequency EMG
fS_KIN = 50; %sampling frequency SIMI
ratio = fS_VICON/fS_robot;
cd ('C:\R-Platform\DATA\2021_02_Robot Stim group\Robot Data')

%% .tdms file save struct
[fname_robot, pname_robot] = uigetfile({'*.tdms','Robot Files (*.tdms)'}, 'Choose Files:','MultiSelect', 'on'); % Select File
filenameR = [pname_robot,fname_robot];
cd([Path,'/v2p5'])
addpath(genpath(Path));
Data = TDMS_getStruct(filenameR);
Data.Recorded_Data.fS_robot = fS_robot; %record at 100 Hz
cd(Path)
en = str2num(filenameR(end-6:end-5));
%% .tdms info file struct
cd('C:\R-Platform\DATA\v2p5')
Info = TDMS_getStruct([filenameR(1:end-5),'_info.tdms']);
cd([Path,'\SubFunctions'])
info = GetInfo_TDMS(Info);
Data.info = info;

%% .c3d file save EMG and trigger
cd([Path,'\adinstruments_convert_adichtFiles'])
%c3d = c3dserver; % Start C3D Server 
%setappdata(0,'UseNativeSystemDialogs',false) % Hack for many Files
[fnameV, pnameV] = uigetfile({'*.adicht','Vicon Files (*.adicht)'}, 'Choose Files:','MultiSelect', 'on'); % Open Files with GUI Window
filename = [pnameV, fnameV];
%pRet = c3d.Open(filename,3); % Load File in C3D Server
f = adi.readFile(filename);

%% Get info .c3d file
name_EMG = {'Biceps', 'Triceps', 'ED', 'trigger'};
emg = [];

for i_emg = 1:2
    channel = f.getChannelByName(name_EMG{i_emg});
    raw_ch_data = channel.getData(en); % get data of the en recording
    emg = [emg, raw_ch_data];
end

channel = f.getChannelByName(name_EMG{end});
trigger = channel.getData(en);

%% .Find LED on kinematic
cd([Path,'\SubFunctions'])
[start_LED, time, Pix2cm] = Extraction_LED_called;
Data.SIMI.trig = start_LED;
Data.SIMI.fS_KIN = fS_KIN;
Data.SIMI.duration = time;
Data.SIMI.Px2cm = Pix2cm;

%%
cd(Path)
% plot of Vicon data for checking
figure;
subplot(4,1,1)
plot(emg(:,1))
ylabel('Biceps')
subplot(4,1,2)
plot(emg(:,2))
ylabel('Triceps')
subplot(4,1,3)
%plot(emg(:,3))
ylabel('Digito')
subplot(4,1,4)
plot(trigger)
ylabel('Trigger')

%% Synch EMG-Robot
cd([Path,'\SubFunctions'])
[start_EMG] = synch_EMGVicon_Robot(trigger);
emg(1:start_EMG-1,:)=[]; % delete starting part
cd(Path)

Data.VICON.EMG = emg;
Data.VICON.trig = start_EMG;
Data.VICON.fS_EMG = fS_VICON;
Data.VICON.name = filename;

figure; % figure, to chek synchro
hold on
plot(emg(1:ratio:end,2)*10)
plot(emg(1:ratio:end,1)*10)
plot(Data.Recorded_Data.Fz.data)
plot(trigger(start_EMG:ratio:end))
legend('Triceps','Biceps','Force z','trigger')

%% Save struct with data

save([filenameR(1:end-5),'.mat'],'Data');

