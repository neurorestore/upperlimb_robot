%% Analysis of the EMG around the force peaks
% 1. Frequency analysis of EMG around force peaks or start movement (Power 
% Spectrum and Spectrogram)
% 2. Graph Force vs EMG

clear all
close all
Path = cd;

INT = 0.2; %seconds to integrate the signal and plot them
in_case = 1; % switch for aligning (1-> force peaks, 2-> start movement)

%% Load
cd ('C:\R-Platform\DATA\2021_02_Robot Stim group\Robot Data')
[filename, pathname] = uigetfile({'*.mat','Matlab AnalysisPeaks (*.mat)'}, 'Choose Files:','MultiSelect', 'on'); % Select File
% Load file data synchronized
load([pathname,filename]);

fSrobot = Data.Recorded_Data.fS_robot;
fSEMG = Data.VICON.fS_EMG;
cd([Path,'\SubFunctions'])

%% EMG around onset of force peaks in zeta direction
file = [filename(1:end-5),'_Fz'];
if strcmp (Data.info.Status,'Pull_active')
    [Spect,rel_max,MeanEnvelope] = PlotEMGaroundonset (Data.Recorded_Data.Analysis.Fzpeaks_contra,Data.SIMI.start,fSrobot,Data.VICON.EMG, Data.VICON.EMGfilt, fSEMG, Data.Recorded_Data.Fz.data,...
    pathname,file, in_case,Data.info.Status);
else
[Spect,rel_max,MeanEnvelope] = PlotEMGaroundonset (Data.Recorded_Data.Analysis.Fzpeaks,Data.SIMI.start,fSrobot,Data.VICON.EMG, Data.VICON.EMGfilt, fSEMG, Data.Recorded_Data.Fz.data,...
    pathname,file, in_case,Data.info.Status);
end
Data.VICON.spect = Spect;
Data.VICON.EnvOnsetF.rel_max = rel_max;
Data.VICON.EnvOnsetF.MeanEnvelope = MeanEnvelope;


%% EMG vs Force graph
cd([Path,'\SubFunctions'])
[R,Rt,Rd,Area_cy] = PlotEMGvsF(Data.Recorded_Data.Fz.data,Data.VICON.EMGfilt,Data.good_trials,fSrobot,fSEMG,INT,rel_max(:,2,:),[pathname filename]);

nEMG = size(Data.VICON.EMGfilt,2);
Data.VICON.Analysis.Corr = [R(2,1),Rt(2,1),Rd']; % correlation coefficient biceps, triceps, DE
Data.VICON.Analysis.Area_cy = Area_cy(:,2:nEMG+1); % mean integration of the EMG activity during every cycle

Data.Recorded_Data.Analysis.Area_cy = Area_cy(:,1); % mean integration of the force z activity during every cycle

save([pathname,filename(1:end-5),'_x04.mat'],'Data');

% file = [filename(1:end-5),'_Fx'];
% PlotEMGaroundonset (Data.Recorded_Data.Analysis.Fxpeaks,Data.Recorded_Data.fS_robot,Data.VICON.EMG, Data.VICON.EMGfilt, Data.VICON.fS_EMG, pathname,file);
% cd([Path,'\SubFunctions'])

% d=1;
%     for n = 1:20%nburst
%         if burst_duration(n)>0.5
%             figure
%             plot(EMG(st_burst(n):en_burst(n),nEMG))
%             [pxx,f] = pwelch(EMG(st_burst(n):en_burst(n),nEMG),100,50,[],fS);
%             %plot(f,10*log10(pxx))
%             %xlabel('Frequency (Hz)')
%             %ylabel('Magnitude (dB)')
%             figure
%             medfreq(pxx,f);
%             p(d) = medfreq(pxx,f);
%             d = d+1;
%         end
%     end
%    mean(p)