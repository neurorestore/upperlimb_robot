%% x02bis Add movement 2d (CAM4)
% 1. Load all the files extracted with DLC (csv)
% 2. Synchronize them thanks to the extraction of the LED performed on x01
% 3. Fix problems in the coordinates when the likelihood is too low
% 4. Plot all data coordinated (and save it)

close all
clear all

% number of points extracted with DLC (in this moment is only the wrist,
% but we can imagine to add other points)
points_DLC = 1; 

Path = cd;

% Load files
cd ('C:\R-Platform\DATA\2021_02_Robot Stim group\Robot Data')
[filename, pathname] = uigetfile({'*.mat','Data struct x02(*.mat)'}, 'Choose Files:','MultiSelect', 'on'); % Select struct matlab x02
load([pathname,filename]);
%cd('C:\R-Platform\DATA\2019_02_Group1\SIMI results\DAY02')
n_files = length(Data.Props.SingleTime); % check how many files are joined together
pos_tot = [];
for i_file = 1 : n_files
    
    [filename_xls, pathname_xls] = uigetfile({'*.csv','Position of DLC (*.csv)';'*.*', 'All Files (*.*)'}, 'Choose Files:','MultiSelect', 'on'); % Select File csv coming from DLC
        pos = csvread([pathname_xls,filename_xls],3,0);
    
    name = readtable([pathname_xls,filename_xls],'ReadVariableNames',false);
    name_matrix = table2array(name(1:3,:));
    clear name

    fS_robot = Data.Recorded_Data.fS_robot;
    fS = Data.SIMI.fS_KIN;
    ratio = fS_robot/fS;
    trigger = Data.SIMI.trig(i_file)*fS;
    fine = round(Data.Props.SingleTime(i_file)*fS);
    pos_syn = pos(trigger:trigger+fine-1,:);
    x = find(contains (name_matrix(3,:),'x'));
    y = find(contains (name_matrix(3,:),'y'));
    like = find(contains (name_matrix(3,:),'likelihood'));
    
    %change direction movement for right forelimb
    if strcmp(Data.info.Paw,'right')
        for n_points = 1 : points_DLC
            pos_syn(:,x(n_points)) = -pos_syn(:,x(n_points));
        end
    end
    
    cd ([Path,'\SubFunctions'])
    %fix the trajectory when the likelihood is very low
    for n_points = 1 : points_DLC
        likelihood = pos_syn(:,like(n_points));
        pos_syn(:,x(n_points)) = Reconstruct_Trajectory(likelihood, pos_syn(:,x(n_points)),fS);
        pos_syn(:,y(n_points)) = Reconstruct_Trajectory(likelihood, pos_syn(:,y(n_points)),fS);
%         index_low = find(likelihood<0.75);
%         ii = find(diff(index_low)>1.5);
%         pos_syn(index_low,x(n_points)) = NaN;
%         pos_syn(:,x(n_points)) = fillgaps(pos_syn(:,x(n_points)),80,3);%fillmissing(pos_syn(:,x(n_points)),'pchip');
%         pos_syn(index_low,y(n_points)) = NaN;
%         pos_syn(:,y(n_points)) = fillmissing(pos_syn(:,y(n_points)),'pchip');
        % figure to check if start and end of points with low likelihood
        % where well identified
%         figure;plot(likelihood)
%         hold on 
%         scatter(st_l,ones(length(st_l),1))
%         scatter(en_l,ones(length(st_l),1))
    end
    

    %filter of the position signal
    for n_points = 1:points_DLC
        pos_syn(:,x(n_points))= sgolayfilt(pos_syn(:,x(n_points)),3,15);
        pos_syn(:,y(n_points))= sgolayfilt(pos_syn(:,y(n_points)),3,15);
    end
    
    pos_tot = [ pos_tot; pos_syn];
    cd(pathname_xls)
end

Data.SIMI.x = pos_tot(:,x(1));
Data.SIMI.y = pos_tot(:,y(1));

save([pathname,filename],'Data')

%% figure to verify if it is synchronized
figure
hold on
xt = [1/fS:1/fS:size(pos_tot)/fS];
plot(xt,(pos_tot(:,x(1))-nanmean(pos_tot(:,x(1))))/60,'LineWidth',1.5)
plot(xt,pos_tot(:,4),'LineWidth',1.5)
fz = Data.Recorded_Data.Fz.data;
xt = [1/fS_robot:1/fS_robot:length(fz)/fS_robot];
plot(xt,fz)

%% figure all synch
PLOT = 1;
if PLOT
    figure
    hold on
    subplot(3,1,1)
    fz = Data.Recorded_Data.Fz.data;
    xt = [1/fS_robot:1/fS_robot:length(fz)/fS_robot];
    plot(xt,fz,'b')
    ylabel('Force (N)')
    %xlim([0 80])
%     subplot(3,1,2)
%     fS_EMG = Data.VICON.fS_EMG;
%     Bi = Data.VICON.EMG(:,1);
%     xt = [1/fS_EMG:1/fS_EMG:length(Bi)/fS_EMG];
%     plot(xt,Bi,'m')
%     %xlim([0 80])
%     ylabel('EMG Biceps (mV)')
    subplot(3,1,3)
    xt = [1/fS:1/fS:size(pos_tot)/fS];
    plot(xt,(pos_tot(:,x(1))),'g','LineWidth',1.5)
    ylabel('Position')
    %xlim([0 80])
    xlabel('Time (s)')
    savefig([pathname,filename(1:end-4),'_SyncData.fig'])
end
%% Add yellow band to a figure
% f = openfig('C:\R-Platform\DATA\Group1_2019_02\Recordings_robot\DAY10_2019_10_01_8DPI\R04\R042019_10_01x8DPIFile 3_File 5x02_SyncData.fig')
% for i=1:3
% hold on
% subplot(3,1,i); hold on
% for i_ac = 1:17
%         xx = [Data.SIMI.start(i_ac) Data.SIMI.start(i_ac)+1 Data.SIMI.start(i_ac)+1 Data.SIMI.start(i_ac)];
%         yy = [-3 -3 1100 1100];
%         p = patch( xx,yy,[255 237 0]./255,'EdgeColor','none');
%         set(p,'FaceAlpha',0.5)
% end
% end