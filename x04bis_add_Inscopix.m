%% Function to synchronize inscopix data with robot data
% Run once for every animal
% IMPORTANT -> Things to do before running this code:
% 1. Run Matlab Codes x01, x03, x04;
% 2. Extract single units activity from inscopix files (inscopix program)

%%REMEMBER ADD deltaF/F for data not analysed on the inscopix software!!

close all
clear all

% Select folder for file mat analysis
data_robot = 'C:\R-Platform\DATA\2019_02_Group1\Recordings_robot\';
% Select folder for inscopix
data_inscopix = 'C:\R-Platform\DATA\2019_02_Group1\SUs-software inscopix\';
%add = 'RegCells'; %folder with the matrix output of GUI CellReg

% Select animal
Rat = 'R04';

% Select days to analyse
DAYS = {'DAY03','DAY04','DAY05','DAY07','DAY09_4DPI','DAY10_8DPI','DAY11_14DPI','DAY12_21DPI','DAY13_35DPI'};%,'DAY04','DAY10_8DPI'};

% Select type of task (half-Push_active, active-active)
Task = {'half','active'};
TaskR = {'Push_active','active'};

% FILT=1 if filtering calcium signals
FILT = 1;

fSins = 20; % recordings 20 Hz 
path = cd;
%% Load data
% INSCOPIX
% initialize a variable where all the cells info are collected:
% for each day: 1. name of the unit (by inscopix software)
% 2. activity of units along time
% 3. position (x,y) of the centroid
% 4. size of the unit

% ATTENTION: this matrix is the same of the one in CellRegistration,
% otherwise numbers (of the cells) in the reg_file don't match! 

allCells = cell(length(DAYS),size(Task,2),8);

for i_gg = 1:length(DAYS)
    for i_tas = 1:size(Task,2)        
    pathname = [data_inscopix,Rat,'\'];
    filename = [DAYS{i_gg},'_',Task{i_tas},'.csv'];
    % Read data from the file of the inscopix program .csv
    inscopix = readtable([pathname,filename],'ReadVariableNames',false);
    name_cells = table2cell(inscopix(1,2:end));
    status_cells = table2cell(inscopix(2,2:end));
    cells = csvread([pathname,filename],2,1);
    dim = size(cells);
    time_cells = csvread([pathname,filename],2,0,[2,0,dim(1),0]);
    clear inscopix
    filename = [DAYS{i_gg},'_',Task{i_tas},'-props.csv'];
    pos = csvread([pathname,filename],1,5,[1,5,dim(2),6])';
    size_cells = csvread([pathname,filename],1,8,[1,8,dim(2),8])';
    % eliminate rejected cells
    %b = ones(10,1)/10;
%    dfil = designfilt('lowpassfir','PassbandFrequency',2, ...
%         'HalfPowerFrequency',0.15,'DesignMethod','butter');
    dfil = designfilt('lowpassiir','FilterOrder',12, ...
        'HalfPowerFrequency',0.15,'DesignMethod','butter');

    for i_c = dim(2):-1:1 % start with the last one, to make sure that the index doesn't change
        if contains (status_cells{i_c},'rejected')
            name_cells(i_c)=[];
            cells(:,i_c)=[];
            pos(:,i_c) = [];
            size_cells(:,i_c) = [];
        else
            if FILT % FILTERING
                if sum(isnan(cells(:,i_c)))>0.5
                    cells(:,i_c) = fillmissing(cells(:,i_c),'movmean',4);
                end
                %cells(:,i_c) = filtfilt(b,1,cells(:,i_c));
                cells(:,i_c) = filtfilt(dfil,cells(:,i_c));
            end
        end
    end
    % Load event inscopix file (SNR, Frequency rate, time events
%     filename = [DAYS{i_gg},'_',Task{i_tas},'-Events.csv'];
%     inscopix_events = readtable([pathname,filename],'ReadVariableNames',false);
%     name_events = table2cell(inscopix_events(2:end,2));
%     dim_ev = size(cells,2);
%     time_events = str2double(table2array(inscopix_events(2:end,1)));
%     clear inscopix_events
%     filename = [filename(1:end-4),'-props.csv'];
%     SNR = csvread([pathname,filename],1,1,[1,1,dim_ev,1]);
%     Frate = csvread([pathname,filename],1,2,[1,2,dim_ev,2]);
%     % Create a matrix where there are all the events divided for cell
%     d = 0;
%     for i_c = 1:dim_ev % find the cell with more events
%         if d<sum(contains(name_events,name_cells(i_c)))
%             d = sum(contains(name_events,name_cells(i_c)));
%         end
%     end
%     t_events = NaN(d,dim_ev);
%     for i_c = 1:dim_ev % fill the matrix
%         index = contains(name_events,name_cells(i_c));
%         t_events(1:sum(index),i_c)= time_events(index ==1);
%     end
%     clear inscopix_events
    % save data in a big matrix
    allCells{i_gg,i_tas,1} = name_cells;
    allCells{i_gg,i_tas,2} = cells;
    allCells{i_gg,i_tas,3} = pos;
    allCells{i_gg,i_tas,4} = size_cells;
    allCells{i_gg,i_tas,5} = time_cells;
%     allCells{i_gg,i_tas,6} = t_events;
%     allCells{i_gg,i_tas,7} = SNR;
%     allCells{i_gg,i_tas,8} = Frate;
    end
end

%% ROBOT 
% load files AnalysisPeak_x04.mat
ListDAY = dir(data_robot);
allData = cell(length(DAYS),size(Task,2));
for i_gg = 1: length(DAYS)
    indF = NaN;
    for i_d = 3:length(ListDAY) % find the correct folder
        if contains(ListDAY(i_d).name,DAYS{i_gg}(1:5))
            ind = i_d;
        end
    end
    ListFiles = dir ([data_robot,ListDAY(ind).name,'\',Rat]);
    for i_tas = 1:size(Task,2) % cycle for all the tasks, to save data
    for i_f = 3:length(ListFiles) % find the correct file
        if contains(ListFiles(i_f).name,'Peak_x04.mat')
            load ([data_robot,ListDAY(ind).name,'\',Rat,'\',ListFiles(i_f).name]);
            if strcmp(TaskR{i_tas},Data.info.Status)
                indF = i_f;
            end
            clear Data
        end
    end
    allData{i_gg,i_tas} = load ([data_robot,ListDAY(ind).name,'\',Rat,'\',ListFiles(indF).name]);
    % check if the number of files are corrected (i.e. same number of files
    % in the inscopix and in the matrix robot)
    n_datafile = size(allData{i_gg,i_tas}.Data.Props.name,1);
    st_files_ins = 1;
    for i_switch = 2: length(allCells{i_gg,i_tas,5})
        if allCells{i_gg,i_tas,5}(i_switch)-allCells{i_gg,i_tas,5}(i_switch-1)>1.5 % not consecutive times
            st_files_ins = [st_files_ins,i_switch];
        end
    end
    if n_datafile~=length(st_files_ins)
        disp('Error, number of files in the inscopix data is not the same than in the recorded matrix')
        return
    end
    % if the same number of files are present for all data, we assumed that
    % the upload data are correct
        
    % !!!ATTENTION if we decide to use the events extracted by the inscopix
    % files it is necessary to modify their time removing the big step
    % found it in the time_cell column of allCells!!!
    
    %% synchronization
    for i_sw = 1:n_datafile
        endFile = round(sum(allData{i_gg,i_tas}.Data.Props.SingleTime(1:i_sw))*fSins);
        if i_sw==n_datafile
            endInsc = length(allCells{i_gg,i_tas,5});
        else
            endInsc = st_files_ins(i_sw+1)-1;
        end
        allCells{i_gg,i_tas,2}(endFile+1:endInsc,:)=[];
        allCells{i_gg,i_tas,5}(endFile+1:endInsc)=[];        
    end
    %modify time of events according to the new time and after removing the
    %big step that the inscopix software, create when it joins multiple
    %files
    for i_r = 1:size(allCells{i_gg,i_tas,6},1) % scroll every raw of the matrix with event time
        for i_c = 1: size(allCells{i_gg,i_tas,6},2) % scroll every column of the matrix
            tempo = round(allCells{i_gg,i_tas,5}*fSins);
            if ~ isnan(allCells{i_gg,i_tas,6}(i_r,i_c))
                in_tempo = find(round(allCells{i_gg,i_tas,6}(i_r,i_c)*fSins)==tempo);
                if ~isempty(in_tempo)
                    allCells{i_gg,i_tas,6}(i_r,i_c) = in_tempo/fSins;
                else
                    allCells{i_gg,i_tas,6}(i_r,i_c) = NaN;
                end
            end
        end
    end
    figure
    hold on
    fS_robot = allData{i_gg,i_tas}.Data.Recorded_Data.fS_robot;
    plot(allData{i_gg,i_tas}.Data.Recorded_Data.Fz.data(1:fS_robot/fSins:end)*6000)
    plot(allCells{i_gg,i_tas,2})
    end

end

%% Save all data synchronized in a folder

for i_gg = 1:length(DAYS)
    for i_tas = 1:size(Task,2)
        
        %folder to save
        for i_d = 3:length(ListDAY) % find the correct folder
            if contains(ListDAY(i_d).name,DAYS{i_gg}(1:5))
                pathsave = [data_robot,ListDAY(i_d).name,'\',Rat,'\'];
                ind = i_d;
            end
        end
        ListFiles = dir ([data_robot,ListDAY(ind).name,'\',Rat]);
        for i_f = 3:length(ListFiles) % find the correct file
            if contains(ListFiles(i_f).name,'Peak_x04.mat')
                A = load ([data_robot,ListDAY(ind).name,'\',Rat,'\',ListFiles(i_f).name]);
                if strcmp(A.Data.info.Status,TaskR{i_tas})
                    filesave = [ListFiles(i_f).name(1:end-4),'_Ins.mat'];
                end
                clear A
            end
        end
        
        %Data to save
        Data = allData{i_gg,i_tas}.Data;
        Data.INSCOPIX.fS_ins = fSins;
        Data.INSCOPIX.cells_name = allCells{i_gg,i_tas,1};
        Data.INSCOPIX.cells_signal = allCells{i_gg,i_tas,2};
        Data.INSCOPIX.cells_pos = allCells{i_gg,i_tas,3};
        Data.INSCOPIX.cells_size = allCells{i_gg,i_tas,4};
        Data.INSCOPIX.time = allCells{i_gg,i_tas,5};
        Data.INSCOPIX.events = allCells{i_gg,i_tas,6};
        Data.INSCOPIX.cells_SNR = allCells{i_gg,i_tas,7};
        Data.INSCOPIX.Freq_rate = allCells{i_gg,i_tas,8};
        
        save([pathsave,filesave],'Data')
        clear filesave pathsave Data  
   end
end

%% figure all synch
% PLOT = 1;
% if PLOT
% fz = allData{i_gg,i_tas}.Data.Recorded_Data.Fz.data(1:fS_robot/fSins:end)*6000;
% fS_EMG = allData{i_gg,i_tas}.Data.VICON.fS_EMG;
% fS = allData{i_gg,i_tas}.Data.SIMI.fS_KIN;
% fS_ins = Data.INSCOPIX.fS_ins;
% cells =Data.INSCOPIX.cells_signal;
% figure
% set(gca,'FontSize',14)
% hold on
% subplot(4,1,1)
% xt = [1/fS_ins:1/fS_ins:length(cells)/fS_ins];
% %xt = [1/fS_robot:1/fS_robot:length(fz)/fS_robot];
% plot(xt(1:end-1),fz,'b')
% ylabel('Force (N)')
% % for j = 1:length(pks{1,3})
% %     in = pks{1,3};
% %     xx = [in(j)-0.1 in(j)+0.4 in(j)+0.4 in(j)-0.1];
% %     yy = [-5 -5 5 5];
% %     p = patch(xx,yy,[236,176,31]/255,'EdgeColor','none');
% %     set(p,'FaceAlpha',0.5)
% % end
% %xlim([0 80])
% subplot(4,1,2)
% hold on
% % for j = 1:length(pks{1,3})
% %     in = pks{1,3};
% %     xx = [in(j)-0.1 in(j)+0.4 in(j)+0.4 in(j)-0.1];
% %     yy = [-5 -5 5 5];
% %     p = patch(xx,yy,[236,176,31]/255,'EdgeColor','none');
% %     set(p,'FaceAlpha',0.5)
% % end
% %fS_EMG = Data.VICON.fS_EMG;
% Bi = Data.VICON.EMG(:,1);
% xt = [1/fS_EMG:1/fS_EMG:length(Bi)/fS_EMG];
% plot(xt,Bi,'m')
% xlim([0 80])
% ylabel('EMG Biceps (mV)')
% % subplot(4,1,3)
% % hold on
% % % for j = 1:length(pks{1,3})
% % %     in = pks{1,3};
% % %     xx = [in(j)-0.1 in(j)+0.4 in(j)+0.4 in(j)-0.1];
% % %     yy = [-4 -4 5 5];
% % %     p = patch(xx,yy,[236,176,31]/255,'EdgeColor','none');
% % %     set(p,'FaceAlpha',0.5)
% % % end
% % xt = [1/fS:1/fS:size(pos_tot)/fS];
% % plot(xt,(pos_tot(:,x(1))-mean(pos_tot(:,x(1))))/60,'g','LineWidth',1.5)
% % ylabel('Position')
% % xlim([0 80])
% subplot(4,1,4)
% %fS_ins = Data.INSCOPIX.fS_ins;
% %cells = Data.INSCOPIX.cells_signal;
% xt = [1/fS_ins:1/fS_ins:length(cells)/fS_ins];
% plot(xt,cells(:,:))
% %xlim([0 80])
% ylabel('Cells')
% xlabel('Time (s)')
% % for j = 1:length(pks{1,3})
% %     in = pks{1,3};
% %     xx = [in(j)-0.1 in(j)+0.4 in(j)+0.4 in(j)-0.1];
% %     yy = [0 0 5 5];
% %     p = patch(xx,yy,[236,176,31]/255,'EdgeColor','none');
% %     set(p,'FaceAlpha',0.5)
% % end
% % end
% 
