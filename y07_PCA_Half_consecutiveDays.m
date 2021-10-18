%% Plot parameters and PCA Push active(during pulling phase) along consecutive days (1 animal)
% To run this code it is necessary to have run until x04 for all animals of
% the days to analyse
% 1.  You can choose the days to visualize, if the mean of the parameters 
% is performed on the single animals or all together, if you want to
% normalize and also if you want to save results


clear all
%close all

DataToShow = {'PULL active'};
%name of the parameters that you want to plot and to use for the pca
NameParam = {'ForcePeaks','Width','Fmax/Fonset',...%'fwhm',
    'Fmax/fwhm','AUCfwhm','AUCmin','derFup','derFdown',...
    'SmoothF','Tpull','Latency','Npeaks',...
    'Frpeaks','MeanForce'... %'peak_amp_attem','fwhm_attem','width_attem','ratio_amp_attem',...'ratio_pk_fwhm_attem','AUC_fwhm_attem','AUC_min_attem','df_up_attem','df_down_attem','smoothness_attem',
    'Nburst','Ampburst','Maxburst','Tburst','Aucburst','MeanEnv'...
    'sub-movements','attempts','Vpeaks',...
    'height','amplitude','LenTra','AUC',...
    'SpeedTra','SmoothTra','MaxSpeedTra','Acc','MaxAcc'};
    
    %'height','amplitude','length_tr','length_x','length_y','AUC',... % These parameters of the kinematic are considered only during the pulling phase
    %'speed','speed_x','speed_y','smoothness_tr','smoothnees_x',...
    %'smoothness_y', 'speed_max','speed_max_x','speed_max_y'};%,...
    %'acceleration','acceleration_max'};
   
% {'peak_amp','fwhm','width','ratio_amp',...
%     'ratio_pk_fwhm','AUC_fwhm','AUC_min','df_up','df_down',...
%     'smoothness','t_retraction','t_firstPeak','n_peaks','std_amp',...
%     'peak_frequency','peak_amp_attem','fwhm_attem','width_attem','ratio_amp_attem',...
%     'ratio_pk_fwhm_attem','AUC_fwhm_attem','AUC_min_attem','df_up_attem','df_down_attem','smoothness_attem'...
%     'n_peaks_Fx','peak_amp_Fx','fwhm_Fx','width_Fx','ratio_amp_Fx',...
%     'smoothness_Fx','peak_frequency_Fx','n_burst_bi','burst_amp_bi',...
%     'burst_max_bi','burst_time_bi','burst_auc_bi','n_burst_tr',...
%     'burst_amp_tr','burst_max_tr','burst_time_tr','burst_auc_tr'...
%     'submov','attempts','speed_peaks','height','height_back','length'...
%     'mean_speed','smoothness_tr'};

NameForce = {'Fz','Fx'};
%select days to plot (IMPORTANT DPI in name of files recorded after injury, all BL days at the beginning!!!)
DAY = {'DAY02_2021_03_02_BL','DAY03_2021_03_03_BL','DAY04_2021_03_04_BL','DAY05_2021_03_05_BL','DAY08_2021_03_08_BL'...
    'DAY13_2021_03_13_05DPI','DAY14_2021_03_14_06DPI','DAY15_2021_03_15_07DPI'...
    'DAY16_2021_03_16_08DPI','DAY17_2021_03_17_09DPI','DAY18_2021_03_18_10DPI'}; 

%DAY = {'DAY02_2021_03_02_BL','DAY03_2021_03_03_BL','DAY04_2021_03_04_BL','DAY05_2021_03_05_BL','DAY08_2021_03_08_BL'...
    %'DAY01_03DPI','DAY02_04DPI','DAY03_07DPI'...
    %'DAY04_08DPI','DAY05_09DPI'};%,'DAY06_10DPI','DAY07_11DPI'}; 

nBL = 5;
dayPlot = {'BL','5DPI','6DPI','7DPI','8DPI','9DPI','10DPI'}; 
mycolor = {[209 211 212]./255,[253 229 205]./255,[241,181,135]./255,[234,110,47]./255,[240,76,2]./255,...
    [153,28,0]./255,[100,10,0]./255};
PCs = {'PC1','PC2','PC3'}; % select the PC that you want to use for the pca plot
NORM = 0; % 1 if you want normalize on BL otherwise 0
MEDrat = 0; % 0 only one animal!!!
SAVE = 0; % save all the figure for every parameter
PLOT_p = 0; % plot all parameters
REM = 1; % remove outliers PCA
PLOT3 = 1; % plot pca 3D
ch_PARAM = 1; %param to plot single peaks

datapath = 'C:\R-Platform\DATA\2021_02_Robot Stim group\Robot Data\';
Path = cd;

NumParam = length(NameParam);
LenDataTS = length(DataToShow);
%matrix to save all data (righe = Par) * (colonne = trattamento)
MatrixTot = cell(NumParam, LenDataTS);

vect_allpars = cell (size(DAY,2),size(DataToShow,2));
matrix_allpars = cell(size(DAY,2),size(DataToShow,2)); %one cell for every modality
info_rat = cell(size(DAY,2),size(DataToShow,2));
info_day = cell(size(DAY,2),size(DataToShow,2));

%% load data

for i_day = 1:size(DAY,2)
    ListRat = dir ([datapath,DAY{1,i_day}]); 
    if ~contains (DAY{1,i_day},'DPI') % we keep track of the day that we are analysing
        info_day{i_day,1} = 1; % all the baseline data can be analyzed together, so they are all zero
        nRAT = [4,6];
    else
        info_day{i_day,1} = str2num(DAY{1,i_day}(18:19)); % number of the day 
        nRAT = 7;
    end
%     if ~contains (DAY{1,i_day},'DPI') % we keep track of the day that we are analysing
%         info_day{i_day,1} = str2num(DAY{1,i_day}(end-4:end-3))-10;%1; % all the baseline data can be analyzed together, so they are all zero
%         datapath = 'C:\R-Platform\DATA\2021_02_Robot Stim group\Robot Data\';
%         nRAT = [4,6];
%     else
%         info_day{i_day,1} = str2num(DAY{1,i_day}(end-4:end-3)); % number of the day 
%         %datapath = 'C:\R-Platform\DATA\2021_02_Robot Stim group\Robot Array\';
%         nRAT = 7;
%     end
    ListRat = dir ([datapath,DAY{1,i_day}]);
    
    for i_rat = nRAT % Choose the animal (#4 -> R05)
        ListFile = dir([datapath,DAY{1,i_day},'\',ListRat(i_rat,1).name]);
        in_files = [];
        
        for i_file = 3:length(ListFile)
            if contains(ListFile(i_file,1).name,'x04.mat') && ~ contains(ListFile(i_file,1).name,'1DOF') && ~contains(ListFile(i_file,1).name,'_stim')&& ~contains(ListFile(i_file,1).name,'ell')
                in_files = [in_files,i_file];
            end             
        end
        
        for i_file_ok = 1:length(in_files)
            load ([datapath,DAY{1,i_day},'\',ListRat(i_rat,1).name,'\',ListFile(in_files(i_file_ok),1).name]);
            if contains(Data.info.Status,'Push_active')
            
                parametri = Data.Recorded_Data.Analysis.Fzpeaks;
                deltaFcy = Data.Recorded_Data.Analysis.Area_cy;
                parametriFx = Data.Recorded_Data.Analysis.Fxpeaks;
                burstBI = Data.VICON.burst{1,1};
                burstTR = Data.VICON.burst{2,1};
                deltaEMGcy = Data.VICON.Analysis.Area_cy;
                if isfield (Data.SIMI,'Analysis')
                    KIN = Data.SIMI.Analysis;
                else
                    KIN = NaN;
                end
                                
                %Data Fz
                for i = 1: length(Data.good_trials(:,1)) % we want to consider only the trials that have been selected as good trials
                    index_pks = find(parametri(:,14)==Data.good_trials(i,1) & parametri(:,15)==2 & parametri(:,16)==1); %we want to analize only data coming from the pulling phase
                    index_att = find(parametri(:,14)==Data.good_trials(i,1) & parametri(:,15)==2 & parametri(:,16)==0);
                    if ~isempty(index_pks)
                        samples = find(Data.Recorded_Data.cicles.data==Data.good_trials(i,1) & Data.Recorded_Data.T_status.data==2);
                        t_tot = length(samples)/Data.Recorded_Data.fS_robot; % total time to complete a trial
                        t_first = (parametri(index_pks(1),1)-samples(1))/Data.Recorded_Data.fS_robot; % time before the first peak
                        n_peaks = length(index_pks)+length(index_att); %total number of peaks for each trial
                        if length(index_pks)>3
                            std_amp = std(parametri(index_pks,4)); % std of the amplitude of different peaks in one cicle
                        else
                            std_amp = 0;
                        end
                        peak_freq = n_peaks/t_tot;
                        %line1 = [mean(parametri(index_pks,4:13),1),t_tot,t_first,n_peaks,std_amp,peak_freq,mean(parametri(index_att,4:13),1),deltaFcy(i)];  
                        line1 = [nanmean(parametri(index_pks,4:13),1),t_tot,t_first,n_peaks,peak_freq,deltaFcy(i)];
                        %line1 = [min(parametri(index_pks,4)),nanmean(parametri(index_pks,5:13),1),t_tot,t_first,n_peaks,peak_freq,deltaFcy(i)];
                    else                        
                        line1 = NaN(1,15);
                    end
                    line1(2) = [];
%                     % Data FX
%                     index = find (Data.Recorded_Data.cicles.data(parametriFx(:,1))==Data.good_trials(i,1)& Data.Recorded_Data.T_status.data(parametriFx(:,1))==2); %we want to analize only data coming from the pulling phase
%                     if ~isempty(index)
%                         samples = find(Data.Recorded_Data.cicles.data==Data.good_trials(i,1) & Data.Recorded_Data.T_status.data==2);
%                         t_tot = length(samples)/Data.Recorded_Data.fS_robot; % total time to complete a trial
%                         n_peaks_Fx = length(index); %total number of peaks for each trial
%                         peak_freq_Fx = n_peaks/t_tot;
%                         line2 = [n_peaks_Fx,mean(parametriFx(index,4:7),1),mean(parametriFx(index,13),1),peak_freq_Fx];                      
%                     else                        
%                         line2 = NaN(1,7);
%                     end
%                     
%                     % Data EMG BI
                    BIst = round(burstBI(:,1)*Data.Recorded_Data.fS_robot);
                    BIst(BIst==0)=[];% delete of the bursts that starts before the beginnning
                    index = find (Data.Recorded_Data.cicles.data(BIst)==Data.good_trials(i,1));%& Data.Recorded_Data.T_status.data(BIst)==2); %we want to analize only data coming from the pulling phase
                    if ~isempty(index)
%                         if index(1)~=1
%                             tt = round(BIst(index(1)-1)+burstBI(index(1)-1,4)*Data.Recorded_Data.fS_robot); % we want to conside also the burst that start before phase 2 and then it finishes during phase 2, because they can be very long
%                             if Data.Recorded_Data.cicles.data(tt)==Data.good_trials(i,1)&& Data.Recorded_Data.T_status.data(tt)==2
%                                 index = [index(1)-1,index];
%                             end
%                         end
                        n_burst = length(index); %total number of burst for each trial
                        line3 = [n_burst,nanmean(burstBI(index,2:5),1),deltaEMGcy(i,1)];                      
                    else                        
                        line3 = NaN(1,6);
                    end
                    
                    % Data EMG TR
%                     TRst = round(burstTR(:,1)*Data.Recorded_Data.fS_robot);
%                     TRst(TRst==0)=[];
%                     index = find (Data.Recorded_Data.cicles.data(TRst)==Data.good_trials(i,1)& Data.Recorded_Data.T_status.data(TRst)==2); %we want to analize only data coming from the pulling phase
%                     if ~isempty(index) && i_rat~=4  %number 2 has the broken triceps
%                         if index(1)~=1
%                             tt = round(TRst(index(1)-1)+burstTR(index(1)-1,4)*Data.Recorded_Data.fS_robot); % we want to conside also the burst that start before phase 2 and then it finishes during phase 2, because they can be very long
%                             if Data.Recorded_Data.cicles.data(tt)==Data.good_trials(i,1)&& Data.Recorded_Data.T_status.data(tt)==2
%                                 index = [index(1)-1,index];
%                             end
%                         end
%                         n_burst = length(index); %total number of burst for each trial
%                         line4 = [n_burst,mean(burstTR(index,2:5),1)];                      
%                     else                        
%                         line4 = NaN(1,5);
%                     end
                    
                    % Data KINEMATIC
                    line5 = [];
                    if isstruct(KIN)
                        list_par = fieldnames(KIN);
                        for i_kin = 1:length(list_par) 
                            val = KIN.(list_par{i_kin})(:,end);
                            %val(isoutlier(val))=NaN;
                            line5 = [line5 , val(i)]; % val good trial
                        end
                    else
                        for i_kin = 1:18
                            line5 = [line5 , NaN];
                        end
                    end
                    line5(17:18) = [];
                    line5(14:15)=[];
                    line5(11:12)=[];
                    line5(7:8)=[];
            
                    %Join all data and insert in the matrix, create a
                    %matrix with the animal info (to normalize) and
                    %the day info (to study the evolution over time)
                    
                    line_tot = [line1,line3,line5];
                    if sum(isnan(line_tot))<length(line_tot)/3
                           matrix_allpars{i_day,1}= [matrix_allpars{i_day,1};line1,line3,line5];
                           info_rat{i_day,1} = [info_rat{i_day};str2num(ListRat(i_rat,1).name(2:3))];
                    end
                    
                end
                
            end    
        end
        
    end    
    
end

%% JOIN all data (all animals, all days) in the same matrix to do the PCA
%(in this moment exist only one raw so it is not useful the other cycle
%for)
info_data = [];
matrix_allparams_allanimals = [];
for i_raw = 1: size(matrix_allpars,1)
    for i_col = 1: size(matrix_allpars,2)
        matrix_allparams_allanimals = [matrix_allparams_allanimals; matrix_allpars{i_raw,i_col}];
        info_data = [info_data; ones(size(matrix_allpars{i_raw,i_col},1),1)*info_day{i_raw,i_col}, info_rat{i_raw,i_col}]; %first column for the number of day, second one for the number of rat
    end
end

%% Normalization and mean inside the same animal

if NORM
    %create a matrix: an animal for every raw, a parameter for every column
    matr_to_norm = NaN(size(ListRat,1)-3,size(NameParam,2));
    for j_rat = 1:size(matr_to_norm,1)
        n_rat = str2num(ListRat(j_rat+2,1).name(2:3)); %number of the animal, saved in the info
        in_norm = find (info_data(:,1)== 0 & info_data(:,2) == n_rat); % 0 -> days Baseline
        matr_to_norm(j_rat,:) = nanmean(matrix_allparams_allanimals(in_norm,:),1);
    end
        matrix_allparams_allanimals_norm = matrix_allparams_allanimals;
       for j_rat = 1:size(ListRat,1)-3
           n_rat =str2num(ListRat(j_rat+2,1).name(2:3));
           in_med = find (info_data(:,2) == n_rat);
           matrix_allparams_allanimals_norm(in_med,:)= matrix_allparams_allanimals(in_med,:)./matr_to_norm(j_rat,:);
       end
       tot_days = unique(info_data(:,1));
       max_cycles = NaN(length(tot_days),1);
       for j_day = 1 : length(tot_days)
           max_cycles(j_day) = length(find(info_data(:,1) == tot_days(j_day)));
       end
       MM = NaN(max(max_cycles),size(NameParam,2),size(DAY,2)-nBL+1);
       for j_day = 1 : length(tot_days)
           in_med = find(info_data(:,1) == tot_days(j_day));
           MM(1:max_cycles(j_day),:,j_day) = matrix_allparams_allanimals_norm(in_med,:);
       end
else
       tot_days = unique(info_data(:,1));
       max_cycles = NaN(length(tot_days),1);
       for j_day = 1 : length(tot_days)
           max_cycles(j_day) = length(find(info_data(:,1) == tot_days(j_day)));
       end
       MM = NaN(max(max_cycles),size(NameParam,2),size(DAY,2)-nBL+1);
       for j_day = 1 : length(tot_days)
           in_med = find(info_data(:,1) == tot_days(j_day));
           MM(1:max_cycles(j_day),:,j_day) = matrix_allparams_allanimals(in_med,:);
       end
%     end
end

%% Figure and stat for every parameter

nt_p = NumParam; % total number of parameter
nt_d = size(DAY,2)-nBL+1; %total number of days (baseline all together)
xbar = [1:1:nt_d];
  
errhigh = NaN (nt_d,nt_p);
dataPlot = NaN (nt_d,nt_p);
for k_day = 1:nt_d
%     if k_day ==1
%         MMBL =[];
%         for i_bl = 1:nBL
%             MMBL = [MMBL;MM(:,:,i_bl)];
%         end
%         dataPlot(1,:) = nanmean(MMBL,1);
%         errhigh(1,:) = nanstd(MMBL,1)/sqrt(size(MMBL,1));
%     else
        dd = MM(:,:,k_day); % because there are more baseline days!
        dd(isnan(dd(:,1)),:)= [];
        errhigh(k_day,:) = nanstd(dd,1)/sqrt(size(dd,1));
        dataPlot(k_day,:) = nanmean(dd,1);
    %end
end

if PLOT_p
    
for j_par = 1: nt_p
    
    figure(j_par+700)    
    hold on
    bb=bar(xbar,dataPlot(:,j_par),'stacked');
    bb.EdgeColor = 'flat';
    bb.FaceColor ='flat';
%     bb.CData(1,:) = mycolor{1};
%     bb.CData(2,:)= mycolor{3};
%     bb.CData(3,:)= mycolor{5};
    ylabel (NameParam{j_par},'FontSize',13,'FontName','Arial')
    hold on
    er = errorbar(xbar,dataPlot(:,j_par),-errhigh(:,j_par),errhigh(:,j_par),'HandleVisibility','off');
    er.Color = [0 0 0];                            
    er.LineStyle = 'none';
    er.LineWidth = 1;
    xlim([0 length(xbar)+1])
    xticks(xbar)
        
    x_stat = []; % in x_stat data for statistic are collected
    group = []; % in group data for statistic are collected
    for jj = 1 :nt_d
%         if jj ==1
%             for i_bl = 1:nBL
%                 x_stat = [x_stat;MM(:,j_par,i_bl)];
%                 group = [group;ones(length(MM(:,j_par,i_bl)),1)];
%             end
%         else
        x_stat = [x_stat;MM(:,j_par,jj)];
        group = [group;ones(length(MM(:,j_par,jj)),1)*jj];
        %end
    end
    xticklabels(dayPlot);
    
    if SAVE
        cd([datapath,'Results'])
        if NORM
            name2save =[NameParam{j_par},'_norm'];
        else
            name2save =[NameParam{j_par}];
        end
        saveas(gcf,name2save)
    end
    if ~NORM
    %statistic
    [p,tbl,stats] = kruskalwallis(x_stat,group);
    c = multcompare(stats);
    title(NameParam{j_par});
    if SAVE
        cd([datapath,'Results'])
        saveas(gcf,[NameParam{j_par},'stat'])
    end
    end
    
end
end

%% Calculate PCA
zscor_xnan = @(x) bsxfun(@rdivide, bsxfun(@minus, x, mean(x,'omitnan')), std(x, 'omitnan'));
matrix_allparams_allanimals_z =  zscor_xnan(matrix_allparams_allanimals);

cd([Path,'\SubFunctions'])
%matrix_allparams_allanimals_z = Reconstruct_matrix_PCA (matrix_allparams_allanimals_z,info_data,REM);

[COEFF, SCORE, LATENT,~,explained] = pca((matrix_allparams_allanimals_z));%,'algorithm','eig','Rows','pairwise');
figure; imagesc(COEFF)

%nanmean(matrix_allparams_allanimals_z,1);
%nanstd(matrix_allparams_allanimals_z,1);
 
 figure(1103)
 hold on
 set(gca,'FontSize',13)
 legend
 xlabel ([PCs{1},'-',num2str(explained(str2num(PCs{1}(3)))),'%'],'FontSize',18,'FontWeight','bold')
 ylabel ([PCs{2},'-',num2str(explained(str2num(PCs{2}(3)))),'%'],'FontSize',18,'FontWeight','bold')

 cl=1;
 nn = 1;
 nt_g = unique(info_data(:,1));
 yp_all = cell(1,length(nt_g));
 xp_all = cell(1,length(nt_g));
 group_stat = [];
 xp_stat = [];
 yp_stat = [];
for pl_day =1: length(nt_g)
     ind = find(info_data(:,1) ==nt_g(pl_day));
     scatter(nanmean(SCORE(ind,str2num(PCs{1}(3)))),nanmean(SCORE(ind,str2num(PCs{2}(3)))),200,'MarkerFaceColor', mycolor{cl},'MarkerEdgeColor','k','DisplayName',dayPlot{nn})
     scatter(SCORE(ind,str2num(PCs{1}(3))),SCORE(ind,str2num(PCs{2}(3))),20,'MarkerFaceColor', mycolor{cl},'MarkerEdgeColor', mycolor{cl}, 'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.5,'HandleVisibility','off')
     cl=cl+1; 
     yp_all {nn} = SCORE(ind,str2num(PCs{2}(3)));
     xp_all {nn} = SCORE(ind,str2num(PCs{1}(3)));
     group_stat = [group_stat; ones(length(yp_all{nn}),1)*pl_day];
     xp_stat = [xp_stat; xp_all{nn}];
     yp_stat = [yp_stat; yp_all{nn}];
     x_point{nn} = ones(length(yp_all{nn}),1)*nn;
     nn = nn+1;
 end
 
[px,~,stat] = kruskalwallis(xp_stat,group_stat);
[py,~,stat] = kruskalwallis(yp_stat,group_stat);

%% Bar plot distances and stat 

errscore = NaN (length(nt_g),2);
scorePlot = NaN (length(nt_g),2);
n = 0;

for k_gg =1:length(nt_g)
    dd = SCORE(info_data(:,1)==nt_g(k_gg),1:2);
    if ~isempty(dd)
        dd(isnan(dd(:,1)),:)= [];
        errscore(k_gg,:) = nanstd(dd,1)/sqrt(size(dd,1));
        scorePlot(k_gg,:) = nanmean(dd,1);
    end
end
    
xbar = [1:nt_d];
for i_pc = 1:2
    figure
    hold on
    bb=bar(xbar,scorePlot(:,i_pc),'stacked');
    bb.EdgeColor = 'flat';
    bb.FaceColor ='flat';
    bb.CData(1,:) = mycolor{1};
    bb.CData(2,:)= mycolor{2};
    bb.CData(3,:)= mycolor{3};
    bb.CData(4,:)= mycolor{4};
    bb.CData(5,:)= mycolor{5};
    bb.CData(6,:)= mycolor{6};
    bb.CData(7,:)= mycolor{7};
    ylabel (PCs{i_pc} ,'FontSize',13,'FontName','Arial')
    hold on
    er = errorbar(xbar,scorePlot(:,i_pc),-errscore(:,i_pc),errscore(:,i_pc),'HandleVisibility','off');
    er.Color = [0 0 0];                            
    er.LineStyle = 'none';
    er.LineWidth = 1;
    scorestat = [];
    scoregroup = [];
    xlim([0 length(xbar)+1])
    for i_dd = 1:length(nt_g)
        x_point = ones(length(SCORE(info_data(:,1)==nt_g(i_dd),i_pc)),1)*i_dd;
        scatter(x_point,SCORE(info_data(:,1)==nt_g(i_dd),i_pc),20,'k','MarkerFaceColor','k','MarkerEdgeColor','none')
        scorestat = [scorestat;SCORE(info_data(:,1)==nt_g(i_dd),i_pc)];
        scoregroup = [scoregroup;ones(length(SCORE(info_data(:,1)==nt_g(i_dd),i_pc)),1)*i_dd];  
    end
    [p,~,stat] = kruskalwallis(scorestat,scoregroup);
    c = multcompare(stat);    
end


%% Fig more important parameters

pc = [COEFF(:,str2num(PCs{1}(3))),COEFF(:,str2num(PCs{2}(3))),COEFF(:,str2num(PCs{3}(3)))];
pc = abs(pc);
xbar = [1:nt_d];
x_point = ones(size(MM,1),1);
for i_pc = 1:size(pc,2)
    figure(1080+i_pc)
     hold on
    [~,in] = max(pc(:,i_pc));
    subplot (1,2,1)
    pc(in,i_pc) = 0;
    bb=bar(xbar,dataPlot(:,in),'stacked');
    hold on
    bb.EdgeColor = 'flat';
    bb.FaceColor ='flat';
    bb.CData(1,:) = mycolor{1};
    bb.CData(2,:)= mycolor{2};
    bb.CData(3,:)= mycolor{3};
    bb.CData(4,:)= mycolor{4};
    bb.CData(5,:)= mycolor{5};
    bb.CData(6,:)= mycolor{6};
    bb.CData(7,:)= mycolor{7};
    
    ylabel (NameParam{in},'FontSize',13,'FontName','Arial')
    hold on
    er = errorbar(xbar,dataPlot(:,in),-errhigh(:,in),errhigh(:,in),'HandleVisibility','off');
    er.Color = [0 0 0];                            
    er.LineStyle = 'none';
    er.LineWidth = 1;
    for i_dd = 1:nt_d
        scatter(x_point*i_dd,MM(:,in,i_dd),20,'k','MarkerFaceColor','k','MarkerEdgeColor','none')
    end
    xlim([0 length(xbar)+1])
    xticks(xbar)
    xticklabels(dayPlot);
    
    [~,in] = max(pc(:,i_pc));
    subplot (1,2,2)
    hold on
    pc(in,i_pc) = 0;
    bb=bar(xbar,dataPlot(:,in),'stacked');
    bb.EdgeColor = 'flat';
    bb.FaceColor ='flat';
    bb.CData(1,:) = mycolor{1};
    bb.CData(2,:)= mycolor{2};
    bb.CData(3,:)= mycolor{3};
    bb.CData(4,:)= mycolor{4};
    bb.CData(5,:)= mycolor{5};
    bb.CData(6,:)= mycolor{6};
    bb.CData(7,:)= mycolor{7};
    
    ylabel (NameParam{in},'FontSize',13,'FontName','Arial')
    hold on
    er = errorbar(xbar,dataPlot(:,in),-errhigh(:,in),errhigh(:,in),'HandleVisibility','off');
    er.Color = [0 0 0];                            
    er.LineStyle = 'none';
    er.LineWidth = 1;
    for i_dd = 1:nt_d
        scatter(x_point*i_dd,MM(:,in,i_dd),20,'k','MarkerFaceColor','k','MarkerEdgeColor','none')
    end
    xlim([0 length(xbar)+1])
    xticks(xbar)
    xticklabels(dayPlot);

end


%% Create color map for all parameters for PCs

figure;
imagesc(abs(COEFF(:,1:3)))
colormap('Copper')

yticks(1:32)
yticklabels(NameParam)
xticks(1:3)
xticklabels(PCs)

%% Scatter single parameter
%for ch_PARAM =1:length(NameParam)
figure;
hold on
interpolation = [];
ind = find(info_data(:,1)==nt_g(1)); % BL indexes
bl_line = nanmean(matrix_allparams_allanimals(ind,ch_PARAM));
st_line = nanstd(matrix_allparams_allanimals(ind,ch_PARAM))/sqrt(length(ind));
xline = [1,length(nt_g)];
yline = [bl_line,bl_line];
X=[xline, fliplr(xline)];
Y=[yline + st_line,fliplr(yline - st_line)];
fill( X,Y,mycolor{1});
alpha(.10)
plot(xline,yline,'Color',mycolor{1})
plot(xline,yline-st_line,'k--')
plot(xline,yline+st_line,'k--')
for i = 2:length(nt_g)
    dd = matrix_allparams_allanimals(info_data(:,1)==nt_g(i),ch_PARAM);
    dd(isoutlier(dd,'quartiles'))=[];
    xpoint = linspace(i-1,i,length(dd));
    scatter(xpoint,dd,40,mycolor{i},'filled')
    if i~=1 %i>5 %~=1
        interpolation = [interpolation;xpoint',dd];
    end
end

p = polyfit(interpolation(:,1),interpolation(:,2),1);
yline = xline*p(1)+p(2);
plot(xline,yline,'k','LineWidth',1.5)

ylabel(NameParam{ch_PARAM})
xlabel('Time (DPI)')
%end

