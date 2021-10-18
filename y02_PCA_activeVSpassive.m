%% PCA Half (during pulling phase) vs Active 
% To run this code it is necessary to have run until x03 for all animals of
% the days to analyse
% It works for healthy or for post injury (just need to select the correct
% days)
% 1.  You can choose the PCs to visualize, the days to visualize and if the
% the mean of the parameters is performed on the single animals or all
% together

clear all
close all

DataToShow = {'PULL task','PASSIVE task'};
%NameParam = {'peak_amp','fwhm','width','ratio_amp',...
 %   'ratio_pk_fwhm','AUC_fwhm','AUC_min','df_up','df_down',...
 %   'smoothness','t_retraction','t_firstPeak','n_peaks'};
 
% NameParam = {'peak_amp','width','ratio_amp',...%'fwhm',
%     'ratio_pk_fwhm','AUC_fwhm','AUC_min','df_up','df_down',...
%     'smoothness','t_retraction','t_firstPeak','n_peaks',...
%     'peak_frequency','deltaFcy','n_burst_bi','burst_amp_bi',...
%     'burst_max_bi','burst_time_bi','burst_auc_bi','deltaBIcy',...    
%     'height','amplitude','length_tr','AUC',...
%     'speed','smoothness_tr','speed_max','acceleration','acceleration_max'};

NameParam = {'ForcePeaks','Width','Fmax/Fonset',...%'fwhm',
    'Fmax/fwhm','AUCfwhm','AUCmin','derFup','derFdown',...
    'SmoothF','Tpull','Latency','Npeaks',...
    'Frpeaks','MeanForce'... %'peak_amp_attem','fwhm_attem','width_attem','ratio_amp_attem',...'ratio_pk_fwhm_attem','AUC_fwhm_attem','AUC_min_attem','df_up_attem','df_down_attem','smoothness_attem',
    'Nburst','Ampburst','Maxburst','Tburst','Aucburst','MeanEnv'...
    'height','amplitude','LenTra','AUC',...
    'SpeedTra','SmoothTra','MaxSpeedTra','Acc','MaxAcc'};
    
    %'height','amplitude','length_tr','length_x','length_y','AUC',... % I choose the parameters for all trajectory, we can also choose them only for one phase
    %'speed','speed_x','speed_y','smoothness_tr','smoothnees_x',...
    %'smoothness_y', 'speed_max','speed_max_x','speed_max_y',...
    %'acceleration','acceleration_max'};

    %,'n_peaks_Fx','peak_amp_Fx','fwhm_Fx','width_Fx','ratio_amp_Fx','smoothness_Fx','peak_frequency_Fx'
    %'n_burst_bi','burst_amp_bi',...
    %'burst_max_bi','burst_time_bi','burst_auc_bi','n_burst_tr',...
    %'burst_amp_tr','burst_max_tr','burst_time_tr','burst_auc_tr'};

NameForce = {'Fz','Fx'};
REM = 0; % remove outliers PCA
MEDrat = 1 ; % if parameters are plot as a mean of every rat
PCs = {'PC1','PC2','PC3'}; % select the PC that you want to use for the pca plot
PLOT3 = 1; % to plot PCA in 3D space
NumParam = length(NameParam);
LenDataTS = length(DataToShow);
%matrice dove salverò tutto (righe = Par) * (colonne = trattamento)
MatrixTot = cell(NumParam, LenDataTS);
DAY = {'DAY02_2019_08_19','DAY03_2019_08_20','DAY04_2019_08_23','DAY05_2019_08_28','DAY06_2019_09_11','DAY07_2019_09_13','DAY08_2019_09_18'};
mycolor = {[250 82 91]./255,[25 211 197]./255};
%mycolor = {[243 146 0]./255,[128 0 128]./255,[100 0 255]./255,[255 237 0]./255,...
%     [255 0 0]./255,[0 255 0]./255,[230 0 230]./255};

vect_allpars = cell (size(DAY,2),size(DataToShow,2));
matrix_allpars = cell(size(DAY,2),size(DataToShow,2)); %one cell for every modality
info_rat = cell(size(DAY,2),size(DataToShow,2));
info_day = cell(size(DAY,2),size(DataToShow,2));


%% load data
datapath = 'C:\R-Platform\DATA\2019_02_Group1\Recordings_robot\';
Path = cd;

for i_day = 1:size(DAY,2)
    ListRat = dir ([datapath,DAY{1,i_day}]);
    %if ~contains (DAY{1,i_day},'DPI') % we keep track of the day that we are analysing
    %    info_day{i_day,1} = 0; % all the baseline data can be analyzed together, so they are all zero
    %    info_day{i_day,2} = 0;
    %else
        info_day{i_day,1} = str2num(DAY{1,i_day}(4:5)); % number of the day 
        info_day{i_day,2} = str2num(DAY{1,i_day}(4:5)); % number of the day 
    %end
    
    for i_rat = 3:10 
        ListFile = dir([datapath,DAY{1,i_day},'\',ListRat(i_rat,1).name]);
        in_files = [];
        for i_file = 3:length(ListFile)
            if contains(ListFile(i_file,1).name,'x04.mat') 
                in_files = [in_files,i_file];
            end             
        end
        
        for i_file_ok = 1:length(in_files)
            load ([datapath,DAY{1,i_day},'\',ListRat(i_rat,1).name,'\',ListFile(in_files(i_file_ok),1).name]);
            
            
                parametri = Data.Recorded_Data.Analysis.Fzpeaks;
                parametriFx = Data.Recorded_Data.Analysis.Fxpeaks;
                deltaFcy = Data.Recorded_Data.Analysis.Area_cy;
                deltaEMGcy = Data.VICON.Analysis.Area_cy;
                burstBI = Data.VICON.burst{1,1};
                burstTR = Data.VICON.burst{2,1};
                if isfield (Data.SIMI,'Analysis')
                    KIN = Data.SIMI.Analysis;
                else
                    KIN = NaN;
                end
                
                %Data Fz
                for i = 1: length(Data.good_trials(:,1)) % we want to consider only the trials that have been selected as good trials
                    if contains(Data.info.Status,'Push_active')%we want to analize only real peaks of the animals
                        %index = find (Data.Recorded_Data.cicles.data(parametri(:,1))==Data.good_trials(i,1)& Data.Recorded_Data.T_status.data(parametri(:,1))==2); 
                        index = find(parametri(:,14)==Data.good_trials(i,1) & parametri(:,15)==2);
                    else
                        %index = find (Data.Recorded_Data.cicles.data(parametri(:,1))==Data.good_trials(i,1));
                        index = find(parametri(:,14)==Data.good_trials(i,1));
                    end
                    if ~isempty(index)
                        samples = find(Data.Recorded_Data.cicles.data==Data.good_trials(i,1) & Data.Recorded_Data.T_status.data==2);
                        t_tot = length(samples)/Data.Recorded_Data.fS_robot; % total time to complete a trial
                        t_first = (parametri(index(1),1)-samples(1))/Data.Recorded_Data.fS_robot; % time before the first peak
                        n_peaks = length(index); %total number of peaks for each trial
                        if length(index)>3
                            std_amp = std(parametri(index,4)); % std of the amplitude of different peaks in one cicle
                        else
                            std_amp = NaN;
                        end
                        peak_freq = n_peaks/t_tot;
                        line1 = [mean(parametri(index,4:13),1),t_tot,t_first,n_peaks,peak_freq,deltaFcy(i)];%std_amp,peak_freq];                      
                    else                        
                        line1 = NaN(1,15);
                        line1(11) = length(samples)/Data.Recorded_Data.fS_robot;
                        line1(13) = 0;
                        line1(14) = 0;
                        %line1(15) = 0;
                    end
                    line1(2) = [];
%                     % Data FX
%                     if contains(Data.info.Status,'Push_active')%we want to analize only real peaks of the animals
%                         index = find (Data.Recorded_Data.cicles.data(parametriFx(:,1))==Data.good_trials(i,1)& Data.Recorded_Data.T_status.data(parametriFx(:,1))==2); 
%                     else
%                         index = find (Data.Recorded_Data.cicles.data(parametriFx(:,1))==Data.good_trials(i,1)); 
%                     end
%                     if ~isempty(index)
%                         samples = find(Data.Recorded_Data.cicles.data==Data.good_trials(i,1) & Data.Recorded_Data.T_status.data==2);
%                         t_tot = length(samples)/Data.Recorded_Data.fS_robot; % total time to complete a trial
%                         n_peaks_Fx = length(index); %total number of peaks for each trial
%                         peak_freq_Fx = n_peaks/t_tot;
%                         line2 = [n_peaks_Fx,mean(parametriFx(index,4:7),1),mean(parametriFx(index,13),1),peak_freq_Fx];                      
%                     else                        
%                         line2 = NaN(1,7);
%                         line1(7) = 0;
%                     end
                    
                    % Data EMG BI
                    BIst = round(burstBI(:,1)*Data.Recorded_Data.fS_robot);
                    BIst(BIst==0)=[];% delete of the bursts that starts before the beginnning
                    if contains(Data.info.Status,'Push_active')%we want to analize only real peaks of the animals
                        index = find (Data.Recorded_Data.cicles.data(BIst)==Data.good_trials(i,1)& Data.Recorded_Data.T_status.data(BIst)==2); 
                    else
                        index = find (Data.Recorded_Data.cicles.data(BIst)==Data.good_trials(i,1)); 
                    end
                    if ~isempty(index)
                        if index(1)~=1
                            tt = round(BIst(index(1)-1)+burstBI(index(1)-1,4)*Data.Recorded_Data.fS_robot); % we want to conside also the burst that start before phase 2 and then it finishes during phase 2, because they can be very long
                            if Data.Recorded_Data.cicles.data(tt)==Data.good_trials(i,1)&& Data.Recorded_Data.T_status.data(tt)==2
                                index = [index(1)-1,index];
                            end
                        end
                        n_burst = length(index); %total number of burst for each trial
                        line3 = [n_burst,mean(burstBI(index,2:5),1),deltaEMGcy(i,1)];                      
                    else                        
                        line3 = NaN(1,6);
                        line3(1) = 0;
                    end
                    
                    % Data EMG TR
                    TRst = round(burstTR(:,1)*Data.Recorded_Data.fS_robot);
                    TRst(TRst==0)=[];
                    if contains(Data.info.Status,'Push_active')%we want to analize only real peaks of the animals
                        index = find (Data.Recorded_Data.cicles.data(TRst)==Data.good_trials(i,1)& Data.Recorded_Data.T_status.data(TRst)==2); 
                    else
                        index = find (Data.Recorded_Data.cicles.data(TRst)==Data.good_trials(i,1)); 
                    end
                    if ~isempty(index) && i_rat~=4  %number 2 has the broken triceps
                        if index(1)~=1
                            tt = round(TRst(index(1)-1)+burstTR(index(1)-1,4)*Data.Recorded_Data.fS_robot); % we want to consider also the burst that start before phase 2 and then it finishes during phase 2, because they can be very long
                            if Data.Recorded_Data.cicles.data(tt)==Data.good_trials(i,1)&& Data.Recorded_Data.T_status.data(tt)==2
                                index = [index(1)-1,index];
                            end
                        end
                        n_burst = length(index); %total number of burst for each trial
                        line4 = [n_burst,mean(burstTR(index,2:5),1)];                      
                    else                        
                        line4 = NaN(1,5);
                        line4(1) = 0;
                    end
                    
                    % Data KINEMATIC
                    line5 = [];
                    if isstruct(KIN)
                    list_par = fieldnames(KIN);
                    for i_kin = 4:length(list_par) % first three parameters are possible only for half task
                        val = KIN.(list_par{i_kin})(:,1);
                        val(isoutlier(val))=NaN;
                        line5 = [line5 , val(i)]; % val good trial
                    end
                    else
                        for i_kin = 1:17
                            line5 = [line5 , NaN];
                        end
                    end
                    line5(14:15) = [];
                    line5(11:12)=[];
                    line5(8:9)=[];
                    line5(4:5)=[];
                    
                    %Join all data and insert in the matrix, create a
                    %matrix with the animal info (to normalize) and
                    %the day info (to study the evolution over time)
                    line_tot = [line1,line3,line5];
                    if sum(isnan(line_tot))<1
                                      
                    if  strcmp(Data.info.Status,'Push_active')
                       matrix_allpars{i_day,1}= [matrix_allpars{i_day,1};line1,line3,line5];%,line2,line3,line4];
                       info_rat{i_day,1} = [info_rat{i_day,1};str2num(ListRat(i_rat,1).name(2:3))];
                    else
                       matrix_allpars{i_day,2}= [matrix_allpars{i_day,2};line1,line3,line5];%,line2,line3,line4];
                       info_rat{i_day,2} = [info_rat{i_day,2};str2num(ListRat(i_rat,1).name(2:3))];
                    end
                    
                    end
                    
                end

        end        
        
    end    
    
end

%% JOIN all data (all animals, all days) in the same matrix to do the PCA

%here we want to study the comparison between active and passive
info_data = [];
matrix_allparams_allanimals = [];
for i_raw = 1: size(matrix_allpars,1)
    for i_col = 1: size(matrix_allpars,2)
        matrix_allparams_allanimals = [matrix_allparams_allanimals; matrix_allpars{i_raw,i_col}];
        info_data = [info_data; ones(size(matrix_allpars{i_raw,i_col},1),1)*info_day{i_raw,i_col}, info_rat{i_raw,i_col},ones(size(matrix_allpars{i_raw,i_col},1),1)*i_col]; %first column for the day, second one for the number of rat, third one for the task
    end
end

%% Fill missing data
% It is necessary to find a way to fill the NaN value because otherwise we
% miss important data, especially for the active task, where it is possible
% to don't have any peaks during the trial
% for i_col = 1 : size(matrix_allparams_allanimals,2) % in the columns there are the different parameters
%     xnan = find(isnan(matrix_allparams_allanimals(:,i_col))==1);   
%     for ix = 1:length(xnan)
%         %find date coming from the same day, same animal, same task
%         index_dateMean = find(info_data(:,1)==info_data(xnan(ix),1) & info_data(:,2)==info_data(xnan(ix),2) & info_data(:,3)==info_data(xnan(ix),3));
% %         if i_col ==1
% %             value = min(matrix_allparams_allanimals(index_dateMean,i_col));
% %         else
% %             value = min(matrix_allparams_allanimals(index_dateMean,i_col));
% %         end
%         pos_randi = randi(length(index_dateMean));
%         value = nanmean(matrix_allparams_allanimals(index_dateMean,i_col));
%         matrix_allparams_allanimals(xnan,i_col) = value;
%     end
% end

%% Calculate PCA
zscor_xnan = @(x) bsxfun(@rdivide, bsxfun(@minus, x, mean(x,'omitnan')), std(x, 'omitnan'));
matrix_allparams_allanimals_z =  zscor_xnan(matrix_allparams_allanimals);

%[COEFF, SCORE, LATENT,~,explained] = pca((matrix_allparams_allanimals_z),'algorithm','als');
%cd([Path,'\SubFunctions'])
%matrix_allparams_allanimals_z = Reconstruct_matrix_PCA (matrix_allparams_allanimals_z,info_data(:,3),REM);

[COEFF, SCORE, LATENT,~,explained] = pca((matrix_allparams_allanimals_z));%,'algorithm','eig','Rows','pairwise');

figure; imagesc(COEFF)
 
 figure(1002)
 hold on
 set(gca,'FontSize',13)
 legend
 xlabel ([PCs{1},'-',num2str(explained(str2num(PCs{1}(3)))),'%'],'FontSize',18,'FontWeight','bold')
 ylabel ([PCs{2},'-',num2str(explained(str2num(PCs{2}(3)))),'%'],'FontSize',18,'FontWeight','bold')

 cl=1;
 yp_all = cell(1,2);
 xp_all = cell(1,2);
 for pl_task = 1:size(DataToShow,2)        
     index = find (info_data(:,3)==pl_task);
     xp = SCORE(index,str2num(PCs{1}(3)));
     yp = SCORE(index,str2num(PCs{2}(3)));
     
%      xp = NaN(size(ListRat,1)-3,1);
%      yp = NaN(size(ListRat,1)-3,1);
%      %animale = [4,5,6,8];
%      for i_rat = 1:size(xp,1)
%         %i_rat = animale(ii);
%         index = find(info_data(:,3)==pl_task & info_data(:,2)== i_rat);
%         xp(i_rat) = nanmean(SCORE(index,str2num(PCs{1}(3))));
%         yp(i_rat) = nanmean(SCORE(index,str2num(PCs{2}(3))));
%      end
     x_r = xp;
     y_r = yp;
     x_r = isoutlier(xp);
     y_r = isoutlier(yp);
     %xp((x_r+y_r)>0.5)=[];
     %yp((x_r+y_r)>0.5)=[];
     scatter(nanmedian(xp),nanmedian(yp),200,'MarkerFaceColor', mycolor{cl},'MarkerEdgeColor','k','DisplayName',DataToShow{pl_task})
     scatter((xp),(yp),20,'MarkerFaceColor', mycolor{cl},'MarkerEdgeColor', mycolor{cl}, 'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.5,'HandleVisibility','off')
     cl=cl+1;
     yp_all {pl_task} = yp;
     xp_all {pl_task} = xp;
 end
 
% % figure and stat on distance PC graph 
% figure
% yp_stat = [];
% xp_stat = [];
% group_stat = [];
% dataPlotx = NaN(1,2);
% errhighx = NaN(1,2);
% dataPloty = NaN(1,2);
% errhighy = NaN(1,2);
%  for pl_day = 1:2
%      yp_stat = [yp_stat;yp_all{pl_day}];
%      xp_stat = [xp_stat;xp_all{pl_day}];
%      group_stat = [group_stat; ones(length(yp_all{pl_day}),1)*pl_day];
%      dataPlotx(pl_day) = nanmean(xp_all{pl_day});
%      errhighx(pl_day) = nanstd(xp_all{pl_day})/sqrt(length(yp_all{pl_day}));
%      dataPloty(pl_day) = nanmean(yp_all{pl_day});
%      errhighy(pl_day) = nanstd(yp_all{pl_day})/sqrt(length(yp_all{pl_day}));
%  end
% xbar = [1:2];
% subplot(1,2,1)
% bb=bar(xbar,dataPlotx,'stacked');
% bb.EdgeColor = 'flat';
% bb.FaceColor ='flat';
% bb.CData(1,:) = mycolor{1};
% bb.CData(2,:)= mycolor{2};
% ylabel ('PC1','FontSize',13,'FontName','Arial')
% hold on
% er = errorbar(xbar,dataPlotx,-3*errhighx,3*errhighx,'HandleVisibility','off');
% er.Color = [0 0 0];                            
% er.LineStyle = 'none';
% er.LineWidth = 1;
% xlim([0 length(xbar)+1])
% xticks(xbar)
% xticklabels(DataToShow);
% 
% subplot(1,2,2)
% bb=bar(xbar,dataPloty,'stacked');
% bb.EdgeColor = 'flat';
% bb.FaceColor ='flat';
% bb.CData(1,:) = mycolor{1};
% bb.CData(2,:)= mycolor{2};
% ylabel ('PC2','FontSize',13,'FontName','Arial')
% hold on
% er = errorbar(xbar,dataPloty,-3*errhighy,3*errhighy,'HandleVisibility','off');
% er.Color = [0 0 0];                            
% er.LineStyle = 'none';
% er.LineWidth = 1;
% xlim([0 length(xbar)+1])
% xticks(xbar)
% xticklabels(DataToShow);
% 
% [cx,~,stat] = kruskalwallis(xp_stat,group_stat);
% [cy,~,stat] = kruskalwallis(yp_stat,group_stat);
%%
 figure(1005)
 hold on
 set(gca,'FontSize',13)
 legend
 xlabel ([PCs{1},'-',num2str(explained(str2num(PCs{1}(3)))),'%'],'FontSize',18,'FontWeight','bold')
 ylabel ([PCs{2},'-',num2str(explained(str2num(PCs{2}(3)))),'%'],'FontSize',18,'FontWeight','bold')

 cl=1;  
  for pl_task = 1: size(DataToShow,2)
     index = find (info_data(:,3)==pl_task);
     xp = SCORE(index,str2num(PCs{1}(3)));
     yp = SCORE(index,str2num(PCs{2}(3)));
     x_r = isoutlier(xp);
     y_r = isoutlier(yp);
     xp((x_r+y_r)>0.5)=[];
     yp((x_r+y_r)>0.5)=[];
                 
     er_x1 = nanstd(xp)/sqrt(length(xp));
     er_y1 = nanstd(yp)/sqrt(length(yp));
          
     x0=nanmedian(xp); % x0,y0 ellipse centre coordinates
     y0=nanmedian(yp);
     t=-pi:0.01:pi;
     xc=x0+2*er_x1*cos(t);
     yc=y0+2*er_y1*sin(t);
     fill(xc,yc,mycolor{cl},'LineStyle', 'none','DisplayName',DataToShow{pl_task})
     alpha(.3)     
     scatter(nanmedian(xp),nanmedian(yp),150,'MarkerFaceColor', mycolor{cl},'MarkerEdgeColor','k','HandleVisibility','off')
     cl=cl+1;
  end

%% Plot PCA in three components
if PLOT3
    figure
    
    cl=1;
    yp_all = cell(1,2);
    xp_all = cell(1,2);
    zp_all = cell(1,2);
    for pl_task = 1:size(DataToShow,2)
        xp = NaN(size(ListRat,1)-3,1);
        yp = NaN(size(ListRat,1)-3,1);
        zp = NaN(size(ListRat,1)-3,1);
        animale = [1,2,3,7];
        for ii = 1:4%size(xp,1)
            i_rat = animale(ii);
            index = find(info_data(:,3)==pl_task & info_data(:,2)== i_rat);
            xp(i_rat) = nanmean(SCORE(index,str2num(PCs{1}(3))));
            yp(i_rat) = nanmean(SCORE(index,str2num(PCs{2}(3))));
            zp(i_rat) = nanmean(SCORE(index,str2num(PCs{3}(3))));
        end
        x_r = isoutlier(xp);
        y_r = isoutlier(yp);
        z_r = isoutlier(zp);
        %xp((x_r+y_r)>0.5)=[];
        %yp((x_r+y_r)>0.5)=[];
        scatter3(nanmedian(xp),nanmedian(yp),nanmedian(zp),200,'MarkerFaceColor', mycolor{cl},'MarkerEdgeColor','k','DisplayName',DataToShow{pl_task})
        hold on
        scatter3((xp),(yp),(zp),20,'MarkerFaceColor', mycolor{cl},'MarkerEdgeColor', mycolor{cl}, 'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.5,'HandleVisibility','off')
        cl=cl+1;
        yp_all {pl_task} = yp;
        xp_all {pl_task} = xp;
        zp_all {pl_task} = zp;
    end
    set(gca,'FontSize',13)
    legend
    explained_R = round(explained);
    xlabel ([PCs{1},'-',num2str(explained_R(str2num(PCs{1}(3)))),'%'],'FontSize',18,'FontWeight','bold')
    ylabel ([PCs{2},'-',num2str(explained_R(str2num(PCs{2}(3)))),'%'],'FontSize',18,'FontWeight','bold')
    zlabel ([PCs{3},'-',num2str(explained_R(str2num(PCs{3}(3)))),'%'],'FontSize',18,'FontWeight','bold')
    
    %Plot distances PCs
figure
zp_stat = [];
yp_stat = [];
xp_stat = [];
pc_stat = [];
group_stat = [];
dataPlotx = NaN(1,2);
errhighx = NaN(1,2);
dataPloty = NaN(1,2);
errhighy = NaN(1,2);
dataPlotz = NaN(1,2);
errhighz = NaN(1,2);
dataPlotPC = NaN(1,2);
errhighPC = NaN(1,2);
pc_statd_s = cell(2,1);
 for pl_day = 1:2
     pc_statd = sqrt((yp_all{pl_day}).^2+(xp_all{pl_day}).^2+...
     +(zp_all{pl_day}).^2);
     pc_statd_s{pl_day} = pc_statd;
     yp_stat = [yp_stat;yp_all{pl_day}];
     xp_stat = [xp_stat;xp_all{pl_day}];
     zp_stat = [zp_stat;zp_all{pl_day}];
     pc_stat = [pc_stat;pc_statd];
     group_stat = [group_stat; ones(length(yp_all{pl_day}),1)*pl_day];
     dataPlotx(pl_day) = nanmean(xp_all{pl_day});
     errhighx(pl_day) = nanstd(xp_all{pl_day})/sqrt(length(yp_all{pl_day}));
     dataPloty(pl_day) = nanmean(yp_all{pl_day});
     errhighy(pl_day) = nanstd(yp_all{pl_day})/sqrt(length(yp_all{pl_day}));
     dataPlotz(pl_day) = nanmean(zp_all{pl_day});
     errhighz(pl_day) = nanstd(zp_all{pl_day})/sqrt(length(zp_all{pl_day}));
     dataPlotPC(pl_day) = nanmean(pc_statd);
     errhighPC(pl_day) = nanstd(pc_statd)/sqrt(length(pc_statd));
 end
xbar = [1:2];
subplot(2,2,2)
bb=bar(xbar,dataPlotx,'stacked');
bb.EdgeColor = 'flat';
bb.FaceColor ='flat';
bb.CData(1,:) = mycolor{1};
bb.CData(2,:)= mycolor{2};
ylabel ('PC1','FontSize',13,'FontName','Arial')
hold on
er = errorbar(xbar,dataPlotx,[0,0],errhighx,'HandleVisibility','off');
er.Color = [0 0 0];                            
er.LineStyle = 'none';
er.LineWidth = 1;
x_point = ones(1,8);
scatter(x_point,(xp_all{1}),20,'k','MarkerFaceColor','k')
scatter(x_point*2,(xp_all{2}),20,'k','MarkerFaceColor','k')
xlim([0 length(xbar)+1])
%xticks(xbar)
%xticklabels(DataToShow);

subplot(2,2,3)
bb=bar(xbar,dataPloty,'stacked');
bb.EdgeColor = 'flat';
bb.FaceColor ='flat';
bb.CData(1,:) = mycolor{1};
bb.CData(2,:)= mycolor{2};
ylabel ('PC2','FontSize',13,'FontName','Arial')
hold on
er = errorbar(xbar,dataPloty,[0,0],errhighy,'HandleVisibility','off');
er.Color = [0 0 0];                            
er.LineStyle = 'none';
er.LineWidth = 1;
scatter(x_point,(yp_all{1}),20,'k','MarkerFaceColor','k')
scatter(x_point*2,(yp_all{2}),20,'k','MarkerFaceColor','k')
xlim([0 length(xbar)+1])
%xticks(xbar)
%xticklabels(DataToShow);

subplot(2,2,4)
bb=bar(xbar,dataPlotz,'stacked');
bb.EdgeColor = 'flat';
bb.FaceColor ='flat';
bb.CData(1,:) = mycolor{1};
bb.CData(2,:)= mycolor{2};
ylabel ('PC3','FontSize',13,'FontName','Arial')
hold on
er = errorbar(xbar,dataPlotz,[0,0],errhighz,'HandleVisibility','off');
er.Color = [0 0 0];                            
er.LineStyle = 'none';
er.LineWidth = 1;
scatter(x_point,(zp_all{1}),20,'k','MarkerFaceColor','k')
scatter(x_point*2,(zp_all{2}),20,'k','MarkerFaceColor','k')
xlim([0 length(xbar)+1])
%xticks(xbar)
%xticklabels(DataToShow);

subplot(2,2,1)
bb=bar(xbar,dataPlotPC,'stacked');
bb.EdgeColor = 'flat';
bb.FaceColor ='flat';
bb.CData(1,:) = mycolor{1};
bb.CData(2,:)= mycolor{2};
ylabel ('PCs','FontSize',13,'FontName','Arial')
hold on
er = errorbar(xbar,dataPlotPC,[0,0],errhighPC,'HandleVisibility','off');
er.Color = [0 0 0];                            
er.LineStyle = 'none';
er.LineWidth = 1;
scatter(x_point,(pc_statd_s{1}),20,'k','MarkerFaceColor','k')
scatter(x_point*2,(pc_statd_s{2}),20,'k','MarkerFaceColor','k')
xlim([0 length(xbar)+1])
%xticks(xbar)
%xticklabels(DataToShow);

[cx,~,stat] = kruskalwallis(xp_stat,group_stat);
[cy,~,stat] = kruskalwallis(yp_stat,group_stat);
[cz,~,stat] = kruskalwallis(zp_stat,group_stat);
[cpc,~,stat] = kruskalwallis(pc_stat,group_stat);

end    
 %% Normalization and mean inside the same animal

if MEDrat
        matrix_all_medrat = NaN(size(ListRat,1)-3,size(NameParam,2),size(DataToShow,2));
       for j_rat = 1:size(matrix_all_medrat,1)
           n_rat =str2num(ListRat(j_rat+2,1).name(2:3));
           for j_task = 1: size(DataToShow,2)
               in_med = find (info_data(:,3)== j_task & info_data(:,2) == n_rat);
               matrix_all_medrat(j_rat,:,j_task) = nanmean(matrix_allparams_allanimals(in_med,:),1);
           end
       end
       MM = matrix_all_medrat;
else
       tot_task = unique(info_data(:,3));
       max_cycles = NaN(length(tot_task),1);
       for j_task = 1 : length(tot_task)
           max_cycles(j_task) = length(find(info_data(:,3) == tot_task(j_task)));
       end
       MM = NaN(max(max_cycles),size(NameParam,2),size(DataToShow,2));
       for j_task = 1 : length(tot_task)
           in_med = find(info_data(:,3) == tot_task(j_task));
           MM(1:max_cycles(j_task),:,j_task) = matrix_allparams_allanimals(in_med,:);
       end
end

  %% Figure and stat for every parameter
 
nt_p = NumParam; % total number of parameter
nt_d = size(DataToShow,2); % different task
    
xbar = [1:1:nt_d];
errhigh = NaN (nt_d,nt_p);
dataPlot = NaN (nt_d,nt_p);
for k_task = 1:nt_d
        dd = MM(:,:,k_task); 
        dd(isnan(dd(:,1)),:)= [];
        errhigh(k_task,:) = nanstd(dd,1)/sqrt(size(dd,1));
        dataPlot(k_task,:) = nanmean(dd,1);
end
% 
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
    xticklabels(DataToShow);
        
    x_stat = []; % in x_stat data for statistic are collected
    group = []; % in group data for statistic are collected
    for jj = 1 :nt_d
        x_stat = [x_stat;MM(:,j_par,jj)];
        group = [group;ones(length(MM(:,j_par,jj)),1)*jj];
    end
     p = kruskalwallis(x_stat,group);
     title(NameParam{j_par});  
end

%% Fig more important parameters

pc = [COEFF(:,str2num(PCs{1}(3))),COEFF(:,str2num(PCs{2}(3))),COEFF(:,str2num(PCs{3}(3)))];
pc = abs(pc);
for i_pc = 1:size(pc,2)
    figure(1012+i_pc)
     hold on
    [~,in] = max(pc(:,i_pc));
    subplot (1,2,1)
    pc(in,i_pc) = 0;
    bb=bar(xbar,dataPlot(:,in),'stacked');
    bb.EdgeColor = 'flat';
    bb.FaceColor ='flat';
    bb.CData(1,:) = mycolor{1};
    bb.CData(2,:)= mycolor{2};
    
    ylabel (NameParam{in},'FontSize',13,'FontName','Arial')
    hold on
    er = errorbar(xbar,dataPlot(:,in),[0,0],errhigh(:,in),'HandleVisibility','off');
    er.Color = [0 0 0];                            
    er.LineStyle = 'none';
    er.LineWidth = 1;
    scatter(x_point,MM(:,in,1),20,'k','MarkerFaceColor','k','MarkerEdgeColor','none')
    scatter(x_point*2,MM(:,in,2),20,'k','MarkerFaceColor','k','MarkerEdgeColor','none')
    xlim([0 length(xbar)+1])
    xticks(xbar)
    xticklabels(DataToShow);
    
    
    [~,in] = max(pc(:,i_pc));
    subplot (1,2,2)
    pc(in,i_pc) = 0;
    bb=bar(xbar,dataPlot(:,in),'stacked');
    bb.EdgeColor = 'flat';
    bb.FaceColor ='flat';
    bb.CData(1,:) = mycolor{1};
    bb.CData(2,:)= mycolor{2};
    
    ylabel (NameParam{in},'FontSize',13,'FontName','Arial')
    hold on
    er = errorbar(xbar,dataPlot(:,in),[0,0],errhigh(:,in),'HandleVisibility','off');
    er.Color = [0 0 0];                            
    er.LineStyle = 'none';
    er.LineWidth = 1;
    scatter(x_point,MM(:,in,1),20,'k','MarkerFaceColor','k','MarkerEdgeColor','none')
    scatter(x_point*2,MM(:,in,2),20,'k','MarkerFaceColor','k','MarkerEdgeColor','none')
    xlim([0 length(xbar)+1])
    xticks(xbar)
    xticklabels(DataToShow);

end

%% Create color map for all parameters for PCs
figure;
imagesc(abs(COEFF(:,1:3)))
colormap('Copper')

yticks(1:32)
yticklabels(NameParam)
xticks(1:3)
xticklabels(PCs)

% %% Figure param PC3 for bipedal and quadrupedal divided
% 
% pc = COEFF(:,str2num(PCs{3}(3)));
% pc = abs(pc);
% an_bi = [1,2,3,7];
% %an_quad = [4,5,6,8];
% x_point = ones(1,4);
%     figure(1021)
%     hold on
%     [~,in1] = max(pc);
%     subplot (1,2,1)
%     pc(in1) = 0;
%     d_plot = NaN(2,1);
%     d_plot(1) = nanmean(MM(an_bi,in1,1));
%     d_plot(2) = nanmean(MM(an_bi,in1,2));
%     bb=bar(xbar,d_plot,'stacked');
%     bb.EdgeColor = 'flat';
%     bb.FaceColor ='flat';
%     bb.CData(1,:) = mycolor{1};
%     bb.CData(2,:)= mycolor{2};
%     
%     ylabel (NameParam{in1},'FontSize',13,'FontName','Arial')
%     hold on
%     er_plot = NaN(2,1);
%     er_plot(1) = nanstd(MM(an_bi,in1,1));
%     er_plot(2) = nanstd(MM(an_bi,in1,2));
%     er = errorbar(xbar,d_plot,[0,0],er_plot,'HandleVisibility','off');
%     er.Color = [0 0 0];                            
%     er.LineStyle = 'none';
%     er.LineWidth = 1;
%     scatter(x_point,MM(an_bi,in1,1),20,'k','MarkerFaceColor','k','MarkerEdgeColor','none')
%     scatter(x_point*2,MM(an_bi,in1,2),20,'k','MarkerFaceColor','k','MarkerEdgeColor','none')
%     xlim([0 length(xbar)+1])
%     xticks(xbar)
%     xticklabels(DataToShow);
%     
%     
%     [~,in2] = max(pc);
%     subplot (1,2,2)
%     pc(in2) = 0;
%     d_plot = NaN(2,1);
%     d_plot(1) = nanmean(MM(an_bi,in2,1));
%     d_plot(2) = nanmean(MM(an_bi,in2,2));
%     bb=bar(xbar,d_plot,'stacked');
%     bb.EdgeColor = 'flat';
%     bb.FaceColor ='flat';
%     bb.CData(1,:) = mycolor{1};
%     bb.CData(2,:)= mycolor{2};
%     
%     ylabel (NameParam{in2},'FontSize',13,'FontName','Arial')
%     hold on
%     er_plot = NaN(2,1);
%     er_plot(1) = nanstd(MM(an_bi,in2,1));
%     er_plot(2) = nanstd(MM(an_bi,in2,2));
%     er = errorbar(xbar,d_plot,[0,0],er_plot,'HandleVisibility','off');
%     er.Color = [0 0 0];                            
%     er.LineStyle = 'none';
%     er.LineWidth = 1;
%     scatter(x_point,MM(an_bi,in2,1),20,'k','MarkerFaceColor','k','MarkerEdgeColor','none')
%     scatter(x_point*2,MM(an_bi,in2,2),20,'k','MarkerFaceColor','k','MarkerEdgeColor','none')
%     xlim([0 length(xbar)+1])
%     xticks(xbar)
%     xticklabels(DataToShow);
% 
