%% Plot parameters and PCA Push active(during pulling phase) along days 
% To run this code it is necessary to have run until x04 for all animals of
% the days to analyse
% 1.  You can choose the days to visualize, if the mean of the parameters 
% is performed on the single animals or all together, if you want to
% normalize and also if you want to save results


clear all
close all

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

NameForce = {'Fz','Fx'};
%select days to plot 
DAY = {'DAY02_2019_08_19','DAY03_2019_08_20','DAY04_2019_08_23','DAY05_2019_08_28','DAY06_2019_09_11','DAY07_2019_09_13','DAY08_2019_09_18'};
    %'DAY09_2019_09_27_4DPI','DAY10_2019_10_01_8DPI','DAY11_2019_10_07_14DPI','DAY12_2019_10_14_21DPI','DAY13_2019_10_29_35DPI'}; 
nBL = 7;
dayPlot = {'Bipedal','Quadrupedal'};
ListBI = [1,2,3,7];
ListQU = [4,5,6,8];
mycolor = {[186,155,201]./255,[109,105,183]./255,[240,76,2]./255,[153,28,0]./255};%[209 211 212]./255,[253 229 205]./255,
PCs = {'PC1','PC2','PC3'}; % select the PC that you want to use for the pca plot
MEDrat = 1; % 1 if you want to make the mean of the animals otherwise 0 
SAVE = 0; % save all the figure for every parameter
PLOT_p = 1; % plot all parameters
REM = 1; % remove outliers PCA
PLOT3 = 1; % plot pca 3D

datapath = 'C:\R-Platform\DATA\2019_02_Group1\Recordings_robot\';
Path = cd;

NumParam = length(NameParam);
LenDataTS = length(DataToShow);
%matrix to save all data (righe = Par) * (colonne = trattamento)
MatrixTot = cell(NumParam, LenDataTS);

vect_allpars = cell (size(DAY,2),size(DataToShow,2));
matrix_allpars = cell(size(DAY,2),size(DataToShow,2)); %one cell for every modality
info_rat = cell(size(DAY,2),size(DataToShow,2));
%info_day = cell(size(DAY,2),size(DataToShow,2));

%% load data

for i_day = 1:size(DAY,2)
    ListRat = dir ([datapath,DAY{1,i_day}]); 
%     if ~contains (DAY{1,i_day},'DPI') % we keep track of the day that we are analysing
%         info_day{i_day,1} = 0; % all the baseline data can be analyzed together, so they are all zero
%     else
%         info_day{i_day,1} = str2num(DAY{1,i_day}(4:5)); % number of the day 
%     end
    for i_rat = 3:size(ListRat,1)-1 % because there is always the file ReadMe in every folder
        ListFile = dir([datapath,DAY{1,i_day},'\',ListRat(i_rat,1).name]);
        in_files = [];
        
        for i_file = 3:length(ListFile)
            if contains(ListFile(i_file,1).name,'x04.mat') 
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
                        line1 = [mean(parametri(index_pks,4:13),1),t_tot,t_first,n_peaks,peak_freq,deltaFcy(i)];
                    else                        
                        line1 = NaN(1,15);
                    end
                    line1(2) = [];
                    % Data EMG BI
                    BIst = round(burstBI(:,1)*Data.Recorded_Data.fS_robot);
                    BIst(BIst==0)=[];% delete of the bursts that starts before the beginnning
                    index = find (Data.Recorded_Data.cicles.data(BIst)==Data.good_trials(i,1));%& Data.Recorded_Data.T_status.data(BIst)==2); %we want to analize only data coming from the pulling phase
                    if ~isempty(index)
                        n_burst = length(index); %total number of burst for each trial
                        line3 = [n_burst,mean(burstBI(index,2:5),1),deltaEMGcy(i,1)];                      
                    else                        
                        line3 = NaN(1,6);
                    end

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

info_data = [];
matrix_allparams_allanimals = [];
for i_raw = 1: size(matrix_allpars,1)
    for i_col = 1: size(matrix_allpars,2)
        matrix_allparams_allanimals = [matrix_allparams_allanimals; matrix_allpars{i_raw,i_col}];
        info_data = [info_data; info_rat{i_raw,i_col}]; %first the number of rat
    end
end

%% Normalization and mean inside the same animal
    %if rat by rat 
if MEDrat
    ListAn = unique(info_data);
    matrix_all_medrat = NaN(length(ListAn),size(NameParam,2));
   for j_rat = 1:size(matrix_all_medrat,1)
       n_rat =ListAn(j_rat); %str2num(ListRat(j_rat+2,1).name(2:3));
       in_med = find (info_data(:,1)== n_rat);
       matrix_all_medrat(j_rat,:) = nanmean(matrix_allparams_allanimals(in_med,:),1);       
   end
   MM = matrix_all_medrat;
else
   max_cycles = NaN(length(tot_days),1);
   for j_day = 1 : length(tot_days)
       max_cycles(j_day) = length(find(info_data(:,1) == tot_days(j_day)));
   end
   MM = NaN(max(max_cycles),size(NameParam,2),size(DAY,2)-nBL+1);
   for j_day = 1 : length(tot_days)
       in_med = find(info_data(:,1) == tot_days(j_day));
       MM(1:max_cycles(j_day),:,j_day) = matrix_allparams_allanimals(in_med,:);
   end
end


%% Figure and stat for every parameter

nt_p = NumParam; % total number of parameter
nt_d = 2; %bipedal vs quadrupedal
xbar = [1,2];
  
errhigh = NaN (nt_d,nt_p);
dataPlot = NaN (nt_d,nt_p);
for k_pos = 1:nt_d
    if k_pos == 1
        dd = MM(ListBI,:);
    else
        dd = MM(ListQU,:);
    end
        dd(isnan(dd(:,1)),:)= [];
        errhigh(k_pos,:) = nanstd(dd,1)/sqrt(size(dd,1));
        dataPlot(k_pos,:) = nanmean(dd,1);
end

if PLOT_p
    
for j_par = 1: nt_p
    figure(j_par+700)    
    hold on
    bb=bar(xbar,dataPlot(:,j_par),'stacked');
    bb.EdgeColor = 'flat';
    bb.FaceColor ='flat';
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
        if jj == 1
            ind = ListBI;
        else
            ind = ListQU;
        end
        x_stat = [x_stat;MM(ind,j_par)];
        group = [group;ones(length(MM(ind,j_par)),1)*jj];
    end
    xticklabels(dayPlot);

    %statistic
%     [p,tbl,stats] = kruskalwallis(x_stat,group);
%     c = multcompare(stats);
    [p,h,stats] = ranksum(x_stat(1:4),x_stat(5:8)) % remember for attempts x_stat(2) is outlier with isoutlier
    title(NameParam{j_par});    
end
end

%% Calculate PCA
zscor_xnan = @(x) bsxfun(@rdivide, bsxfun(@minus, x, mean(x,'omitnan')), std(x, 'omitnan'));
matrix_allparams_allanimals_z =  zscor_xnan(matrix_allparams_allanimals);

cd([Path,'\SubFunctions'])
matrix_allparams_allanimals_z = Reconstruct_matrix_PCA (matrix_allparams_allanimals_z,info_data,REM);

[COEFF, SCORE, LATENT,~,explained] = pca((matrix_allparams_allanimals_z));%,'algorithm','eig','Rows','pairwise');
figure; imagesc(COEFF)

%nanmean(matrix_allparams_allanimals_z,1);
%nanstd(matrix_allparams_allanimals_z,1);
 
 figure(1108)
 hold on
 set(gca,'FontSize',13)
 legend
 xlabel ([PCs{1},'-',num2str(explained(str2num(PCs{1}(3)))),'%'],'FontSize',18,'FontWeight','bold')
 ylabel ([PCs{2},'-',num2str(explained(str2num(PCs{2}(3)))),'%'],'FontSize',18,'FontWeight','bold')

 cl=1;
 yp_all = cell(1,nt_d);
 xp_all = cell(1,nt_d);
xp = NaN(length(ListAn),1);
yp = NaN(length(ListAn),1);
for i_rat = 1:size(xp,1)
    index = find(info_data(:,1)==i_rat);
    xp(i_rat) = nanmean(SCORE(index,str2num(PCs{1}(3))));
    yp(i_rat) = nanmean(SCORE(index,str2num(PCs{2}(3))));
end
vect = NaN(size(info_data,1),1);
for i = 1:length(vect)
    b = info_data(i,1);
    if ismember(b,ListBI)
        vect(i) =1;
    else
        vect(i) =0;
    end
end
 for pl_day = 1:nt_d
     if pl_day==1 % we keep track of the position of the animals
       ind1 = ListBI;
       ind = find(vect==1);
    else
        ind1 = ListQU;
        ind = find(vect==0);
     end    
     scatter(nanmean(xp(ind1)),nanmean(yp(ind1)),200,'MarkerFaceColor', mycolor{cl},'MarkerEdgeColor','k','DisplayName',dayPlot{pl_day})
     %scatter(xp(ind),yp(ind),20,'MarkerFaceColor', mycolor{cl},'MarkerEdgeColor', mycolor{cl}, 'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.5,'HandleVisibility','off')
     scatter(SCORE(ind,str2num(PCs{1}(3))),SCORE(ind,str2num(PCs{2}(3))),20,'MarkerFaceColor', mycolor{cl},'MarkerEdgeColor', mycolor{cl}, 'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.5,'HandleVisibility','off')
     cl=cl+1; 
     yp_all {pl_day} = yp(ind1);
     xp_all {pl_day} = xp(ind1);
 end
 
%% Plot PCA in three components
if PLOT3
    figure
    cl=1;
    yp_all = cell(1,nt_d);
    xp_all = cell(1,nt_d);
    zp_all = cell(1,nt_d);
    xp = NaN(length(ListAn),1);
    yp = NaN(length(ListAn),1);
    zp = NaN(length(ListAn),1);
    for i_rat = 1:size(xp,1)
        index = find(info_data(:,1)==ListAn(i_rat));
        xp(i_rat) = nanmean(SCORE(index,str2num(PCs{1}(3))));
        yp(i_rat) = nanmean(SCORE(index,str2num(PCs{2}(3))));
        zp(i_rat) = nanmean(SCORE(index,str2num(PCs{3}(3))));
    end
    for pl_day = 1:nt_d
        if pl_day==1 % we keep track of the day that we are analysing
            ind = ListBI;
        else
            ind = ListQU;
        end 
        scatter3(nanmean(xp(ind)),nanmean(yp(ind)),nanmean(zp(ind)),200,'MarkerFaceColor', mycolor{cl},'MarkerEdgeColor','k','DisplayName',dayPlot{pl_day})
        hold on
        scatter3((xp(ind)),(yp(ind)),zp(ind),20,'MarkerFaceColor', mycolor{cl},'MarkerEdgeColor', mycolor{cl}, 'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.5,'HandleVisibility','off')
        cl=cl+1;
        zp_all {pl_day} = zp(ind);
        yp_all {pl_day} = yp(ind);
        xp_all {pl_day} = xp(ind);
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
dataPlotx = NaN(1,nt_d);
errhighx = NaN(1,nt_d);
dataPloty = NaN(1,nt_d);
errhighy = NaN(1,nt_d);
dataPlotz = NaN(1,nt_d);
errhighz = NaN(1,nt_d);
dataPlotPC = NaN(1,nt_d);
errhighPC = NaN(1,nt_d);
pc_statd_s = cell(nt_d,1);
 for pl_day = 1:nt_d 
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
xbar = [1:nt_d];
x_point = ones(1,4);
subplot(2,2,2)
bb=bar(xbar,dataPlotx,'stacked');
bb.EdgeColor = 'flat';
bb.FaceColor ='flat';
bb.CData(1,:) = mycolor{1};
bb.CData(2,:)= mycolor{2};
ylabel ('PC1','FontSize',13,'FontName','Arial')
hold on
er = errorbar(xbar,dataPlotx,-errhighx,errhighx,'HandleVisibility','off');
er.Color = [0 0 0];                            
er.LineStyle = 'none';
er.LineWidth = 1;
for i_dd = 1:nt_d
    if i_dd ==1
        x_point = ones(length(ListBI),1);
    else
        x_point = ones(length(ListQU),1);
    end
    scatter(x_point*i_dd,(xp_all{i_dd}),20,'k','MarkerFaceColor','k','MarkerEdgeColor','none')
end
xlim([0 length(xbar)+1])
%xticks(xbar)
%xticklabels(dayPlot);

subplot(2,2,3)
bb=bar(xbar,dataPloty,'stacked');
bb.EdgeColor = 'flat';
bb.FaceColor ='flat';
bb.CData(1,:) = mycolor{1};
bb.CData(2,:)= mycolor{2};
ylabel ('PC2','FontSize',13,'FontName','Arial')
hold on
er = errorbar(xbar,dataPloty,-errhighy,errhighy,'HandleVisibility','off');
er.Color = [0 0 0];                            
er.LineStyle = 'none';
er.LineWidth = 1;
for i_dd = 1:nt_d
    if i_dd ==1
        x_point = ones(length(ListBI),1);
    else
        x_point = ones(length(ListQU),1);
    end
    scatter(x_point*i_dd,(yp_all{i_dd}),20,'k','MarkerFaceColor','k','MarkerEdgeColor','none')
end
xlim([0 length(xbar)+1])
%xticks(xbar)
%xticklabels(dayPlot);

subplot(2,2,4)
bb=bar(xbar,dataPlotz,'stacked');
bb.EdgeColor = 'flat';
bb.FaceColor ='flat';
bb.CData(1,:) = mycolor{1};
bb.CData(2,:)= mycolor{2};
ylabel ('PC3','FontSize',13,'FontName','Arial')
hold on
er = errorbar(xbar,dataPlotz,-errhighz,errhighz,'HandleVisibility','off');
er.Color = [0 0 0];                            
er.LineStyle = 'none';
er.LineWidth = 1;
for i_dd = 1:nt_d
    if i_dd ==1
        x_point = ones(length(ListBI),1);
    else
        x_point = ones(length(ListQU),1);
    end
    scatter(x_point*i_dd,(zp_all{i_dd}),20,'k','MarkerFaceColor','k','MarkerEdgeColor','none')
end
xlim([0 length(xbar)+1])
%xticks(xbar)
%xticklabels(dayPlot);

subplot(2,2,1)
bb=bar(xbar,dataPlotPC,'stacked');
bb.EdgeColor = 'flat';
bb.FaceColor ='flat';
bb.CData(1,:) = mycolor{1};
bb.CData(2,:)= mycolor{2};
ylabel ('PCs','FontSize',13,'FontName','Arial')
hold on
er = errorbar(xbar,dataPlotPC,-errhighPC,errhighPC,'HandleVisibility','off');
er.Color = [0 0 0];                            
er.LineStyle = 'none';
er.LineWidth = 1;
for i_dd = 1:nt_d
    if i_dd ==1
        x_point = ones(length(ListBI),1);
    else
        x_point = ones(length(ListQU),1);
    end
    scatter(x_point*i_dd,(pc_statd_s{i_dd}),20,'k','MarkerFaceColor','k','MarkerEdgeColor','none')
end
xlim([0 length(xbar)+1])
%xticks(xbar)
%xticklabels(dayPlot);

% [p,~,stat] = kruskalwallis(xp_stat,group_stat);
% cx = multcompare(stat);
% [p,~,stat] = kruskalwallis(yp_stat,group_stat);
% cy = multcompare(stat);
% [p,~,stat] = kruskalwallis(zp_stat,group_stat);
% cz = multcompare(stat);
% [p,~,stat] = kruskalwallis(pc_stat,group_stat);
% cpc = multcompare(stat);

[px,h,stats] = ranksum(xp_stat(1:4),xp_stat(5:8));
[py,h,stats] = ranksum(yp_stat(1:4),yp_stat(5:8));
[pz,h,stats] = ranksum(zp_stat(1:4),zp_stat(5:8));
[pp,h,stats] = ranksum(pc_stat(1:4),pc_stat(5:8));

end    

%% Fig more important parameters

pc = [COEFF(:,str2num(PCs{1}(3))),COEFF(:,str2num(PCs{2}(3))),COEFF(:,str2num(PCs{3}(3)))];
pc = abs(pc);
xbar = [1:nt_d];
for i_pc = 1:size(pc,2)
    figure(1080+i_pc)
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
    er = errorbar(xbar,dataPlot(:,in),-errhigh(:,in),errhigh(:,in),'HandleVisibility','off');
    er.Color = [0 0 0];                            
    er.LineStyle = 'none';
    er.LineWidth = 1;
    x_point = ones(length(ListBI),1);
    scatter(x_point*1,MM(ListBI,in),20,'k','MarkerFaceColor','k','MarkerEdgeColor','none')
    x_point = ones(length(ListQU),1);
    scatter(x_point*2,MM(ListQU,in),20,'k','MarkerFaceColor','k','MarkerEdgeColor','none')
    xlim([0 length(xbar)+1])
    xticks(xbar)
    xticklabels(dayPlot);
    
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
    er = errorbar(xbar,dataPlot(:,in),-errhigh(:,in),errhigh(:,in),'HandleVisibility','off');
    er.Color = [0 0 0];                            
    er.LineStyle = 'none';
    er.LineWidth = 1;
    x_point = ones(length(ListBI),1);
    scatter(x_point*1,MM(ListBI,in),20,'k','MarkerFaceColor','k','MarkerEdgeColor','none')
    x_point = ones(length(ListQU),1);
    scatter(x_point*2,MM(ListQU,in),20,'k','MarkerFaceColor','k','MarkerEdgeColor','none')
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

