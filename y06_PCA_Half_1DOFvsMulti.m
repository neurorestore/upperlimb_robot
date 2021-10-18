%% Plot parameters and PCA Push active(during pulling phase) along days to compare 1DOF vs multi DOF
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
DAY = {'DAY02_2021_03_02_BL','DAY03_2021_03_03_BL','DAY04_2021_03_04_BL','DAY05_2021_03_05_BL','DAY08_2021_03_08_BL'}; 
dayPlot = {'1DOF','MultiDOF'};
mycolor = {[186,155,201]./255,[240,76,2]./255,[153,28,0]./255};
PCs = {'PC1','PC2','PC3'}; % select the PC that you want to use for the pca plot
SAVE = 0; % save all the figure for every parameter
PLOT_p = 1; % plot all parameters
REM = 1; % remove outliers PCA
TRIALS = 0; % 1 if stat and bar with all trials, 0 if mean of the day 
n_rat = 2;
%PLOT3 = 1; % plot pca 3D

datapath = 'C:\R-Platform\DATA\2021_02_Robot Stim group\Robot Data\';
Path = cd;

NumParam = length(NameParam);
LenDataTS = length(DataToShow);
%matrix to save all data (righe = Par) * (colonne = trattamento)
MatrixTot = cell(NumParam, LenDataTS);

vect_allpars = cell (size(DAY,2),size(DataToShow,2));
matrix_allpars = cell(size(DAY,2),size(DataToShow,2)); %one cell for every modality
info_rat = cell(size(DAY,2),size(DataToShow,2));
info_task = cell(size(DAY,2),size(DataToShow,2)); % 1= 1 DOF, 4 = Multi DOF

%% load data

for i_day = 1:size(DAY,2)
    ListRat = dir ([datapath,DAY{1,i_day}]); 
    for i_rat = 3:size(ListRat,1)-1 % because there is always the file ReadMe in every folder
        cd([datapath,DAY{1,i_day},'\',ListRat(i_rat,1).name]);
        [filename, pathname] = uigetfile({'*.mat','Matlab AnalysisPeaks (*.mat)'}, 'Choose Files:','MultiSelect', 'on'); % Select File
        while ischar(filename)       
            load ([pathname filename]);
                        
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
                       if contains(filename,'1DOF')
                           info_task{i_day,1} = [info_task{i_day};1];
                       else
                           info_task{i_day,1} = [info_task{i_day};4];
                       end
                end

            end
            [filename, pathname] = uigetfile({'*.mat','Matlab AnalysisPeaks (*.mat)'}, 'Choose Files:','MultiSelect', 'on'); % Select File

        end
        
    end    
    
end

%% JOIN all data (all animals, all days) in the same matrix to do the PCA

info_data = [];
matrix_allparams_allanimals = [];
for i_raw = 1: size(matrix_allpars,1)
    for i_col = 1: size(matrix_allpars,2)
        matrix_allparams_allanimals = [matrix_allparams_allanimals; matrix_allpars{i_raw,i_col}];
        info_data = [info_data; info_rat{i_raw,i_col},info_task{i_raw,i_col},ones(length(info_rat{i_raw,i_col}),1)*i_raw]; %first the number of rat, then the type of task, then the day
    end
end

%% Figure and stat for every parameter

nt_p = NumParam; % total number of parameter
nt_d = [1,4]; %single vs multi DOF
nt_g = unique (info_data(:,3));
an_num = unique(info_data(:,1));
xbar = [1,2];

n = 0;
% Option 1 -> plot all trials
if TRIALS
    errhigh = NaN (length(nt_d),nt_p);
    dataPlot = NaN (length(nt_d),nt_p);
    for k_pos = nt_d
        n = n+1;
        dd = matrix_allparams_allanimals(info_data(:,2)==k_pos,:);
        dd(isnan(dd(:,1)),:)= [];
        errhigh(n,:) = nanstd(dd,1)/sqrt(size(dd,1));
        dataPlot(n,:) = nanmean(dd,1);
    end
else %Option 2 -> plot mean animal/day
    errhigh = NaN (length(nt_d),nt_p,length(nt_g)*n_rat);
    dataPlot = NaN (length(nt_d),nt_p,length(nt_g)*n_rat);
    for k_pos = nt_d
        n = n+1;
        for k_gg =1:length(nt_g)
            for k_an = 1:n_rat
                dd = matrix_allparams_allanimals(info_data(:,1) == an_num(k_an) & info_data(:,2)==k_pos & info_data(:,3) ==k_gg,:);
                if ~isempty(dd)
                    dd(isnan(dd(:,1)),:)= [];
                    errhigh(n,:,k_gg+(k_an-1)*length(nt_g)) = nanstd(dd,1)/sqrt(size(dd,1));
                    dataPlot(n,:,k_gg+(k_an-1)*length(nt_g)) = nanmean(dd,1);
                end
            end
        end
    end
end

if PLOT_p    
    for j_par = 1: nt_p
        x_stat = []; % in x_stat data for statistic are collected
        group = []; % in group data for statistic are collected
        figure(j_par+700)    
        hold on        
        if TRIALS
            bb=bar(xbar,dataPlot(:,j_par),'stacked');
            er = errorbar(xbar,dataPlot(:,j_par),-errhigh(:,j_par),errhigh(:,j_par),'HandleVisibility','off');
            for jj = nt_d
                x_stat = [x_stat;matrix_allparams_allanimals(info_data(:,2)==jj,j_par)];
                group = [group;ones(length(matrix_allparams_allanimals(info_data(:,2)==jj,j_par)),1)*jj];
            end
        else
            bb=bar(xbar,nanmean(dataPlot(:,j_par,:),3),'stacked');
            er = errorbar(xbar,nanmean(dataPlot(:,j_par,:),3),-nanmean(errhigh(:,j_par,:),3),nanmean(errhigh(:,j_par,:),3),'HandleVisibility','off');
            n = 1;
            for jj = nt_d
                x_stat = [x_stat;permute(dataPlot(n,j_par,:),[3,1,2])];
                group = [group;ones(length(dataPlot(n,j_par,:)),1)*jj];
                n = n+1;
            end
        end        
        bb.EdgeColor = 'flat';
        bb.FaceColor ='flat';
        ylabel (NameParam{j_par},'FontSize',13,'FontName','Arial')
        hold on
        er.Color = [0 0 0];                            
        er.LineStyle = 'none';
        er.LineWidth = 1;
        xlim([0 length(xbar)+1])
        xticks(xbar)
        xticklabels(dayPlot);

        %statistic
        [p,tbl,stats] = kruskalwallis(x_stat,group);
        c = multcompare(stats);
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
yp_all = cell(1,length(nt_d));
xp_all = cell(1,length(nt_d));
nn = 1;
xp_stat = [];
yp_stat = [];
group_stat = [];
 for pl_day = nt_d
     ind = find(info_data(:,2) ==pl_day);
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
errscore = NaN (length(nt_d),length(nt_g)*n_rat,2);
scorePlot = NaN (length(nt_d),length(nt_g)*n_rat,2);
n = 0;
for k_pos = nt_d
    n = n+1;
    for k_gg =1:length(nt_g)
        for k_an = 1:n_rat
            dd = SCORE(info_data(:,1)==an_num(k_an) & info_data(:,2)==k_pos & info_data(:,3) ==k_gg,1:2);
            if ~isempty(dd)
                dd(isnan(dd(:,1)),:)= [];
                errscore(n,k_gg+(k_an-1)*length(nt_g),:) = nanstd(dd,1)/sqrt(size(dd,1));
                scorePlot(n,k_gg+(k_an-1)*length(nt_g),:) = nanmean(dd,1);
            end
        end
    end
end
    
xbar = [1:2];
for i_pc = 1:2
    figure
    hold on
    bb=bar(xbar,nanmean(scorePlot(:,:,i_pc),2),'stacked');
    bb.EdgeColor = 'flat';
    bb.FaceColor ='flat';
    bb.CData(1,:) = mycolor{1};
    bb.CData(2,:)= mycolor{2};
    ylabel (PCs{i_pc} ,'FontSize',13,'FontName','Arial')
    hold on
    er = errorbar(xbar,nanmean(scorePlot(:,:,i_pc),2),-nanmean(errscore(:,:,i_pc),2),nanmean(errscore(:,:,i_pc),2),'HandleVisibility','off');
    er.Color = [0 0 0];                            
    er.LineStyle = 'none';
    er.LineWidth = 1;
    scorestat = [];
    scoregroup = [];
    for i_dd = 1:length(nt_d)
        x_point = ones(size(scorePlot(:,:,i_pc),2),1)*i_dd;
        scatter(x_point,scorePlot(i_dd,:,i_pc),20,'k','MarkerFaceColor','k','MarkerEdgeColor','none')
        scorestat = [scorestat,scorePlot(i_dd,:,i_pc)];
        scoregroup = [scoregroup;ones(length(scorePlot(i_dd,:,i_pc)),1)*i_dd];  
    end
    [p,~,stat] = kruskalwallis(scorestat,scoregroup);
    xlim([0 length(xbar)+1])
end

%% Fig more important parameters

pc = [COEFF(:,str2num(PCs{1}(3))),COEFF(:,str2num(PCs{2}(3)))]; %,COEFF(:,str2num(PCs{3}(3)))];
pc = abs(pc);
for i_pc = 1:size(pc,2)
    figure(1080+i_pc)     
    [~,in] = max(pc(:,i_pc));
    subplot (1,2,1)
    hold on
    pc(in,i_pc) = 0;
    if TRIALS
        bb=bar(xbar,dataPlot(:,in),'stacked');
        er = errorbar(xbar,dataPlot(:,in),-errhigh(:,in),errhigh(:,in),'HandleVisibility','off');
        for gr =1:2
            y_point = matrix_allparams_allanimals(info_data(:,2)==nt_d(gr),in);
            scatter(x_point{gr},y_point,20,'k','MarkerFaceColor','k','MarkerEdgeColor','none')
        end
    else
        bb=bar(xbar,nanmean(dataPlot(:,in,:),3),'stacked');
        er = errorbar(xbar,nanmean(dataPlot(:,in,:),3),-nanmean(errhigh(:,in,:),3),nanmean(errhigh(:,in,:),3),'HandleVisibility','off');
        for gr =1:2
            x_point = ones(size(dataPlot,3),1)*gr;
            y_point = dataPlot(gr,in,:);
            scatter(x_point,y_point,20,'k','MarkerFaceColor','k','MarkerEdgeColor','none')
        end
    end
    bb.EdgeColor = 'flat';
    bb.FaceColor ='flat';
    bb.CData(1,:) = mycolor{1};
    bb.CData(2,:)= mycolor{2};
    
    ylabel (NameParam{in},'FontSize',13,'FontName','Arial')
    er.Color = [0 0 0];                            
    er.LineStyle = 'none';
    er.LineWidth = 1;
    xlim([0 length(xbar)+1])
    xticks(xbar)
    xticklabels(dayPlot);
    
    [~,in] = max(pc(:,i_pc));
    subplot (1,2,2)
    hold on
    pc(in,i_pc) = 0;
    if TRIALS
        bb=bar(xbar,dataPlot(:,in),'stacked');
        er = errorbar(xbar,dataPlot(:,in),-errhigh(:,in),errhigh(:,in),'HandleVisibility','off');
        for gr =1:2
            y_point = matrix_allparams_allanimals(info_data(:,2)==nt_d(gr),in);
            scatter(x_point{gr},y_point,20,'k','MarkerFaceColor','k','MarkerEdgeColor','none')
        end
    else
        bb=bar(xbar,nanmean(dataPlot(:,in,:),3),'stacked');
        er = errorbar(xbar,nanmean(dataPlot(:,in,:),3),-nanmean(errhigh(:,in,:),3),nanmean(errhigh(:,in,:),3),'HandleVisibility','off');
        for gr =1:2
            x_point = ones(size(dataPlot,3),1)*gr;
            y_point = dataPlot(gr,in,:);
            scatter(x_point,y_point,20,'k','MarkerFaceColor','k','MarkerEdgeColor','none')
        end
    end
    bb.EdgeColor = 'flat';
    bb.FaceColor ='flat';
    bb.CData(1,:) = mycolor{1};
    bb.CData(2,:)= mycolor{2};
    
    ylabel (NameParam{in},'FontSize',13,'FontName','Arial')
    er.Color = [0 0 0];                            
    er.LineStyle = 'none';
    er.LineWidth = 1;
    
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