%% Plot parameters and PCA Push active(during pulling phase) along days 
%% SPONTANEOUS vs ROB
% To run this code it is necessary to have run until x04 for all animals of
% the days to analyse
% 1.  You can choose the days to visualize, if the mean of the parameters 
% is performed on the single animals or all together, if you want to
% normalize and also if you want to save results


clear all
close all

n_treat = 2; % number of different treatments/groups of animals

%name of the parameters that you want to plot and to use for the pca
NameParam = {'ForcePeaks','Width','Fmax/Fonset',...%'fwhm',
    'Fmax/fwhm','AUCfwhm','AUCmin','derFup','derFdown',...
    'SmoothF','Tpull','Latency','Npeaks',...
    'Frpeaks','MeanForce'... %'peak_amp_attem','fwhm_attem','width_attem','ratio_amp_attem',...'ratio_pk_fwhm_attem','AUC_fwhm_attem','AUC_min_attem','df_up_attem','df_down_attem','smoothness_attem',
    'Nburst','Ampburst','Maxburst','Tburst','Aucburst','MeanEnv'...
    'sub-movements','attempts','Vpeaks',...
    'height','amplitude','LenTra','AUC',...
    'SpeedTra','SmoothTra','MaxSpeedTra','Acc','MaxAcc'};
    

NameForce = {'Fz'};
%select days to plot (IMPORTANT DPI in name of files recorded after injury, all BL days at the beginning!!!)
DAY = {'DAY01_BL_2020_07_27','DAY02_Wk2_2020_08_11_14DPI','DAY03_Wk4_2020_08_28_35DPI'};%,'DAY04_Wk8_2020_09_23_56DPI'}; %,'DAY04_Wk8_2020_09_23_56DPI'
n_days = size(DAY,2);

dayPlot = {'BL','2 Weeks','4 Weeks'}; 
mycolor = {[209 211 212]./255,[253,227,186]./255,[255,192,139]./255,... % Color reha %[255,100,87]./255,
    [136,141,145]./255,[185,202,233]./255,[60,89,164]./255,... % Color Rob %[40,65,120]./255,
    [63,66,67]./255,[206,232,215]./255,[107,187,136]./255}; % Color Rob+reha %,[57,137,86]./255

PCs = {'PC1','PC2','PC3'}; % select the PC that you want to use for the pca plot
NORM = 0; % 1 if you want normalize on BL otherwise 0
MEDrat = 1; % 1 if you want to make the mean of the animals otherwise 0 
SAVE = 0; % save all the figure for every parameter
PLOT_p = 0; % plot all parameters
REM = 1; % remove outliers PCA
PLOT3 = 1; % plot pca 3D

datapath_ROB = 'C:\R-Platform\DATA\2020_08_RRR_group\Robot Data\';
datapath_gr = 'C:\R-Platform\DATA\2020_08_RRR_group\Code_yRRRgroup\Groups\';
filenames_gr = {'Rehab_group.txt','Rob_group.txt'};%,'Rob+Rehab_group.txt'};
Path = cd;

NumParam = length(NameParam);
%matrix to save all data (righe = Par) * (colonne = trattamento)
MatrixTot = cell(NumParam, 1);

vect_allpars = cell (size(DAY,2)*n_treat,1);
matrix_allpars = cell(size(DAY,2)*n_treat,1); 
info_rat = cell(size(DAY,2)*n_treat,1); %the number of each group is +100
info_day = cell(size(DAY,2)*n_treat,1);

%% load data
for i_gr = 1:n_treat
    fpt = fopen([datapath_gr,filenames_gr{i_gr}]);  % Open programming students grade file
    namesROB = fscanf(fpt, '%s',[1 inf]); %

    for i_day = 1:n_days
        ListRat = dir ([datapath_ROB,DAY{1,i_day}]); 
        if ~contains (DAY{1,i_day},'DPI') % we keep track of the day that we are analysing
            info_day{i_day+(i_gr-1)*n_days,1} = 0; % all the baseline data can be analyzed together, so they are all zero
        else
            info_day{i_day+(i_gr-1)*n_days,1} = str2num(DAY{1,i_day}(end-4:end-3)); % number of the day 
        end
        for i_rat = 3:size(ListRat,1)-1 % because there is always the file ReadMe in every folder
            % check rat in the list               
            if contains (namesROB,ListRat(i_rat,1).name)
                ListFile = dir([datapath_ROB,DAY{1,i_day},'\',ListRat(i_rat,1).name]);
                in_files = [];

                for i_file = 3:length(ListFile)
                    if contains(ListFile(i_file,1).name,'x04.mat') 
                        in_files = [in_files,i_file];
                    end             
                end

                for i_file_ok = 1:length(in_files)
                    load ([datapath_ROB,DAY{1,i_day},'\',ListRat(i_rat,1).name,'\',ListFile(in_files(i_file_ok),1).name]);
                    if contains(Data.info.Status,'Push_active')

                        parametri = Data.Recorded_Data.Analysis.Fzpeaks;
                        deltaFcy = Data.Recorded_Data.Analysis.Area_cy;
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

        %                     % Data EMG BI
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
                                   matrix_allpars{i_day+(i_gr-1)*n_days,1}= [matrix_allpars{i_day+(i_gr-1)*n_days,1};line1,line3,line5];
                                   info_rat{i_day+(i_gr-1)*n_days,1} = [info_rat{i_day+(i_gr-1)*n_days};str2num(ListRat(i_rat,1).name(2:3))+(i_gr-1)*100]; %to distinguish the animals of the second group we had 100
                            end

                        end

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
        info_data = [info_data; ones(size(matrix_allpars{i_raw,i_col},1),1)*info_day{i_raw,i_col}, info_rat{i_raw,i_col}]; %first column for the number of day, second one for the number of rat
    end
end

%% Normalization and mean inside the same animal
ListRat = unique(info_data(:,2));
if NORM
    %create a matrix: an animal for every raw, a parameter for every column
    matr_to_norm = NaN(length(ListRat),size(NameParam,2));
    for j_rat = 1:size(matr_to_norm,1)
        n_rat = ListRat(j_rat); %number of the animal, saved in the info
        in_norm = find (info_data(:,1)== 0 & info_data(:,2) == n_rat); % 0 -> days Baseline
        matr_to_norm(j_rat,:) = nanmean(matrix_allparams_allanimals(in_norm,:),1);
    end
    
    %if rat by rat I normalize after the mean of the animal
    if MEDrat
        matrix_all_medrat = NaN(length(ListRat),size(NameParam,2),n_days);
       for j_rat = 1:size(matrix_all_medrat,1)
           n_rat = ListRat(j_rat);
           for j_day = 1 : n_days
               if ~contains (DAY{1,j_day},'DPI') % we keep track of the day that we are analysing
                   n_day = 0; % all the baseline data can be analyzed together, so they are all zero
               else
                   n_day = str2num(DAY{1,j_day}(end-4:end-3)); % number of the day
               end
               in_med = find (info_data(:,1)== n_day & info_data(:,2) == n_rat);
               matrix_all_medrat(j_rat,:,j_day) = nanmean(matrix_allparams_allanimals(in_med,:),1);
           end
       end
       MM = matrix_all_medrat./matr_to_norm;
    else
        matrix_allparams_allanimals_norm = matrix_allparams_allanimals;
       for j_rat = 1:length(ListRat)
           n_rat = ListRat(j_rat);
           in_med = find (info_data(:,2) == n_rat);
           matrix_allparams_allanimals_norm(in_med,:)= matrix_allparams_allanimals(in_med,:)./matr_to_norm(j_rat,:);
       end
       tot_days = unique(info_data(:,1));
       max_cycles = NaN(n_days,1);
       for j_day = 1 : n_days
           max_cycles(j_day) = length(find(info_data(:,1) == tot_days(j_day)));
       end
       MM = NaN(max(max_cycles),size(NameParam,2),n_days);
       for j_day = 1 : n_days
           in_med = find(info_data(:,1) == tot_days(j_day));
           MM(1:max_cycles(j_day),:,j_day) = matrix_allparams_allanimals_norm(in_med,:);
       end
    end
else %if NORM
    %if rat by rat 
    if MEDrat
        matrix_all_medrat = NaN(length(ListRat),size(NameParam,2),n_days);
       for j_rat = 1:size(matrix_all_medrat,1)
           n_rat = ListRat(j_rat);
           for j_day = 1: n_days
               if ~contains (DAY{1,j_day},'DPI') % we keep track of the day that we are analysing
                   n_day = 0; % all the baseline data can be analyzed together, so they are all zero
               else
                   n_day = str2num(DAY{1,j_day}(end-4:end-3)); % number of the day
               end
               in_med = find (info_data(:,1)== n_day & info_data(:,2) == n_rat);
               matrix_all_medrat(j_rat,:,j_day) = nanmean(matrix_allparams_allanimals(in_med,:),1);
           end
       end
       MM = matrix_all_medrat;
    else
       tot_days = unique(info_data(:,1));
       max_cycles = NaN(n_days,1);
       for j_day = 1 : n_days
           max_cycles(j_day) = length(find(info_data(:,1) == tot_days(j_day)));
       end
       MM = NaN(max(max_cycles),size(NameParam,2),n_days);
       for j_day = 1 : n_days
           in_med = find(info_data(:,1) == tot_days(j_day));
           MM(1:max_cycles(j_day),:,j_day) = matrix_allparams_allanimals(in_med,:);
       end
    end
end

%% Figure and stat for every parameter

nt_p = NumParam; % total number of parameter
nt_d = n_days; %total number of days (baseline all together)
n_gr = [sum(ListRat<100); sum(ListRat>100 & ListRat<200);sum(ListRat>200 & ListRat<300)];
xbar = [1:1:nt_d];
  
errhigh = NaN (nt_d,nt_p,n_treat);
dataPlot = NaN (nt_d,nt_p,n_treat);
for k_day = 1:nt_d
        dd = MM(:,:,k_day); 
        %dd(isnan(dd(:,1)),:)= [];
        for k_gr = 1:n_treat
            en = sum(n_gr(1:k_gr));
            st = en-n_gr(k_gr)+1;
            errhigh(k_day,:,k_gr) = nanstd(dd(st:en,:),1)/sqrt(size(dd(st:en,:),1));
            dataPlot(k_day,:,k_gr) = nanmean(dd(st:en,:),1);
        end
end

if PLOT_p    
    for j_par = 1: nt_p

        figure(j_par+700)    
        hold on
        datapl = permute(dataPlot(:,j_par,:),[1,3,2]);
        err = permute (errhigh(:,j_par,:),[1,3,2]);
        bb=bar(datapl);
        ylabel (NameParam{j_par},'FontSize',13,'FontName','Arial')
        hold on
        % Find the number of groups and the number of bars in each group
        ngroups = size(datapl, 1);
        nbars = size(datapl, 2);
        % Calculate the width for each bar group
        groupwidth = min(0.8, nbars/(nbars + 1.5));
        % Set the position of each error bar in the centre of the main bar
        % Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
        for i = 1:nbars
            % Calculate center of each bar
            x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
            errorbar(x, datapl(:,i), err(:,i), 'k', 'linestyle', 'none');
        end
        hold off
        xlim([0 length(xbar)+1])
        xticks(xbar)
        xticklabels(dayPlot);
        title(NameParam{j_par});

        if SAVE
            cd([datapath,'Results'])
            if NORM
                name2save =[NameParam{j_par},'_norm'];
            else
                name2save =[NameParam{j_par}];
            end
            saveas(gcf,name2save)
        end
        x_stat = [];
        group = [];
        if ~NORM
            for jj = 1 :nt_d
                x_stat = [x_stat;MM(1:5,j_par,jj)];
                group = [group;ones(5,1)*jj];
            end
            %statistic
            [p,tbl,stats] = kruskalwallis(x_stat,group);
            c1 = multcompare(stats);
            title(NameParam{j_par});
            x_stat = [];
            group = [];
            for jj = 1 :nt_d
                x_stat = [x_stat;MM(6:10,j_par,jj)];
                group = [group;ones(5,1)*jj];
            end
            %statistic
            [p,tbl,stats] = kruskalwallis(x_stat,group);
            c2 = multcompare(stats);
            title(NameParam{j_par});
%             x_stat = [];
%             group = [];
%             for jj = 1 :nt_d
%                 x_stat = [x_stat;MM(11:end,j_par,jj)];
%                 group = [group;ones(5,1)*jj];
%             end
%             %statistic
%             [p,tbl,stats] = kruskalwallis(x_stat,group);
%             c3 = multcompare(stats);
            
        end

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

 yp_all = cell(1,nt_d);
 xp_all = cell(1,nt_d);
 gr_all = cell(1,nt_d);
%  for pl_day = 1:nt_d
%      if pl_day==1 % we keep track of the day that we are analysing
%            n_day = 0; % all the baseline data can be analyzed together, so they are all zero
%      else
%            n_day = str2num(DAY{1,pl_day}(end-4:end-3)); % number of the day
%      end
%     xp = NaN(length(ListRat),1);
%     yp = NaN(length(ListRat),1);
%     gr = NaN(length(ListRat),1);
%     for i_rat = 1:size(xp,1)
%         animale = ListRat(i_rat);
%         index = find(info_data(:,1)==n_day & info_data(:,2)== animale);
%         xp(i_rat) = nanmean(SCORE(index,str2num(PCs{1}(3))));
%         yp(i_rat) = nanmean(SCORE(index,str2num(PCs{2}(3))));
%         if animale<100
%             gr(i_rat)=1;
%         elseif animale>100 && animale<200
%             gr(i_rat)=2;
%         else
%             gr(i_rat)=3;
%         end
%     end
%      x_r = xp;
%      y_r = yp;
%      scatter(nanmean(x_r(gr==1)),nanmean(y_r(gr==1)),200,'MarkerFaceColor', mycolor{pl_day},'MarkerEdgeColor','k','DisplayName',[dayPlot{pl_day},'-REHA'])
%      scatter((x_r(gr==1)),(y_r(gr==1)),20,'MarkerFaceColor', mycolor{pl_day},'MarkerEdgeColor', mycolor{pl_day}, 'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.5,'HandleVisibility','off')
%      scatter(nanmean(x_r(gr==2)),nanmean(y_r(gr==2)),200,'MarkerFaceColor', mycolor{pl_day+nt_d},'MarkerEdgeColor','k','DisplayName',[dayPlot{pl_day},'-ROB'])
%      scatter((x_r(gr==2)),(y_r(gr==2)),20,'MarkerFaceColor', mycolor{pl_day+nt_d},'MarkerEdgeColor', mycolor{pl_day+nt_d}, 'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.5,'HandleVisibility','off')
%      scatter(nanmean(x_r(gr==3)),nanmean(y_r(gr==3)),200,'MarkerFaceColor', mycolor{pl_day+2*nt_d},'MarkerEdgeColor','k','DisplayName',[dayPlot{pl_day},'-REHA+ROB'])
%      scatter((x_r(gr==3)),(y_r(gr==3)),20,'MarkerFaceColor', mycolor{pl_day+2*nt_d},'MarkerEdgeColor', mycolor{pl_day+2*nt_d}, 'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.5,'HandleVisibility','off')
%      yp_all {pl_day} = yp;
%      xp_all {pl_day} = xp;
%      gr_all {pl_day} = gr;
%  end
 % Plot with dots for trial
for pl_day = 1:nt_d
     if pl_day==1 % we keep track of the day that we are analysing
           n_day = 0; % all the baseline data can be analyzed together, so they are all zero
     else
           n_day = str2num(DAY{1,pl_day}(end-4:end-3)); % number of the day
     end
       for i_gr = 1:n_treat
           if i_gr ==1
               index = find(info_data(:,1)==n_day & info_data(:,2)< 100);
               cl = pl_day;
           elseif i_gr ==2
               index = find(info_data(:,1)==n_day & info_data(:,2)> 100& info_data(:,2)< 200);
               cl = pl_day+nt_d;
           else
               index = find(info_data(:,1)==n_day & info_data(:,2)> 200);
               cl = pl_day+nt_d;
           end
             xp = SCORE(index,str2num(PCs{1}(3)));
             yp = SCORE(index,str2num(PCs{2}(3)));
             x_r = xp;
             y_r = yp;
         scatter(nanmean(x_r),nanmean(y_r),200,'MarkerFaceColor', mycolor{cl},'MarkerEdgeColor','k')
         scatter((x_r),(y_r),20,'MarkerFaceColor', mycolor{cl},'MarkerEdgeColor', mycolor{cl}, 'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.5,'HandleVisibility','off')
       end
 end
 if ~NORM
     if SAVE
            cd([datapath,'Results'])
            saveas(gcf,'PCA')
     end
 end

%% Plot PCA in three components
if PLOT3
    figure
    
    yp_all = cell(1,nt_d);
    xp_all = cell(1,nt_d);
    zp_all = cell(1,nt_d);
    gr_all = cell(1,nt_d);
    for pl_day = 1:nt_d
        if pl_day==1 % we keep track of the day that we are analysing
            n_day = 0; % all the baseline data can be analyzed together, so they are all zero
        else
          	n_day = str2num(DAY{1,pl_day}(end-4:end-3));  % number of the day
        end
        xp = NaN(length(ListRat),1);
        yp = NaN(length(ListRat),1);
        zp = NaN(length(ListRat),1);
        gr = NaN(length(ListRat),1);
        for i_rat = 1:size(xp,1)
            animale = ListRat(i_rat);
            index = find(info_data(:,1)==n_day & info_data(:,2)== animale);
            xp(i_rat) = nanmean(SCORE(index,str2num(PCs{1}(3))));
            yp(i_rat) = nanmean(SCORE(index,str2num(PCs{2}(3))));
            zp(i_rat) = nanmean(SCORE(index,str2num(PCs{3}(3))));
            if animale<100
                gr(i_rat)=1;
            elseif animale>100 && animale<200
                gr(i_rat)=2;
            else
                gr(i_rat)=3;
            end
        end
        x_r = xp;
        y_r = yp;
        z_r = zp;
        scatter3(nanmean(x_r(gr==1)),nanmean(y_r(gr==1)),nanmean(z_r(gr==1)),200,'MarkerFaceColor', mycolor{pl_day},'MarkerEdgeColor','k','DisplayName',[dayPlot{pl_day},'-REHA'])
        hold on
        scatter3(x_r(gr==1),y_r(gr==1),z_r(gr==1),20,'MarkerFaceColor', mycolor{pl_day},'MarkerEdgeColor', mycolor{pl_day}, 'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.5,'HandleVisibility','off')
        scatter3(nanmean(x_r(gr==2)),nanmean(y_r(gr==2)),nanmean(z_r(gr==2)),200,'MarkerFaceColor', mycolor{pl_day+nt_d},'MarkerEdgeColor','k','DisplayName',[dayPlot{pl_day},'-ROB'])
        hold on
        scatter3(x_r(gr==2),y_r(gr==2),z_r(gr==2),20,'MarkerFaceColor', mycolor{pl_day+nt_d},'MarkerEdgeColor', mycolor{pl_day+nt_d}, 'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.5,'HandleVisibility','off')
        scatter3(nanmean(x_r(gr==3)),nanmean(y_r(gr==3)),nanmean(z_r(gr==3)),200,'MarkerFaceColor', mycolor{pl_day+2*nt_d},'MarkerEdgeColor','k','DisplayName',[dayPlot{pl_day},'-REHA+ROB'])
        hold on
        scatter3(x_r(gr==3),y_r(gr==3),z_r(gr==3),20,'MarkerFaceColor', mycolor{pl_day+2*nt_d},'MarkerEdgeColor', mycolor{pl_day+2*nt_d}, 'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.5,'HandleVisibility','off')
        zp_all {pl_day} = zp;
        yp_all {pl_day} = yp;
        xp_all {pl_day} = xp;
        gr_all {pl_day} = gr;
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
dataPlotx = NaN(nt_d,n_treat);
errhighx = NaN(nt_d,n_treat);
dataPloty = NaN(nt_d,n_treat);
errhighy = NaN(nt_d,n_treat);
dataPlotz = NaN(nt_d,n_treat);
errhighz = NaN(nt_d,n_treat);
dataPlotPC = NaN(nt_d,n_treat);
errhighPC = NaN(nt_d,n_treat);
pc_statd_s = cell(nt_d,1);
 for pl_day = 1:nt_d     
%      pc_statd = sqrt((yp_all{pl_day}-yp_all{1}).^2+(xp_all{pl_day}-xp_all{1}).^2+...
%      +(zp_all{pl_day}-zp_all{1}).^2);
     pc_statd = sqrt((yp_all{pl_day}).^2+(xp_all{pl_day}).^2+...
     +(zp_all{pl_day}).^2);
    pc_statd_s{pl_day} = pc_statd;
%      pc_statd_s{pl_day} = pc_statd;
%      yp_stat = [yp_stat;yp_all{pl_day}-yp_all{1}];
%      xp_stat = [xp_stat;xp_all{pl_day}-xp_all{1}];
%      zp_stat = [zp_stat;zp_all{pl_day}-zp_all{1}];
%      yp_stat = [yp_stat;yp_all{pl_day}];
%      xp_stat = [xp_stat;xp_all{pl_day}];
%      zp_stat = [zp_stat;zp_all{pl_day}];
%      dataPlotx(pl_day) = nanmean(xp_all{pl_day}-xp_all{1});
    for pl_gr = 1:n_treat
         dataPlotx(pl_day,pl_gr) = nanmean(xp_all{pl_day}(gr==pl_gr));
         errhighx(pl_day,pl_gr) = nanstd(xp_all{pl_day}(gr==pl_gr))/sqrt(length(yp_all{pl_day}(gr==pl_gr)));
         dataPloty(pl_day,pl_gr) = nanmean(yp_all{pl_day}(gr==pl_gr));
         errhighy(pl_day,pl_gr) = nanstd(yp_all{pl_day}(gr==pl_gr))/sqrt(length(yp_all{pl_day}(gr==pl_gr)));
         dataPlotz(pl_day,pl_gr) = nanmean(zp_all{pl_day}(gr==pl_gr));
         errhighz(pl_day,pl_gr) = nanstd(zp_all{pl_day}(gr==pl_gr))/sqrt(length(zp_all{pl_day}(gr==pl_gr)));
         dataPlotPC(pl_day,pl_gr) = nanmean(pc_statd(gr==pl_gr));
         errhighPC(pl_day,pl_gr) = nanstd(pc_statd(gr==pl_gr))/sqrt(length(pc_statd(gr==pl_gr)));         
    end
 end

subplot(2,2,2)
bb=bar(dataPlotx);
ylabel (PCs{1},'FontSize',13,'FontName','Arial')
hold on
% Find the number of groups and the number of bars in each group
ngroups = size(dataPlotx, 1);
nbars = size(dataPlotx, 2);
% Calculate the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
% Set the position of each error bar in the centre of the main bar
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
for i = 1:nbars
    % Calculate center of each bar
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, dataPlotx(:,i), errhighx(:,i), 'k', 'linestyle', 'none');
    x_point = ones(1,n_gr(i));
    for i_dd = 1:nt_d
        scatter(x_point*x(i_dd),(xp_all{i_dd}(gr==i)),20,'k','MarkerFaceColor','k','MarkerEdgeColor','none')
    end        
end
xlim([0 nt_d+1])

subplot(2,2,3)
bb=bar(dataPloty);
ylabel (PCs{2},'FontSize',13,'FontName','Arial')
hold on
% Find the number of groups and the number of bars in each group
ngroups = size(dataPloty, 1);
nbars = size(dataPloty, 2);
% Calculate the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
% Set the position of each error bar in the centre of the main bar
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
for i = 1:nbars
    % Calculate center of each bar
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, dataPloty(:,i), errhighy(:,i), 'k', 'linestyle', 'none');
    x_point = ones(1,n_gr(i));
    for i_dd = 1:nt_d
        scatter(x_point*x(i_dd),(yp_all{i_dd}(gr==i)),20,'k','MarkerFaceColor','k','MarkerEdgeColor','none')
    end        
end
xlim([0 nt_d+1])

subplot(2,2,4)
bb=bar(dataPlotz);
ylabel (PCs{3},'FontSize',13,'FontName','Arial')
hold on
% Find the number of groups and the number of bars in each group
ngroups = size(dataPlotz, 1);
nbars = size(dataPlotz, 2);
% Calculate the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
% Set the position of each error bar in the centre of the main bar
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
for i = 1:nbars
    % Calculate center of each bar
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, dataPlotz(:,i), errhighz(:,i), 'k', 'linestyle', 'none');
    x_point = ones(1,n_gr(i));
    for i_dd = 1:nt_d
        scatter(x_point*x(i_dd),(zp_all{i_dd}(gr==i)),20,'k','MarkerFaceColor','k','MarkerEdgeColor','none')
    end        
end
xlim([0 nt_d+1])

subplot(2,2,1)
bb=bar(dataPlotPC);
ylabel ('PCs','FontSize',13,'FontName','Arial')
hold on
% Find the number of groups and the number of bars in each group
ngroups = size(dataPlotPC, 1);
nbars = size(dataPlotPC, 2);
% Calculate the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
% Set the position of each error bar in the centre of the main bar
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
for i = 1:nbars
    % Calculate center of each bar
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, dataPlotPC(:,i), errhighPC(:,i), 'k', 'linestyle', 'none');
    x_point = ones(1,n_gr(i));
    for i_dd = 1:nt_d
        scatter(x_point*x(i_dd),(pc_statd_s{i_dd}(gr==i)),20,'k','MarkerFaceColor','k','MarkerEdgeColor','none')
    end        
end
xlim([0 nt_d+1])

%% Stats for distances 
% 2-ways Anova repeated measures

% stats = rm_anova2(Y,S,F1,F2,FACTNAMES);
% NON FUNZIONA perchè due gruppi di animali
%    Y          dependent variable (numeric) in a column vector
%    S          grouping variable for SUBJECT
%    F1         grouping variable for factor #1
%    F2         grouping variable for factor #2
%    FACTNAMES  a cell array w/ two char arrays: {'factor1', 'factor2'}

% x variable
FACTNAMES = {'Days','Group'};
Y = [];
S = [];
F1 = [];
F2 = [];
vect_an = [];
for f2 = 1:n_treat %factor 2 treatment
    vect_an = [vect_an; ones(n_gr(f2),1)*f2];        
end
for f1 = 1:nt_d %factor 1 days
    Y = [Y;xp_all{f1}];
    F2 = [F2;vect_an];
    S = [S;ListRat];
    F1 = [F1;ones(length(ListRat),1)*f1];
end
% results to copy in sigmaplot for repeated measures anova 2 ways

% kruskal wallis day by day
for f1 = 1:nt_d
   [p,~,stat] = kruskalwallis(xp_all{f1},vect_an);
end
% kruskal wallis for groups
for f2 = 1:n_treat
    x_stat = [];
    g_stat = [];
    for f1 = 1:nt_d
        en = sum(n_gr(1:f2));
        st = en-n_gr(f2)+1;
        x_stat = [x_stat;xp_all{f1}(st:en)];
        g_stat = [g_stat;ones(n_gr(f2),1)*f1];
    end
   [p,~,stat] = kruskalwallis(x_stat,g_stat);
   c = multcompare(stat);
end

end    
%% Fig more important parameters

pc = [COEFF(:,str2num(PCs{1}(3))),COEFF(:,str2num(PCs{2}(3))),COEFF(:,str2num(PCs{3}(3)))];
pc = abs(pc);
xbar = [1:nt_d];
for i_pc = 1:size(pc,2)
    figure(580+i_pc)
    [~,in] = max(pc(:,i_pc));
    subplot (1,2,1)
    pc(in,i_pc) = 0;
    hold on
    datapl = permute(dataPlot(:,in,:),[1,3,2]);
    err = permute (errhigh(:,in,:),[1,3,2]);
    bb=bar(datapl);
    ylabel (NameParam{in},'FontSize',13,'FontName','Arial')
    hold on
    % Find the number of groups and the number of bars in each group
    ngroups = size(datapl, 1);
    nbars = size(datapl, 2);
    % Calculate the width for each bar group
    groupwidth = min(0.8, nbars/(nbars + 1.5));
    % Set the position of each error bar in the centre of the main bar
    % Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
    for i = 1:nbars
        % Calculate center of each bar
        x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
        errorbar(x, datapl(:,i), err(:,i), 'k', 'linestyle', 'none');
        x_point = ones(1,n_gr(i));
        list = find(gr==i);
        for i_dd = 1:nt_d
            scatter(x_point*x(i_dd),MM(list,in,i_dd),20,'k','MarkerFaceColor','k','MarkerEdgeColor','none')
        end
    end
    hold off
    xlim([0 length(xbar)+1])
    xticks(xbar)
    xticklabels(dayPlot);
    
    [~,in] = max(pc(:,i_pc));
    subplot (1,2,2)
    pc(in,i_pc) =0;
    hold on
    datapl = permute(dataPlot(:,in,:),[1,3,2]);
    err = permute (errhigh(:,in,:),[1,3,2]);
    bb=bar(datapl);
    ylabel (NameParam{in},'FontSize',13,'FontName','Arial')
    hold on
    % Find the number of groups and the number of bars in each group
    ngroups = size(datapl, 1);
    nbars = size(datapl, 2);
    % Calculate the width for each bar group
    groupwidth = min(0.8, nbars/(nbars + 1.5));
    % Set the position of each error bar in the centre of the main bar
    % Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
    for i = 1:nbars
        % Calculate center of each bar
        x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
        errorbar(x, datapl(:,i), err(:,i), 'k', 'linestyle', 'none');
        x_point = ones(1,n_gr(i));
        list = find(gr==i);
        for i_dd = 1:nt_d
            scatter(x_point*x(i_dd),MM(list,in,i_dd),20,'k','MarkerFaceColor','k','MarkerEdgeColor','none')
        end
    end
    hold off
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