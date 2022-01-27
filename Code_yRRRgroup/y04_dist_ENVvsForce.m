%% take all distances envelope-force peaks in half task
clear all
close all


Path = cd;
datapath = 'C:\R-Platform\DATA\2020_08_RRR_group\Robot Data';

ListDay = dir(datapath);
DataToShow = {'BL','2 Weeks','4 Weeks','8 Weeks'}; 
gr_names = {'REHA','ROB','REH+ROB'};
fS = 1000;
dist = NaN (18,2); % n ratti, n_emg
dist_ALL = cell(4,1);

mycolor = {[209 211 212]./255,[253,227,186]./255,[255,192,139]./255,[255,100,87]./255,... % Color reha %[255,100,87]./255,
    [136,141,145]./255,[185,202,233]./255,[60,89,164]./255,[40,65,120]./255,... % Color Rob %[40,65,120]./255,
    [63,66,67]./255,[206,232,215]./255,[107,187,136]./255,[57,137,86]./255}; % Color Rob+reha %,[57,137,86]./255

for i_d = 4:length(ListDay) % scroll days
        ListRat = dir([datapath,'\',ListDay(i_d,1).name]);
        
        for i_rat = 3:length(ListRat)-1 % scroll rats, the last file is the Read me
            
            filepath = [datapath,'\',ListDay(i_d,1).name,'\',ListRat(i_rat,1).name];
            ListFiles = dir(filepath);
            for i_file = 3: length(ListFiles)
                
                if contains(ListFiles(i_file,1).name,'AnalysisPeak_x04')
                    load([filepath,'\',ListFiles(i_file,1).name]);
                    if strcmp(Data.info.Status,'Push_active')
                        dist(i_rat-2,1) = nanmean(Data.VICON.EnvOnsetF.rel_max(:,2,1))/fS; %biceps
                        dist(i_rat-2,2) = nanmean(Data.VICON.EnvOnsetF.rel_max(:,2,2))/fS; %triceps
                    end                    
                end
            end
            
        end
        dist_ALL{i_d-3} = dist;
        dist = NaN (18,2);
end

%stat
group = [ones(1,18),ones(1,18)*2];
for i_dpi = 1:length(dist_ALL)
    x_stat = [dist_ALL{i_dpi}(:,1);dist_ALL{i_dpi}(:,2)];%n rats, n _emg
    [p,~,stat] = kruskalwallis(x_stat,group);
end

%% Figures all groups
figure
hold on
xbar = [1:2];
x_point = ones(1,18);
for i_day = 1:length(DataToShow)
    subplot(2,2,i_day)
    hold on
    dataPlot = [nanmean(dist_ALL{i_day}(:,1)),nanmean(dist_ALL{i_day}(:,2))];
    errhigh = [std(dist_ALL{i_day}(:,1)/sqrt(15)),nanstd(dist_ALL{i_day}(:,2)/sqrt(15))];
    bb=bar(xbar,dataPlot,'stacked');
    bb.EdgeColor = 'flat';
    bb.FaceColor ='flat';
    bb.CData= mycolor{i_day};
    hold on
    er = errorbar(xbar,dataPlot,[0,0],errhigh,'HandleVisibility','off');
    er.Color = [0 0 0];                            
    er.LineStyle = 'none';
    er.LineWidth = 1;
    scatter(x_point,(dist_ALL{i_day}(:,1)),20,'k','MarkerFaceColor','k','MarkerEdgeColor','none')
    scatter(x_point*2,(dist_ALL{i_day}(:,2)),20,'k','MarkerFaceColor','k','MarkerEdgeColor','none')
    xlim([0 length(xbar)+1])
    xticks(xbar)
end


%% Divide for groups (REHA, ROB, REH+ROB)
datapath_gr = 'C:\R-Platform\DATA\2020_08_RRR_group\Code_yRRRgroup\Groups\';
filenames_gr = {'Rehab_group.txt','Rob_group.txt','Rob+Rehab_group.txt'};
groups = zeros(18,1);
n_treat = length(filenames_gr);
for i_gr = 1:n_treat
    fpt = fopen([datapath_gr,filenames_gr{i_gr}]);  % Open programming students grade file
    namesROB = fscanf(fpt, '%s',[1 inf]);
    for i_an = 1:length(namesROB)/3
        numROB = str2num(namesROB((i_an-1)*3+2:(i_an-1)*3+3));
        groups(numROB)=i_gr;
    end
end
n_gr = [sum(groups ==1); sum(groups ==2);sum(groups ==3)]; 

%stat for groups
group = [ones(1,5),ones(1,5)*2];
for i_dpi = 1:length(dist_ALL)
    for i_gr = 1:n_treat
        rows = find(groups==i_gr);
        x_stat = [dist_ALL{i_dpi}(rows,1);dist_ALL{i_dpi}(rows,2)];%n rats, n _emg
        [p,~,stat] = kruskalwallis(x_stat,group);
    end
end

%figure for groups
cd([Path,'\Subfunctions'])
for i_gr = 1:n_treat
    figure
    hold on
    xbar = [1:2];
    x_point = ones(1,5);
    rows = find(groups==i_gr);
    for i_day = 1:length(DataToShow)
        subplot(2,2,i_day)
        hold on
        dataPlot = [nanmean(dist_ALL{i_day}(rows,1)),nanmean(dist_ALL{i_day}(rows,2))];
        errhigh = [std(dist_ALL{i_day}(rows,1)/sqrt(15)),nanstd(dist_ALL{i_day}(rows,2)/sqrt(15))];
        bb=bar(xbar,dataPlot,'stacked');
        bb.EdgeColor = 'flat';
        bb.FaceColor ='flat';
        bb.CData= mycolor{i_day};
        hold on
        er = errorbar(xbar,dataPlot,[0,0],errhigh,'HandleVisibility','off');
        er.Color = [0 0 0];                            
        er.LineStyle = 'none';
        er.LineWidth = 1;
        scatter(x_point,(dist_ALL{i_day}(rows,1)),20,'k','MarkerFaceColor','k','MarkerEdgeColor','none')
        scatter(x_point*2,(dist_ALL{i_day}(rows,2)),20,'k','MarkerFaceColor','k','MarkerEdgeColor','none')
        xlim([0 length(xbar)+1])
        xticks(xbar)
    end
    mtit(gr_names{i_gr})
end

