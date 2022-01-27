%% take all correlation coefficient Force-Bicep and compare them in the two tasks
clear all
close all


Path = cd;
datapath = 'C:\R-Platform\DATA\2020_08_RRR_group\Robot Data';
ListDay = dir(datapath);
%%
% Along Days

DataToShow = {'BL','2 Weeks','4 Weeks','8 Weeks'}; 
nt_d = length(DataToShow);

ListPullDPI = NaN(4,18);
ListPullDPI_p = NaN(4,18);

mycolor = {[209 211 212]./255,[253,227,186]./255,[255,192,139]./255,[255,100,87]./255,... % Color reha %[255,100,87]./255,
    [136,141,145]./255,[185,202,233]./255,[60,89,164]./255,[40,65,120]./255,... % Color Rob %[40,65,120]./255,
    [63,66,67]./255,[206,232,215]./255,[107,187,136]./255,[57,137,86]./255}; % Color Rob+reha %,[57,137,86]./255

 
for i_d = 4:length(ListDay) % scroll days
        ListRat = dir([datapath,'\',ListDay(i_d,1).name]);
        
        for i_rat = 3:length(ListRat)-1 % scroll rats, the last file is the Read me
            
            filepath = [datapath,'\',ListDay(i_d,1).name,'\',ListRat(i_rat,1).name];
            ListFiles = dir(filepath);
            for i_file = 3: length(ListFiles)
                
                if contains(ListFiles(i_file,1).name,'BicForceRelation')
                    uiopen([filepath,'\',ListFiles(i_file,1).name]);
                    fig = gca;
                    in = strfind(fig.Title.String,',');
                    p = str2num(fig.Title.String(10:in-1));
                    pp = str2num(fig.Title.String(in+6:end));
                    
                    if ~isnan(p)
                        ind_st = strfind(ListFiles(i_file,1).name,'Bic');
                        name = ListFiles(i_file,1).name(1:ind_st-1);
                        for i_f = 3:length(ListFiles) % find the type of the task
                            if contains(ListFiles(i_f,1).name,[name,'AnalysisPeaks.mat'])
                                load ([filepath,'\',ListFiles(i_f,1).name]);
                                task = Data.info.Status;
                            end
                        end
                                                   
                        if strcmp(task,'Push_active')
                            ListPullDPI(i_d-3,i_rat-2) = p;
                            ListPullDPI_p(i_d-3,i_rat-2) = pp;
                        end
                        
                    end                   
                    
                end
            end
            
        end
end
save ('ListPullDPI','ListPullDPI')        
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
%% Plot Correl Coeff Pass vs Act and stat

% remove outliers
for i_row = 1:size(ListPullDPI,1)
    ind_out = isoutlier(ListPullDPI(i_row,:));
    ListPullDPI(i_row,ind_out) = NaN;
    ListPullDPI_p(i_row,ind_out) = NaN;
end
% stat along DAYS
x_stat = [ListPullDPI(1,groups==1), ListPullDPI(2,groups==1), ListPullDPI(3,groups==1), ListPullDPI(4,groups==1)];
group = [ones(1,5),ones(1,5)*2,ones(1,5)*3,ones(1,5)*4];
[p,~,stat] = kruskalwallis(x_stat,group);
c = multcompare(stat);

x_stat = [ListPullDPI(1,groups==2), ListPullDPI(2,groups==2), ListPullDPI(3,groups==2), ListPullDPI(4,groups==2)];
group = [ones(1,5),ones(1,5)*2,ones(1,5)*3,ones(1,5)*4];
[p,~,stat] = kruskalwallis(x_stat,group);
c = multcompare(stat);

x_stat = [ListPullDPI(1,groups==3), ListPullDPI(2,groups==3), ListPullDPI(3,groups==3), ListPullDPI(4,groups==3)];
group = [ones(1,5),ones(1,5)*2,ones(1,5)*3,ones(1,5)*4];
[p,~,stat] = kruskalwallis(x_stat,group);
c = multcompare(stat);

%  BAR plot
figure
hold on
xbar = [1:1:nt_d];
dataPlot = [nanmean(ListPullDPI(:,groups==1),2),nanmean(ListPullDPI(:,groups==2),2),nanmean(ListPullDPI(:,groups==3),2)];
err = [nanstd(ListPullDPI(:,groups==1),1,2)/sqrt(n_gr(1)),nanstd(ListPullDPI(:,groups==2),1,2)/sqrt(n_gr(2)),nanstd(ListPullDPI(:,groups==3),1,2)/sqrt(n_gr(3))];
bb=bar(dataPlot);
ylabel ('Correlation Coefficient','FontSize',13,'FontName','Arial')
hold on
% Find the number of groups and the number of bars in each group
ngroups = size(dataPlot, 1);
nbars = size(dataPlot, 2);
% Calculate the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
% Set the position of each error bar in the centre of the main bar
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
for i = 1:nbars
    % Calculate center of each bar
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, dataPlot(:,i), err(:,i), 'k', 'linestyle', 'none');
    x_point = ones(1,n_gr(i));
    list = find(groups==i);
    for i_dd = 1:nt_d
        scatter(x_point*(x(i_dd)),ListPullDPI(i_dd,list),20,'k','MarkerFaceColor','k','MarkerEdgeColor','none')
    end
end
hold off
xlim([0 length(xbar)+1])
xticks(xbar)
xticklabels(DataToShow);

%Plot all together
figure
hold on
dataPlot = nanmean(ListPullDPI,2);
err = nanstd(ListPullDPI,1,2)/sqrt(15);
bar(xbar,dataPlot)
errorbar(xbar,dataPlot,-err,err,'k', 'linestyle', 'none')