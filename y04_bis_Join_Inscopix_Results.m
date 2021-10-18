%% Join inscopix results

clear all
close all


Rat = {'R04','R08'};
data_inscopix = 'C:\R-Platform\DATA\2019_02_Group1\SUs-software inscopix\';
add = 'RegCells'; %folder with the matrix output of GUI CellReg
DAYSsel = {'BL','DAY09_4DPI','DAY10_8DPI','DAY11_14DPI','DAY12_21DPI','DAY13_35DPI'};
%Weeks = {'Week 1','Week 2','Week 3','Week 5'};
Weeks = {'4DPI','8DPI','14DPI','21DPI','35DPI'};
n_BL = 4; % number of days of baseline
n_day_tot = n_BL +length(DAYSsel)-1; % number of total days
n_day = length(Weeks)+1;

% Sample frequency
fSrobot = 100; %robot
fSemg = 1000; %VICON
fSkin = 50; %SIMI
fScalcium = 20; %inscopix

CodePath = cd;

%% Load data all inscopix animals
n_an = size(Rat,2);
allData = cell(n_an,1);
for i_rat = 1:n_an
    allData{i_rat} = load([data_inscopix,Rat{i_rat},add,'\Inscopix_results.mat']);
end

%% Analysis Force peaks vs Calcium peaks

R_BL = NaN(n_an*n_BL,1);
P_BL = NaN(n_an*n_BL,1);
R_DPI = NaN(n_an,length(Weeks));
P_DPI = NaN(n_an,length(Weeks));

for i_rat = 1:n_an
    
    for i_day = 1:n_BL
        RP = allData{i_rat}.Ins_results.RP_tot(i_day,:);
        R_BL(i_day+(i_rat-1)*n_BL) = RP(1);
        P_BL(i_day+(i_rat-1)*n_BL) = RP(2);
    end
    for i_day = n_BL+1:n_day_tot
        RP = allData{i_rat}.Ins_results.RP_tot(i_day,:);
        R_DPI(i_rat,i_day-n_BL) = RP(1);
        P_DPI(i_rat,i_day-n_BL) = RP(2);
    end
    
end
% DPI version
%R_mean = [mean(R_BL),mean(R_DPI)]; % mean correlation coefficient for multiple animals
% week version
% R_mean = [mean(R_BL),mean([R_DPI(:,1);R_DPI(:,2)]), mean(R_DPI(:,3:end))]; % mean correlation coefficient for multiple animals
% R_std = [std(R_BL)/sqrt(length(R_BL)),std([R_DPI(:,1);R_DPI(:,2)])/sqrt(4), std(R_DPI(:,3:end))/sqrt(2)];

R_mean = [mean(R_BL), mean(R_DPI(:,1:end))]; % mean correlation coefficient for multiple animals
R_std = [std(R_BL)/sqrt(length(R_BL)),std(R_DPI(:,1:end))/sqrt(2)];

%P_mean = [mean(P_BL),mean([P_DPI(:,1);P_DPI(:,2)]), mean(P_DPI(:,3:end))]; % mean correlation coefficient for multiple animals
%P_std = [std(P_BL)/sqrt(length(P_BL)),std([P_DPI(:,1);P_DPI(:,2)])/sqrt(4), std(P_DPI(:,3:end))/sqrt(2)];

P_mean = [mean(P_BL), mean([P_DPI(1,:),P_DPI(2,:)])];
P_std = [std(P_BL)/sqrt(length(P_BL)), std([P_DPI(1,:),P_DPI(2,:)])/sqrt(10)];
% stat on p value
p_stat = [P_BL;P_DPI(1,:)';P_DPI(2,:)'];
group_stat = [ones(1, length(P_BL)),ones(1,2*size(P_DPI,2))*2];
p = kruskalwallis(p_stat,group_stat);
%% Plot Correlation Coefficient mean
mycolor = {[243 146 0]./255,[0,51,102]./255,[0,128,128]./255,[0,168,0]./255,[130,167,0]./255,[255 237 0]./255}; %

%xbar = [1.1,1.9;4.1,4.9;7.1,7.9;10.1,10.9;13.1,13.9];
%xbar2 = [1.5,4.5,7.5,10.5,13.5];
xbar = [1:n_day];
figure
hold on
%bb=bar(xbar2,[R_mean;P_mean]');
bb=bar(xbar,R_mean);
for i_c = 1:1
    bb(i_c).EdgeColor = 'flat';
    bb(i_c).FaceColor ='flat';
    bb(i_c).CData(1,:) = mycolor{1};
    bb(i_c).CData(2,:) = mycolor{2};
    bb(i_c).CData(3,:)= mycolor{3};
    bb(i_c).CData(4,:)= mycolor{4};
    bb(i_c).CData(5,:)= mycolor{5};
    bb(i_c).CData(6,:)= mycolor{6};
end
        
ylabel ('Correlation Coefficient','FontSize',13,'FontName','Arial')
hold on
%er = errorbar(xbar,[R_mean;P_mean]',[-R_std;-P_std]',[+R_std;+P_std]','HandleVisibility','off');
er = errorbar(xbar,R_mean,-R_std,R_std,'HandleVisibility','off');
for i_c = 1:1
    er(i_c).Color = [0 0 0];                            
    er(i_c).LineStyle = 'none';
    er(i_c).LineWidth = 1;
end
xlim([0 length(xbar)+1])
xticks(xbar)
xticklabels(['BL',Weeks]);

figure
hold on
ylabel ('P value','FontSize',13,'FontName','Arial')
b = bar([1,2],P_mean);
b.FaceColor = 'flat';
b.CData = [100 100 100]/255;
er = errorbar([1,2],P_mean,-P_std,P_std,'HandleVisibility','off');
er.Color = [0 0 0];
er.LineStyle = 'none';
er.LineWidth = 1;
xlim([0 3])
xticks([1,2])
xticklabels({'BL','Injured'});


%% Plot mean graph of the evolution of units

labels = {'active movement','silent movement','not responding','no detected unit'};

t1col = [153 194 181]/255;          
t2col = [179 223 127]/255;   
t3col = [253 192 145]/255;
% t1col = [191 208 197]/255;          
% t2col = [225 234 151]/255;   
% t3col = [254 230 170]/255;
piecolor = [t1col; t2col; t3col];

barcolor = [[103 165 153]/255; [140 205 96]/255; [252 167 115]/255; [100 100 100]/255;[255 255 255]/255];

n_ct_all = NaN([size(allData{1}.Ins_results.n_cell_type),2]);
n_eBL_all = NaN([size(allData{1}.Ins_results.n_evolution_BL),2]);
n_eDPI_all = NaN([size(allData{1}.Ins_results.n_evolution_DPI),2]);

% join data different animals
for i = 1:n_an
    n_ct_all(:,:,i) = allData{i}.Ins_results.n_cell_type./sum(allData{i}.Ins_results.n_cell_type,2)*100;
    n_eBL_all(:,:,i) = allData{i}.Ins_results.n_evolution_BL./sum(allData{i}.Ins_results.n_evolution_BL,2)*100;
    n_eDPI_all(:,:,i) = allData{i}.Ins_results.n_evolution_DPI./sum(allData{i}.Ins_results.n_evolution_DPI,2)*100;
end
n_eBL_all (1,:,:) = [];
n_eDPI_all (1,:,:) = [];
% Plot pie chart
figure(2101);hold on
for i_day = 1:n_day
    if i_day ==1 %baseline
        subplot(2,n_day,i_day)
        min_pie = flip(min(n_ct_all(1,:,:),[],3));
        max_pie = flip(max(n_ct_all(1,:,:),[],3));
        names = {[num2str(round(min_pie(1))),'%-',num2str(round(max_pie(1))),'%'],[num2str(round(min_pie(2))),'%-',num2str(round(max_pie(2))),'%'],[num2str(round(min_pie(3))),'%-',num2str(round(max_pie(3))),'%']};
        pie(flip(mean(n_ct_all(1,:,:),3)),[0 0 1],names)
        colormap(gca,flip(piecolor))
        legend(flip(labels(1:3)))
        title('Baseline')
%     elseif i_day ==2
%         subplot(2,n_day,i_day+n_day)
%         min_pie = min(min(n_ct_all(2:3,:,:),[],3),[],1);
%         max_pie = max(max(n_ct_all(2:3,:,:),[],3),[],1);
%         names = {[num2str(round(min_pie(1))),'%-',num2str(round(max_pie(1))),'%'],[num2str(round(min_pie(2))),'%-',num2str(round(max_pie(2))),'%'],[num2str(round(min_pie(3))),'%-',num2str(round(max_pie(3))),'%']};
%         pie(mean(mean(n_ct_all(2:3,:,:),3),1),[0 0 0],names)
%         colormap(gca,piecolor)
%         legend('off')
%         title(Weeks{i_day-1})
    else
        subplot(2,n_day,i_day+n_day)
        min_pie = flip(min(n_ct_all(i_day,:,:),[],3));
        max_pie = flip(max(n_ct_all(i_day,:,:),[],3));
        names = {[num2str(round(min_pie(1))),'%-',num2str(round(max_pie(1))),'%'],[num2str(round(min_pie(2))),'%-',num2str(round(max_pie(2))),'%'],[num2str(round(min_pie(3))),'%-',num2str(round(max_pie(3))),'%']};
        pie(flip(mean(n_ct_all(i_day,:,:),3)),[0 0 1],names)
        colormap(gca,flip(piecolor))
        legend('off')
        title(Weeks{i_day-1})
    end
end

% Plot bar how is units change category after injury
% for i_day = 2:n_day
%     if i_day ==2
%         subplot(2,n_day,i_day)
%         bb = bar (mean(mean(n_eBL_all(1:2,:,:),3),1),'Stacked','FaceColor','flat');        
%     else
%         subplot(2,n_day,i_day)
%         bb = bar (mean(n_eBL_all(i_day,:,:),3));
%     end
%     bb.EdgeColor = 'flat';
%     bb.FaceColor ='flat';
%     bb.CData(1,:) = barcolor(1,:);
%     bb.CData(2,:) = barcolor(2,:);
%     bb.CData(3,:) = barcolor(3,:);
%     bb.CData(4,:) = barcolor(4,:);
%     ylim([0,60])
%     
% end

%% Plot Pie chart external
figure
n_altre = mean(n_ct_all(1,:,:),3);
n_altre_tot = n_altre(2)+n_altre(3); % # not responding and quiescente neurons
for i_day = 1:n_day-1
    subplot(1,n_day,i_day)
    pie(flip([mean(n_eBL_all(i_day,:,:),3)*n_altre(1)/100,n_altre_tot]))
    colormap(gca,flip(barcolor))
%     bb = bar (mean(n_eBL_all(i_day,:,:),3));
%     bb.EdgeColor = 'flat';
%     bb.FaceColor ='flat';
%     bb.CData(1,:) = barcolor(1,:);
%     bb.CData(2,:) = barcolor(2,:);
%     bb.CData(3,:) = barcolor(3,:);
%     bb.CData(4,:) = barcolor(4,:);
%     ylim([0,60])
    
end

figure
for i_day = 1:n_day-1
    n_altre = mean(n_ct_all(i_day+1,:,:),3);
    n_altre_tot = n_altre(2)+n_altre(3); % # not responding and quiescente neurons
    subplot(1,n_day,i_day)
    pie(flip([mean(n_eDPI_all(i_day,:,:),3)*n_altre(1)/100,n_altre_tot]))
    colormap(gca,flip(barcolor))
%     bb = bar (mean(n_eBL_all(i_day,:,:),3));
%     bb.EdgeColor = 'flat';
%     bb.FaceColor ='flat';
%     bb.CData(1,:) = barcolor(1,:);
%     bb.CData(2,:) = barcolor(2,:);
%     bb.CData(3,:) = barcolor(3,:);
%     bb.CData(4,:) = barcolor(4,:);
%     ylim([0,60])
    
end

%% Plot pie graph week5 and what they were in baseline
% big difference in the two animals
figure(2200)
subplot(1,2,1)
min_pie = flip(min(n_ct_all(end,:,:),[],3));
max_pie = flip(max(n_ct_all(end,:,:),[],3));
names = {[num2str(round(min_pie(1))),'%-',num2str(round(max_pie(1))),'%'],[num2str(round(min_pie(2))),'%-',num2str(round(max_pie(2))),'%'],[num2str(round(min_pie(3))),'%-',num2str(round(max_pie(3))),'%']};
pie(flip(mean(n_ct_all(end,:,:),3)),[0 0 1],names)
colormap(gca,flip(piecolor))
legend(flip(labels(1:3)))
title(Weeks{end})
subplot(1,2,2)
bb = bar (mean(n_eDPI_all(end,:,:),3));
bb.EdgeColor = 'flat';
bb.FaceColor ='flat';
bb.CData(1,:) = barcolor(1,:);
bb.CData(2,:) = barcolor(2,:);
bb.CData(3,:) = barcolor(3,:);
bb.CData(4,:) = barcolor(4,:);
