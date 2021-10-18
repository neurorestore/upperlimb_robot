%% Function to plot a pie chart for every day after injury to verify how cells
% that are active in baseline change their activity after injury

function [number_cell_type, number_evolution_BL, number_evolution_DPI] = PieChartActive_vDays(cell_type,DAYS)

col_DPI = contains(DAYS,'DPI');
ind_DPIc = find(col_DPI==1); % index DPI for name cells
n_days = sum(col_DPI)+1; %total number of days (all BL count as 1 day)

title_plot = {'BL',DAYS{ind_DPIc}};

labels = {'active movement','silent movement','not responding','no detected unit'};

t1col = [153 194 181]/255;   %red        
t2col = [179 223 127]/255;   %green
t3col = [253 192 145]/255;   %blue
t4col = [1 0 1];   %magenta
tilecolor = [t1col; t2col; t3col; t4col]; 

ind_BL = find(cell_type(:,1)==1); % index of all active units in baseline
ind_BL_in = find(cell_type(:,1)==2);
ind_BL_no = find(cell_type(:,1)==0);

figure(2100); hold on
subplot(2,n_days,1)
BL = [length(ind_BL),length(ind_BL_in),length(ind_BL_no)];
pie(BL,[1 0 0])
legend(labels(BL>0))
colormap(gca,tilecolor(BL>0,:))
title(title_plot{1})

ind_DPI_all = NaN(n_days-1,3);
number_evolution_BL = NaN(n_days-1,4);
number_evolution_DPI = NaN(n_days-1,4);

for i_gg = 2:size(cell_type,2)
    BLact_act = length(find(cell_type(ind_BL,i_gg)==1)); % number of units that are still active after injury
    BLact_nan = sum(isnan(cell_type(ind_BL,i_gg))); % number of units that are not found anymore after injury
    BLact_zero = length(find(cell_type(ind_BL,i_gg)==0)); % number of units that are not significant after injury
    BLact_in = length(find(cell_type(ind_BL,i_gg)==2)); % number of units that become silent after injury

    number_evolution_BL(i_gg-1,:) = [BLact_act,BLact_in,BLact_zero,BLact_nan];
   ind_DPI = find(cell_type(:,i_gg)==1); % index of all active units xDPI
   ind_DPI_all(i_gg-1,1) = length(ind_DPI);
   ind_DPI_in = find(cell_type(:,i_gg)==2);
   ind_DPI_no = find(cell_type(:,i_gg)==0);
   ind_DPI_all(i_gg-1,2) = length(ind_DPI_in);
   ind_DPI_all(i_gg-1,3) = length(ind_DPI_no);
   DPIact_act = length(find(cell_type(ind_DPI,1)==1)); % number of units that are still active after injury
   DPIact_nan = sum(isnan(cell_type(ind_DPI,1))); % number of units that are not found anymore after injury
   DPIact_zero = length(find(cell_type(ind_DPI,1)==0)); % number of units that are not significant after injury
   DPIact_in = length(find(cell_type(ind_DPI,1)==2)); % number of units that become silent after injury
   
   number_evolution_DPI(i_gg-1,:) = [DPIact_act,DPIact_in,DPIact_zero,DPIact_nan];
   
   figure(2100);hold on
   subplot(2,n_days,n_days+i_gg)
   DPI = [length(ind_DPI),length(ind_DPI_in),length(ind_DPI_no)];
   pie(DPI)
   colormap(gca,tilecolor(DPI>0,:))
   legend('off')
   title(title_plot{i_gg}(7:end))
   
   figure(2100);hold on
   subplot(2,n_days,i_gg)
   BLtoDPI = [BLact_act,BLact_in, BLact_zero,BLact_nan]/length(ind_BL);
   bb = bar(BLtoDPI);
   bb.EdgeColor = 'flat';
   bb.FaceColor ='flat';
   bb.CData(1,:) = tilecolor(1,:);
   bb.CData(2,:) = tilecolor(2,:);
   bb.CData(3,:) = tilecolor(3,:);
   bb.CData(4,:) = tilecolor(4,:);
    
%    figure; hold on
%    subplot(1,2,1)
%    BL = [BLact_act,BLact_in, BLact_zero,BLact_nan ]/length(ind_BL);
%    pie(BL)
%    legend(labels(BL>0))
%    colormap(gca,tilecolor(BL>0,:))
%    title(title_plot{1})
% 
%    subplot(1,2,2)
%    DPI = [DPIact_act,DPIact_nan, DPIact_zero, DPIact_in]/length(ind_DPI);
%    pie(DPI);
%    legend(labels(DPI>0))
%    colormap(gca,tilecolor(DPI>0,:))
%    title(title_plot{i_gg})

end

number_cell_type = [length(ind_BL),length(ind_BL_in),length(ind_BL_no);ind_DPI_all];

end