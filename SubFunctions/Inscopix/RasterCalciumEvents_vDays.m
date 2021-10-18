%% Function to find the spikes of the single units and to plot them file by
% file while resting and activity phases are marked

function [activity_in_ToT, RP_tot] = RasterCalciumEvents_vDays(indexMap,Data,cell_type,fScalcium,fSrobot,fSkin,in_case,DAYS)
INT = 1; % one seconds intervals for the activity to consider

col_DPI = contains(DAYS,'DPI');
ind_DPI = find(col_DPI==1);
n_days = sum(col_DPI)+1;
n_BL = length(find(col_DPI==0));

%reconstruct cell_type for every file
if n_BL>1
    for i_type = 2: n_BL
        cell_type = [cell_type(:,1),cell_type];
    end
end

% initialization of a variable where to save all the calcium events file by
% file
ras_in_tot = cell(length(DAYS),1); 
ras_re_tot = cell(length(DAYS),1);
activity_in_ToT = cell(1,size(indexMap,2));
RP_tot = NaN(size(indexMap,2),2);

%% Plot raster for every file and create a matrix with the counted events
for i_file = 1:size(indexMap,2)
    figure(400+i_file); hold on
    i_cellFile = 0;
    locs_tot = [];
    pks_tot = [];
    info_cell = [];
    for i_cell = 1:size(cell_type,1)
        if indexMap(i_cell,i_file)>0.5 % condition existing cell
            i_cellFile = i_cellFile+1;
            signal = zscore(Data{4,i_file}.cells_signal(:,indexMap(i_cell,i_file)));
            th_pr = 2*std(signal);
            [pks,locs,w,prom] = findpeaks(signal,fScalcium,'MinPeakProminence',th_pr,'Annotate','extents');
            scatter(locs,ones(length(locs),1)*i_cellFile,25,'k','filled')
            locs_tot = [locs_tot;locs];
            pks_tot = [pks_tot;pks];
            info_cell = [info_cell;ones(length(locs),1)*i_cell];            
        end        
    end
    switch in_case
        case 1
            index = Data{1,i_file}.Analysis.Fzpeaks(:,1); %index of the force peak
            index = round(index/fSrobot*fScalcium);
            peaks_corr = -Data{1,i_file}.Analysis.Fzpeaks(:,4);
        case 2
            index = Data{3,i_file}.start;
            peaks_corr = -Data{3,1}.pks{1}(Data{3,1}.pks{6}==1);
            index = round(index*fScalcium);            
    end
    rest = Data{1,i_file}.Analysis.OnsetRest/fSrobot;
    index = index/fScalcium;
    raster_in = NaN(length(index),1);
    activity_in = NaN(length(index),2);
    %matrix_activation = zeros(i_cellFile,length(index));
    for i_ac = 1:length(index)
        ind_ac = find(locs_tot>index(i_ac) & locs_tot<(index(i_ac)+INT));
        raster_in(i_ac) = length(ind_ac);
        activity_in(i_ac,1) = sum(pks_tot(ind_ac));
        activity_in(i_ac,2) = peaks_corr(i_ac);
        %matrix_activation(info_cell(ind_ac),i_ac)=1;
        xx = [index(i_ac) index(i_ac)+INT index(i_ac)+INT index(i_ac)];
        yy = [0 0 i_cellFile i_cellFile];
        p = patch( xx,yy,'r','EdgeColor','none');
        set(p,'FaceAlpha',0.3)
    end
    raster_rest = NaN(length(rest),1);
    for i_ac = 1:length(rest)
        raster_rest(i_ac) = length(find(locs_tot>rest(i_ac) & locs_tot<(rest(i_ac)+INT)));
        xx=[rest(i_ac) rest(i_ac)+INT rest(i_ac)+INT rest(i_ac)];
        yy=[0 0 i_cellFile i_cellFile];
        p = patch( xx,yy,'b','EdgeColor','none');
        set(p,'FaceAlpha',0.3)
    end
    ylim([0 i_cellFile])
    ras_in_tot{i_file} = raster_in;
    ras_re_tot{i_file} = raster_rest;  
    activity_in_ToT{i_file} = activity_in;
    
    figure
    set(gca,'FontSize',14)
    hold on
    scatter(activity_in(:,2),activity_in(:,1),'ko','LineWidth',1)
    [R,P] = corrcoef(activity_in(:,1:2));
    RP_tot(i_file,:) = [R(1,2),P(1,2)];
    [p,S] = polyfit(activity_in(:,2),activity_in(:,1),1);
    x = [0:0.01:max(activity_in(:,2))+0.05];
    [y_fit,delta] = polyval(p,x,S);
    plot(x,y_fit,'Color',[0 113 188]/255,'LineWidth',1.5)
    plot(x,y_fit+2*delta,'b--',x,y_fit-2*delta,'b--')
    xlim ([0,max(activity_in(:,2))+0.05])
    xlabel('Force Peaks')
    ylabel('Fluorescence active neurons')
    title (['CorrCoeff ',num2str(R(1,2)),' p=',num2str(P(1,2))])
end

%% Evaluate activity in resting and peaks phase in different days

title_plot = {'BL',DAYS{ind_DPI}};

ras_in_pl = [];
ras_re_pl = [];
data_in = [];
data_re = [];
group_in = [];
group_re = [];

for i_file = 1:size(indexMap,2)
    if i_file < n_BL
        ras_in_pl = [ras_in_pl;ras_in_tot{i_file}];
        ras_re_pl = [ras_re_pl;ras_re_tot{i_file}];
    else
        if i_file == n_BL
            ras_in_pl = [ras_in_pl;ras_in_tot{i_file}];
            ras_re_pl = [ras_re_pl;ras_re_tot{i_file}];
            tot_cell = mean(sum(indexMap(:,1:n_BL)~=0));
        else
            ras_in_pl = ras_in_tot{i_file};
            ras_re_pl = ras_re_tot{i_file};
            tot_cell = size(cell_type,1)-sum(isnan(cell_type(:,i_file)));
        end        
        figure; hold on
        ras_in_pl = ras_in_pl/tot_cell*100;
        ras_re_pl = ras_re_pl/tot_cell*100;
        bar ([mean(ras_in_pl) mean(ras_re_pl)])
        errorbar([mean(ras_in_pl) mean(ras_re_pl)],[std(ras_in_pl)/sqrt(length(ras_in_pl)) std(ras_re_pl)/sqrt(length(ras_re_pl))],'k.')
        ylabel ('Spike counts (%)')
        xticks([1,2])
        xticklabels({'Action','Rest'})
        title(title_plot{i_file-n_BL+1})
        data_in = [data_in;ras_in_pl];
        data_re = [data_re;ras_re_pl];
        group_in = [group_in;ones(length(ras_in_pl),1)*i_file];
        group_re = [group_re;ones(length(ras_re_pl),1)*(i_file+100)];
        p = kruskalwallis([ras_in_pl; ras_re_pl],[ones(length(ras_in_pl),1);ones(length(ras_re_pl),1)*2])        
    end
end

%[p,~,stat] = kruskalwallis([data_in;data_re],[group_in;group_re])
%c=multcompare(stat);

% figure all raster together
dd_in = unique(group_in);
dd_re = unique(group_re);
figure; hold on
xbar = [1:3:(length(dd_in)-1)*3+1];
for i_d = 1: length(dd_in)
    single_in = data_in(group_in==dd_in(i_d));
    single_re = data_re(group_re==dd_re(i_d));
    bb = bar([xbar(i_d) xbar(i_d)+1],[mean(single_in),mean(single_re)])
    errorbar([xbar(i_d) xbar(i_d)+1],[mean(single_in),mean(single_re)],[std(single_in)/sqrt(length(single_in)),std(single_re)/sqrt(length(single_re))],'k.') 
end
xticks(xbar+0.5)
xticklabels({'BL',DAYS{col_DPI}})
%bb.EdgeColor = 'flat';
%bb.FaceColor ='flat';
ylabel ('Spike counts (%)')
% bb.CData(1,:) = mycolor{1};
% bb.CData(2,:)= mycolor{2};
% bb.CData(3,:) = mycolor{3};
% bb.CData(4,:)= mycolor{4};
% bb.CData(5,:) = mycolor{5};
% bb.CData(6,:)= mycolor{6};
[p,~,stat] = kruskalwallis(data_in,group_in)
c=multcompare(stat);

end




