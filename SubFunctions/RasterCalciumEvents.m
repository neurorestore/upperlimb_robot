%% Function to find the spikes of the single units and to plot them file by
% file while resting and activity phases are marked

function RasterCalciumEvents(indexMap,Data,cell_type,fScalcium,fSrobot,in_case)
INT = 1; % one seconds intervals for the activity to consider

ras_in_tot = [];
ras_re_tot = [];
mycolor = {[243 146 0]./255,[128 0 128]./255};
%% Plot raster for every file and create a matrix with the counted events
for i_file = 1:size(indexMap,2)
    figure(400+i_file); hold on
    i_cellFile = 1;
    locs_tot = [];
    pks_tot = [];
    info_cell = [];
    for i_cell = 1:length(cell_type)
        if cell_type(i_cell)==1 && indexMap(i_cell,i_file)>0.5 % condition excitatory cell and existing cell
            signal = zscore(Data{4,i_file}.cells_signal(:,indexMap(i_cell,i_file)));
            th_pr = 2*std(signal);
            [pks,locs,w,prom] = findpeaks(signal,fScalcium,'MinPeakProminence',th_pr,'Annotate','extents');
            scatter(locs,ones(length(locs),1)*i_cellFile,25,'k','filled')
            locs_tot = [locs_tot;locs];
            pks_tot = [pks_tot;pks];
            info_cell = [info_cell;ones(length(locs),1)*i_cell];
            i_cellFile = i_cellFile+1;
        end        
    end
    switch in_case
        case 1
            index = Data{1,i_file}.Analysis.Fzpeaks(:,1); %index of the force peak
            index = round(index/fSrobot*fScalcium);
        case 2
            index = Data{3,i_file}.start;
            index = round(index*fScalcium);
%                     pos_tot = Data{1,i_file}.pos_Spindle_drive.data;
%                     speed = derivative(pos_tot',1/fSrobot);
%                     acc = derivative(speed,1/fSrobot);
%                         th = mean(speed)+std(speed);
%                         [~,index] = findpeaks_GUI(-speed,'MINPEAKHEIGHT',th,'MINPEAKDISTANCE',10);
%                         onset = NaN(1,length(index));
%                         for n_p = 1:length(index)
%                             punti_neg = find(acc(1:index(n_p))<0 & speed(1:index(n_p))<th);
%                             if ~isempty(punti_neg)
%                                 onset(n_p) = punti_neg(end);
%                             end
%                         end
%             if mod(i_file,2)>0.5
%                 status = Data{1,i_file}.T_status.data;
%                 inin = find(status==2);
%                 index = inin(find (diff(inin)>1.5)+1);
%                 index = [inin(1),index];
%                 index = round(index/fSrobot*fScalcium);
%             else
%                 index = Data{3,i_file}.pks{2};
%                 index = round(index*fScalcium);
%             end
    end
    rest = Data{1,i_file}.Analysis.OnsetRest/fSrobot;
    index = index/fScalcium;
    raster_in = NaN(length(index),1);
    activity_in = NaN(length(index),2);
%     activity_in = [];
    for i_ac = 1:length(index)
        ind_ac = find(locs_tot>index(i_ac) & locs_tot<(index(i_ac)+INT));
        raster_in(i_ac) = length(ind_ac);
        activity_in(i_ac,1) = sum(pks_tot(ind_ac));
        activity_in(i_ac,2) = -Data{1,i_file}.Analysis.Fzpeaks(i_ac,4);
%         vv = Data{1,i_file}.Analysis.Fzpeaks(i_ac,4);
%         activity_in = [activity_in; pks_tot(ind_ac),ones(length(ind_ac),1)*vv,info_cell(ind_ac)];
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
    ras_in_tot = [ras_in_tot; raster_in,ones(length(raster_in),1)*i_file];
    ras_re_tot = [ras_re_tot; raster_rest,ones(length(raster_rest),1)*i_file];
    %p = kruskalwallis([raster_in(:,1); raster_rest(:,1)],[ones(length(raster_in),1);ones(length(raster_rest),1)*2]);
    
%     for iii_cell = 1:length(cell_type)
%     ind_cell = find(activity_in(:,3)==iii_cell);
%     if ~isempty(ind_cell)
        figure
        set(gca,'FontSize',14)
        hold on
        activity_in((activity_in(:,1)==0),:)= [];
        scatter(activity_in(:,2),activity_in(:,1),'ko','LineWidth',1)
        [R,P] = corrcoef(activity_in(:,1:2));
        [p,S] = polyfit(activity_in(:,2),activity_in(:,1),1);
        x = [0:0.01:max(activity_in(:,2))+0.05];
        [y_fit,delta] = polyval(p,x,S);
        plot(x,y_fit,'Color',[0 113 188]/255,'LineWidth',1.5)
        plot(x,y_fit+2*delta,'b--',x,y_fit-2*delta,'b--')
        xlim ([0,max(activity_in(:,2))+0.05])
        xlabel('Force Peaks')
        ylabel('Fluorescence active neurons')
        title (['CorrCoeff ',num2str(R(1,2)),' p=',num2str(P(1,2))])
%         title (['Cell # ',num2str(iii_cell)])
%     end
%     end
end

%% Evaluate activity in resting and peaks phase
figure; hold on
bar ([mean(ras_in_tot(:,1)) mean(ras_re_tot(:,1))])
errorbar([mean(ras_in_tot(:,1)) mean(ras_re_tot(:,1))],[std(ras_in_tot(:,1))/sqrt(length(ras_in_tot)) std(ras_re_tot(:,1))/sqrt(length(ras_re_tot))],'k.')
ylabel ('Spike counts (in 1 s)')
xticks([1,2])
xticklabels({'Action','Rest'})

p = kruskalwallis([ras_in_tot(:,1); ras_re_tot(:,1)],[ones(length(ras_in_tot),1);ones(length(ras_re_tot),1)*2])

%% Evaluate activity in peaks phase during active and half task
ras_in_half = [];
ras_re_half = [];
ras_in_act = [];
ras_re_act = [];
for i_file = 1:size(indexMap,2)
    if mod(i_file,2)>0.5
        index = find(ras_in_tot(:,2)==i_file);
        ras_in_act = [ras_in_act;ras_in_tot(index,1)];
        index = find(ras_re_tot(:,2)==i_file);
        ras_re_act = [ras_re_act;ras_re_tot(index,1)];
    else
        index = find(ras_in_tot(:,2)==i_file);
        ras_in_half = [ras_in_half;ras_in_tot(index,1)];
        index = find(ras_re_tot(:,2)==i_file);
        ras_re_half = [ras_re_half;ras_re_tot(index,1)];
    end
end
data = [ras_in_act;ras_re_act;ras_in_half;ras_re_half];
group = [ones(length(ras_in_act),1)*1;ones(length(ras_re_act),1)*2;ones(length(ras_in_half),1)*3;ones(length(ras_re_half),1)*4];
[p,~,stat] = kruskalwallis(data,group);
c = multcompare(stat);

figure; hold on
xbar =([1,2,4,5]);
bb = bar (xbar,[mean(ras_in_half(:,1)) mean(ras_re_half(:,1)) mean(ras_in_act(:,1)) mean(ras_re_act(:,1))])
errorbar(xbar,[mean(ras_in_half(:,1)) mean(ras_re_half(:,1)) mean(ras_in_act(:,1)) mean(ras_re_act(:,1))]...
    ,[std(ras_in_half(:,1))/sqrt(length(ras_in_half)) std(ras_re_half(:,1))/sqrt(length(ras_re_half)) ...
    std(ras_in_act(:,1))/sqrt(length(ras_in_act)) std(ras_re_act(:,1))/sqrt(length(ras_re_act))],'k.')
ylabel ('Spike counts (in 1 s)','FontSize',13,'FontName','Arial')
x = ones(length(ras_in_half),1);
scatter(x,ras_in_half,10,'k','filled');
x = ones(length(ras_re_half),1)*2;
scatter(x,ras_re_half,10,'k','filled');
x = ones(length(ras_in_act),1)*4;
scatter(x,ras_in_act,10,'k','filled');
x = ones(length(ras_re_act),1)*5;
scatter(x,ras_re_act,10,'k','filled');
xticks([1,2,4,5])
xticklabels({'Action','Rest','Action','Rest'})
bb.EdgeColor = 'flat';
bb.FaceColor ='flat';
bb.CData(1,:) = mycolor{1};
bb.CData(2,:)= mycolor{1};
bb.CData(3,:) = mycolor{2};
bb.CData(4,:)= mycolor{2};

end