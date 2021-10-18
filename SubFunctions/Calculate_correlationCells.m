% function to calculate the correlation coefficient between different cells
function R = Calculate_correlationCells(Half,Active,indexMap)

n_task = size(indexMap,2); % type of task (active, half)

%find only the cells that are common to the two tasks
for i_tas = 1: n_task
    indexMap(indexMap(:,i_tas)==0,:)=[];
end
% eliminate cells that have less than 5 activation in both tasks, we can
% imagine that they are not significative
% del_index = [];
% for i_c = 1:size(indexMap,1)
%     ev_ac = Active.INSCOPIX.events(:,indexMap(i_c,1));
%     ev_ac(isnan(ev_ac)==1)= [];
%     ev_ha = Half.INSCOPIX.events(:,indexMap(i_c,2));
%     ev_ha(isnan(ev_ha)==1)= [];
%     if length(ev_ac)<5 && length(ev_ha)<5
%         del_index = [del_index,i_c];
%     end
% end
% indexMap(del_index,:)=[];

% Find events cells
d = size(Half.INSCOPIX.cells_signal,2);
for i_c = 1:d
    signal = Half.INSCOPIX.cells_signal(:,i_c);
    %signal(index)=0;
    med = mean(signal);
    low = min(signal);
    dev = std(signal);
    th = low+2*dev;
    [pks,locs] = findpeaks_GUI(signal,'MINPEAKHEIGHT',th,'MINPEAKDISTANCE',20);
%     figure
%     hold on
%     plot(signal)
%     xlim([0 5000])
%     plot(Half.Recorded_Data.Fz.data(1:5:end)*500)
%     scatter(locs,pks)
%     title(num2str(dev))
    
%     figure
%     plot(Half.Recorded_Data.T_status.data)
%     hold on
%     stem(locs*4,ones(length(locs),1)*1.5)
%     ylim([0.5,4])
%     xlim([0 20000])
    if length(locs)>30
    matrix_F = NaN(length(locs),101);
    for i_cc = 1:length(locs)
        matrix_F(i_cc,:) = Half.Recorded_Data.Fz.data(locs(i_cc)*4-30:locs(i_cc)*4+70);
    end
    x = [-0.3:0.01:0.7];
    y = nanmean(matrix_F);
    y_dev = nanstd(matrix_F)/sqrt(length(locs));
    figure;hold on; set(gca, 'FontSize', 14)
    hold on
    plot(x,y,'LineWidth',2)
    X=[x, fliplr(x)];
    Y=[y + y_dev,fliplr(y - y_dev)];
    fill( X,Y,'b');
    alpha(.10)
    xlim([-0.3,0.7])
    ylim(get(gca, 'YLim'))
    line([0 0], get(gca, 'YLim'), 'Color', [0 0 0], 'LineStyle', '--')
    else
        disp('no peaks')
    end
end


R = NaN(size(indexMap,1)+1,size(indexMap,1)+1,n_task);
% calculate and plot correlation coefficient
for i_tas = 1:n_task
    if i_tas == 1
        signals = Active.INSCOPIX.cells_signal;
        Fz = Active.Recorded_Data.Fz.data(1:5:end);
        Fz = resample(Fz,length(signals),length(Fz));
    else
        signals = Half.INSCOPIX.cells_signal;
        Fz = Half.Recorded_Data.Fz.data(1:5:end);
        Fz = resample(Fz,length(signals),length(Fz));
    end
    matrix = NaN (size(signals,1),size(indexMap,1)+1); % in matrix common cells are aligned between different tasks
    for i_cell = 1:size(indexMap,1)
        matrix(:,i_cell) = signals(:,indexMap(i_cell,i_tas));
    end
    matrix(:,end)= Fz;
    R(:,:,i_tas) = corrcoef(matrix);
    figure
    imagesc(R(:,:,i_tas))
end

% estimate for each task the total correlation coefficient with the other
% cells and check if it is statistically different

sumCorr = NaN(size(R,2),n_task);
for i_tas = 1:n_task
    for i_cell = 1:size(R,2)
        sumCorr(i_cell,i_tas) = sum(R(:,i_cell,i_tas))-R(i_cell,i_cell,i_tas);        
    end
end
dim = size(sumCorr,1);
[p,tbl,stats] = kruskalwallis([sumCorr(:,1);sumCorr(:,2)],[ones(dim,1);ones(dim,1)*2]);