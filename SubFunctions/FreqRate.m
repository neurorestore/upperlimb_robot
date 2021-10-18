% Function to find the spikes of the single units 

function  [p,f_rate] = FreqRate(indexMap,Data,cell_type,fScalcium,fSrobot)

% select same units for the files
f_rate = NaN(size(indexMap,1),2); %first column active and second column passive
for i_cell = 1:size(indexMap,1)
    indexCell = indexMap(i_cell,:);
    % check if at least one active and one passive have this unit
    indexCellAct = indexCell(1:2:end);
    indexCellHal = indexCell(2:2:end);
    if ~isempty(find(indexCellAct)~=0) && ~isempty(find(indexCellHal)~=0) && cell_type(i_cell)~=0
        clear indexCellAct indexCellHal
        
        spikes = cell(1,2);
        time = [0,0];
        for i_file = 1:2:size(Data,2) %scroll files
            % find common threshold
%             if indexCell(i_file)~=0
%                 act = Data{4,i_file}.cells_signal(:,indexCell(i_file));
%             else
%                 act = 0;
%             end
%             if indexCell(i_file+1)~=0
%                 half = Data{4,i_file+1}.cells_signal(:,indexCell(i_file+1));
%             else
%                 half = 0;
%             end
%             th = mean([act;half])+std([act;half]);
%             th_pr = std([act;half]);
            for i_tas = 1:2 % scroll task
                if indexCell(i_file+i_tas-1)~=0
                    
                    signal = zscore(Data{4,i_file+i_tas-1}.cells_signal(:,indexCell(i_file+i_tas-1)));
                    th_pr = 2*std(signal);
                    %[pks,locs] = findpeaks_GUI(signal,'MINPEAKHEIGHT',th,'MINPEAKDISTANCE',10);
                    [pks,locs,w,prom] = findpeaks(signal,fScalcium,'MinPeakProminence',th_pr,'Annotate','extents');
                    % check if it is a good way to find the activity of the
                    % units
%                     figure;hold on
%                     plot(signal)
%                     scatter(locs*fScalcium,pks)
%                     %xlim([0 1600])
%                     hold off
                    
%                     consider peaks in phase 2
%                     status = Data{1,i_file+i_tas-1}.T_status.data;
%                     index_pull = find(status ==2);
%                     for i_p = length(locs):-1:1
%                         if ~ismember(round(locs(i_p)*fSrobot),index_pull)
%                             locs(i_p)=[];
%                             pks(i_p) = [];
%                         end
%                     end
                    spikes{i_tas} = [spikes{i_tas};pks,locs];
                    time(i_tas) = time(i_tas)+length(signal)/fScalcium;
                    %time(i_tas) = time(i_tas)+length(index_pull)/fScalcium;
    
                end
            end  
        end
        
        % sum of spikes in different files and calculation of the frequency rate
        f_rate(i_cell,1) = size(spikes{1},1)/time(1);
        f_rate(i_cell,2) = size(spikes{2},1)/time(2);
    end
    figure(2001)
    hold on
    plot([1,2],f_rate(i_cell,:),'ko-','LineWidth',1.5)
end
xticks([1,2])
xticklabels({'Active','Half'})
ylabel('Frequency rate')
xlim([0,3])
i_row = size(f_rate,1);
p = kruskalwallis([f_rate(:,1);f_rate(:,2)],[ones(i_row,1);ones(i_row,1)*2]);%,'off');
end
    