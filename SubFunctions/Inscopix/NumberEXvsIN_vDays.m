%% Function to compare the activity (in term of excitatory or inhibitory) of the units in different days
% 1-> excitatory; 2-> inhibitory; 0-> not linked
function cell_type = NumberEXvsIN_vDays(Data,indexMap,MENO,LITTLEPIU,in_case,LITTLE,fSrobot,fScalcium,DAYS)

p = NaN(size(indexMap,1),size(indexMap,2)); % save the p value for every unit for every day for the selected task
cell_type = NaN(size(indexMap,1),size(indexMap,2));

%cycle for every file
for i_file = 1:size(indexMap,2)
    for i_cell = 1:size(indexMap,1)
        if indexMap(i_cell,i_file) ~= 0
            %find index to align neural activity
            switch in_case
                case 1 % index = force peaks
                   index = Data{1,i_file}.Analysis.Fzpeaks(:,1); %index of the force peak
                   phase = Data{1,i_file}.Analysis.Fzpeaks(:,15);
                   index(phase<1.5)=[]; % select only index during pulling phase
                                        % because we are working in the
                                        % half tasks
                   index = round(index/fSrobot*fScalcium);
                case 2 % index = start movement
                    index = Data{3,i_file}.start;
                    index = round(index*fScalcium);
            end
            %create matrix with signal of i_cell
            M = NaN(length(index),round(MENO*fScalcium+LITTLEPIU*fScalcium+1));
            n_syn_cell = indexMap(i_cell,i_file);
            %cycle for all the index
            for i_ind = 1:length(index)
                st = round(index(i_ind)-MENO*fScalcium);
                en = round(index(i_ind)+LITTLEPIU*fScalcium);
                if st>1 && en<size(Data{4,i_file}.cells_signal,1)
                    lineall = zscore(Data{4,i_file}.cells_signal(:,n_syn_cell));
                    M(i_ind,:) = lineall(st:en);                 
                end
            end
            %integration of the signal before and after the index
            Area = NaN(size(M,1),2);
            for i_row = 1:size(M,1)
                middle = MENO*fScalcium;
                fine = (MENO+LITTLE)*fScalcium+1;
                Area(i_row,1) = trapz(M(i_row,1:middle))/MENO;
                Area(i_row,2) = trapz(M(i_row,middle+1:fine))/LITTLE;
            end
            if sum(~isnan(Area(:,1)))>0.5 && sum(~isnan(Area(:,2)))>0.5 % check for units that have something different from nan value
                p(i_cell,i_file) = kruskalwallis([Area(:,1);Area(:,2)],[ones(i_row,1);ones(i_row,1)*2],'off');
                if p(i_cell,i_file) < 0.05 % statistical significant
                    if nanmean(Area(:,1))<nanmean(Area(:,2))
                        cell_type(i_cell,i_file) = 1; % excitatory
                    else
                        cell_type(i_cell,i_file) = 2; %inhibitory
                    end
                else
                    cell_type(i_cell,i_file) = 0; % not involved cell
                end
            else
                cell_type(i_cell,i_file) = 0;
            end

        else
            M=NaN(1,round(MENO*fScalcium+LITTLEPIU*fScalcium+1));
            index = 1;
        end        
    end
end

% % count number of cell that are active excitatory and active inhibitory
% find number days of baseline
col_DPI = contains(DAYS,'DPI');
number_group = NaN (2,sum(col_DPI)+1);
ind_DPI = find(col_DPI==1);
cell_typeBL = NaN(size(cell_type,1),1);
for i_gr = 1:sum(col_DPI)+1
    if i_gr == 1 % baseline recording
        n_bl = size(indexMap,2)-sum(col_DPI);
        tot = size(cell_type,1);
        for i_c = 1:size(cell_type,1)
            if sum(~isnan(cell_type(i_c,1:n_bl)))==0
                tot = tot-1;
            else
                cell_typeBL(i_c)= max(cell_type(i_c,1:n_bl));
            end
        end
        number_group(1,1) = length(find(cell_typeBL==1))/tot;
        number_group(2,1) = length(find(cell_typeBL==2))/tot;
    else % days DPI
        tot = size(cell_type,1)-sum(isnan(cell_type(:,ind_DPI(i_gr-1)))); % total number of units for that day
        number_group(1,i_gr) = length(find(cell_type(:,ind_DPI(i_gr-1))==1))/tot;
        number_group(2,i_gr) = length(find(cell_type(:,ind_DPI(i_gr-1))==2))/tot;
    end
end

cell_type = [cell_typeBL,cell_type(:,n_bl+1:end)];
figure
set(gca,'FontSize',14)
hold on
bb = bar(number_group'*100)
xticks([1:sum(col_DPI)+1])
xticklabels({'BL',DAYS{ind_DPI}})
ylabel ('Percentuage of cells involved in the task')

end