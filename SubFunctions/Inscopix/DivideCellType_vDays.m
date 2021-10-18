%% Function to find the type of every cell
% 1-> excitatory; 2-> inhibitory; 0-> not linked
function [Matrix_allCells,cell_type] = DivideCellType_vDays(Data,indexMap,MENO,LITTLEPIU,in_case,LITTLE,fSrobot,fScalcium,fSkin,DAYS,cell_type)

col_DPI = contains(DAYS,'DPI');
ind_DPI = find(col_DPI==1);
n_days = sum(col_DPI)+1;

Matrix_allCells = cell(size(indexMap,1),n_days);
%Info_days_allCells = cell(size(indexMap,1),1);
% Calculate the force peak
%[Mforce_allFiles,days_force] = CalculateMeanForcePeak(Data,in_case,MENO,LITTLEPIU,fSrobot,fScalcium,2);
% Calculate the movement
%[Mmov_allFiles,days_mov] = CalculateMeanMovement(Data,in_case,MENO,LITTLEPIU,fSrobot,fSkin,fScalcium,2);

% cycle for all the remaing cells
for i_cell = 1:size(indexMap,1)
    Matrix_allFiles = [];
    info_days = [];
    for i_file = 1 : size(indexMap,2) % cycle for every files
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
                    linec = lineall(st:en);
                    M(i_ind,:) = linec;                    
                end
            end
        else
            M=NaN(1,round(MENO*fScalcium+LITTLEPIU*fScalcium+1));
        end
        if i_file < ind_DPI(1)
            Matrix_allFiles = [Matrix_allFiles;M];
            Matrix_allCells{i_cell,1} = Matrix_allFiles;
        else
            Matrix_allCells{i_cell,i_file-ind_DPI(1)+2} = M;
        end
    end
end

% %% cycle for every day to evaluate the type of that unit
% 
% p = NaN(size(indexMap,1),n_days);
% cell_type = NaN(size(indexMap,1),n_days);
% 
% for i_cell = 1: size(indexMap,1)
%     for i_gg = 1: n_days
%         n_events = size(Matrix_allCells{i_cell,i_gg},1);
%         if n_events >1.5 % when it is only one row it means that the cell is not present in that file
%             Area = NaN(n_events,2);
%             for i_row = 1:n_events
%                 middle = MENO*fScalcium;
%                 fine = (MENO+LITTLE)*fScalcium+1;
%                 Area(i_row,1) = trapz(Matrix_allCells{i_cell,i_gg}(i_row,1:middle))/MENO;
%                 Area(i_row,2) = trapz(Matrix_allCells{i_cell,i_gg}(i_row,middle+1:fine))/LITTLE;
%             end
%             if sum(~isnan(Area(:,1)))>0.5 && sum(~isnan(Area(:,2)))>0.5
%                 p(i_cell,i_gg) = kruskalwallis([Area(:,1);Area(:,2)],[ones(i_row,1);ones(i_row,1)*2],'off');
%                 if p(i_cell,i_gg) < 0.05 % statistical significant
%                     if nanmean(Area(:,1))<nanmean(Area(:,2))
%                         cell_type(i_cell,i_gg) = 1; % excitatory
%                     else
%                         cell_type(i_cell,i_gg) = 2; %inhibitory
%                     end
%                 else
%                     cell_type(i_cell,i_gg) = 0; % not involved cell
%                 end
%             else
%                 cell_type(i_cell,i_gg) = NaN;
%             end
%             
%         end        
%     end
% end


% plot the cells activity around index for type
title_plot = {'BL',DAYS{ind_DPI}};
for i_gg = 1:n_days
    for i_type = 1:2
        ind_ty = find(cell_type(:,1)==i_type);
        %ind_ty = [1: size(indexMap,1)]';
        PlotImagesc(Matrix_allCells(:,i_gg),ind_ty,title_plot{i_gg},MENO,LITTLEPIU,cell_type(:,i_gg))
    end
end

end