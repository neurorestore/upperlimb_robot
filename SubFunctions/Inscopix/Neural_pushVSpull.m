%% Analysis to check activity of single neurons when there is a push instead that a pull movement

function Neural_pushVSpull(Data,indexMap,MENO,LITTLEPIU,in_case,LITTLE,fSrobot,fScalcium)


p = NaN(size(indexMap,1),1);
cell_type = NaN(size(indexMap,1),2);
Matrix_allCells = cell(size(indexMap,1),1);
Matrix_contra_allCells = cell(size(indexMap,1),1);

% Calculate the force peak
[~,~] = CalculateMeanForcePeak(Data,in_case,MENO,LITTLEPIU,fSrobot,fScalcium,2);

% cycle for all the remaing cells
for i_cell = 1:size(indexMap,1)
    Matrix_allFiles = [];
    Matrix_contra_allFiles = [];
    for i_file = 2:2:size(indexMap,2)
        if indexMap(i_cell,i_file) ~= 0
            %find index to align neural activity
            index = Data{1,i_file}.Analysis.Fzpeaks(:,1); %index of the force peak
            phase = Data{1,i_file}.Analysis.Fzpeaks(:,15);
            index(phase<1.5)=[]; % select only index during pulling phase
                                        % because we are working in the
                                        % half tasks
            index = round(index/fSrobot*fScalcium);
            
            % find index pks contra
            index_contra = Data{1,i_file}.Analysis.Fzpeaks_contra(:,1);
            %index_contra = Data{1,i_file}.Analysis.Fxpeaks_contra(:,1);
            index_contra = round(index_contra/fSrobot*fScalcium);

            %create matrix with signal of i_cell for pks
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
            
            %create matrix with signal of i_cell for pks contra
            M_contra = NaN(length(index_contra),round(MENO*fScalcium+LITTLEPIU*fScalcium+1));
            n_syn_cell = indexMap(i_cell,i_file);
            %cycle for all the index
            for i_ind = 1:length(index_contra)
                st = round(index_contra(i_ind)-MENO*fScalcium);
                en = round(index_contra(i_ind)+LITTLEPIU*fScalcium);
                if st>1 && en<size(Data{4,i_file}.cells_signal,1)
                    lineall = zscore(Data{4,i_file}.cells_signal(:,n_syn_cell));
                    linec = lineall(st:en);
                    M_contra(i_ind,:) = linec;
                end
            end
        else
            M=NaN(1,round(MENO*fScalcium+LITTLEPIU*fScalcium+1));
            M_contra=NaN(1,round(MENO*fScalcium+LITTLEPIU*fScalcium+1));
        end
        Matrix_allFiles = [Matrix_allFiles;M];
        Matrix_contra_allFiles = [Matrix_contra_allFiles;M_contra];
    end
    %integration of the signal before and after the index
    Area = NaN(size(Matrix_allFiles,1),2);
    for i_row = 1:size(Matrix_allFiles,1)
        middle = MENO*fScalcium;
        fine = (MENO+LITTLE)*fScalcium+1;
        Area(i_row,1) = trapz(Matrix_allFiles(i_row,1:middle))/MENO;
        Area(i_row,2) = trapz(Matrix_allFiles(i_row,middle+1:fine))/LITTLE;
    end
    if sum(~isnan(Area(:,1)))>0.5 && sum(~isnan(Area(:,2)))>0.5
    p(i_cell) = kruskalwallis([Area(:,1);Area(:,2)],[ones(i_row,1);ones(i_row,1)*2],'off');
    %p1(i_cell) = signrank(Area(:,2)./Area(:,1),1);
    if p(i_cell) < 0.01 % statistical significant
        if nanmean(Area(:,1))<nanmean(Area(:,2))
            cell_type(i_cell,1) = 1; % excitatory
        else
            cell_type(i_cell,1) = 2; %inhibitory
        end
    else
        cell_type(i_cell,1) = 0; % not involved cell
    end
    else
        cell_type(i_cell,1) = 0;
    end
    
    %integration of the signal before and after the index_contra
    Area = NaN(size(Matrix_contra_allFiles,1),2);
    for i_row = 1:size(Matrix_contra_allFiles,1)
        middle = MENO*fScalcium;
        fine = (MENO+LITTLE)*fScalcium+1;
        Area(i_row,1) = trapz(Matrix_contra_allFiles(i_row,1:middle))/MENO;
        Area(i_row,2) = trapz(Matrix_contra_allFiles(i_row,middle+1:fine))/LITTLE;
    end
    if sum(~isnan(Area(:,1)))>0.5 && sum(~isnan(Area(:,2)))>0.5
    p(i_cell) = kruskalwallis([Area(:,1);Area(:,2)],[ones(i_row,1);ones(i_row,1)*2],'off');
    if p(i_cell) < 0.01 % statistical significant
        if nanmean(Area(:,1))<nanmean(Area(:,2))
            cell_type(i_cell,2) = 1; % excitatory
        else
            cell_type(i_cell,2) = 2; %inhibitory
        end
    else
        cell_type(i_cell,2) = 0; % not involved cell
    end
    else
        cell_type(i_cell,2) = 0;
    end
    
    
    Matrix_allCells{i_cell} = Matrix_allFiles;
    Matrix_contra_allCells{i_cell} = Matrix_contra_allFiles;
end


% plot the cells activity around index for type
title_plot = {'Half task-Excitatory','Half task-Inhibitory'};
for i_type = 1:2
    ind_ty = find(cell_type(:,1)==i_type);
    PlotImagesc(Matrix_allCells,ind_ty,title_plot{i_type},MENO,LITTLEPIU,cell_type(:,1))
end

% plot the cells activity around index for type
title_plot = {'Half task-Excitatory','Half task-Inhibitory'};
for i_type = 1:2
    ind_ty = find(cell_type(:,2)==i_type);
    PlotImagesc(Matrix_contra_allCells,ind_ty,title_plot{i_type},MENO,LITTLEPIU,cell_type(:,2))
end

end