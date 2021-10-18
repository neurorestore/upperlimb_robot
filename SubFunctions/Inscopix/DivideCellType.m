%% Function to find the type of every cell
% 1-> excitatory; 2-> inhibitory; 0-> not linked
function [Matrix_allCells,Matrix_allCellsAct,cell_type,cell_type_act] = DivideCellType(Data,indexMap,MENO,LITTLEPIU,in_case,LITTLE,fSrobot,fScalcium,fSkin)

p = NaN(size(indexMap,1),1);
pact = NaN(size(indexMap,1),1);
%p1 = NaN(size(indexMap,1),1);
cell_type = NaN(size(indexMap,1),1);
cell_type_act = NaN(size(indexMap,1),1);
Matrix_allCells = cell(size(indexMap,1),1);
Info_days_allCells = cell(size(indexMap,1),1);
Distance_allCells = cell(size(indexMap,1),1);
% Calculate the force peak
[Mforce_allFiles,days_force] = CalculateMeanForcePeak(Data,in_case,MENO,LITTLEPIU,fSrobot,fScalcium,2);
% Calculate the movement
%[Mmov_allFiles,days_mov] = CalculateMeanMovement(Data,in_case,MENO,LITTLEPIU,fSrobot,fSkin,fScalcium,2);
% cycle for all the remaing cells
for i_cell = 1:size(indexMap,1)
    Matrix_allFiles = [];
    info_days = [];
    Distance_allFiles = [];
    for i_file = 2:2:size(indexMap,2)
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
%                    index = Data{3,i_file}.pks{3};
%                    index = round(index*fScalcium);
                
            end
            %create matrix with signal of i_cell
            M = NaN(length(index),round(MENO*fScalcium+LITTLEPIU*fScalcium+1));
            dist = NaN(length(index),1);
            n_syn_cell = indexMap(i_cell,i_file);
            %cycle for all the index
            for i_ind = 1:length(index)
                st = round(index(i_ind)-MENO*fScalcium);
                en = round(index(i_ind)+LITTLEPIU*fScalcium);
                if st>1 && en<size(Data{4,i_file}.cells_signal,1)
                    lineall = zscore(Data{4,i_file}.cells_signal(:,n_syn_cell));
                    linec = lineall(st:en);
                    %linec = Data{4,i_file}.cells_signal(st:en,n_syn_cell);
                    % bring the activity to zero level to better calculate
                    % the integrate of the force
                    mm = min(linec);
                    %linec = linec - mm;
                    M(i_ind,:) = linec;
                    [~,dist(i_ind)] = max(linec);
                    dist(i_ind) = dist(i_ind)/fScalcium-MENO;
                end
            end
        else
            M=NaN(1,round(MENO*fScalcium+LITTLEPIU*fScalcium+1));
            index = 1;
            dist = NaN;
        end
        Matrix_allFiles = [Matrix_allFiles;M];
        line_day = ones(length(index),1)*i_file;
        info_days = [info_days;line_day];
        Distance_allFiles = [Distance_allFiles;dist];
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
            cell_type(i_cell) = 1; % excitatory
        else
            cell_type(i_cell) = 2; %inhibitory
        end
    else
        cell_type(i_cell) = 0; % not involved cell
    end
    else
        cell_type(i_cell) = 0;
    end
    Matrix_allCells{i_cell} = Matrix_allFiles;
    Info_days_allCells{i_cell} = info_days;
    Distance_allCells{i_cell} = Distance_allFiles;
end
% plot the cells activity around index for type
title_plot = {'Half task-Excitatory','Half task-Inhibitory'};
for i_type = 1:2
    ind_ty = find(cell_type==i_type);
    PlotImagesc(Matrix_allCells,ind_ty,title_plot{i_type},MENO,LITTLEPIU,cell_type)
end

% Graph force vs calcium
%PlotGraphFC (Matrix_allCells,cell_type,Info_days_allCells,Mforce_allFiles,days_force,MENO,LITTLE,fScalcium,fSrobot)

%% Calculate if is there is significant differences in the time of
% activation
dist_stat = [];
dist_group = [];
for i_cell = 1:size(indexMap,1)
    if cell_type(i_cell)==1
        dist_stat = [dist_stat;Distance_allCells{i_cell}];
        dist_group = [dist_group;ones(length(Distance_allCells{i_cell}),1)*i_cell];
    end
end

[p,~,stat] = kruskalwallis(dist_stat,dist_group);
c = multcompare(stat);

% NO REGULAR PATTERN OF ACTIVATION! (NOT useful)
% % study the pattern of activity

% index_cells_active = find(cell_type==1);
% for i_rem = length(index_cells_active):-1:1
%     if length(Distance_allCells{index_cells_active(i_rem)})~=114
%         index_cells_active(i_rem)=[];
%     end
% end
% 
% figure; hold on
% for i_pks = 100:114
%     dist_pks = NaN(length(index_cells_active),1);
%     for i_cell = 1:length(index_cells_active)
%         dist_pks(i_cell) = Distance_allCells{index_cells_active(i_cell)}(i_pks);
%     end
%     [B,I] = sort(dist_pks);
%     plot(I)
% end

%% Calculate in the active file
% Calculate the force peak
[Mforce_allFiles,days_force] = CalculateMeanForcePeak(Data,in_case,MENO,LITTLEPIU,fSrobot,fScalcium,1);
% Calculate the movement
[Mmov_allFiles,days_mov] = CalculateMeanMovement(Data,in_case,MENO,LITTLEPIU,fSrobot,fSkin,fScalcium,1);
% cycle for all the remaing cells
Matrix_allCellsAct = cell(size(indexMap,1),1);
Info_days_allCells = cell(size(indexMap,1),1);
Distance_allCells = cell(size(indexMap,1),1);
for i_cell = 1:size(indexMap,1)
    Matrix_allFiles = [];
    info_days = [];
    Distance_allFiles = [];
    for i_file = 1:2:size(indexMap,2)
        if indexMap(i_cell,i_file) ~= 0
            %find index to align neural activity
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
%                         status = Data{1,i_file}.T_status.data;
%                         inin = find(status==2);
%                         index = inin(find (diff(inin)>1.5)+1);
%                         index = [inin(1),index];
%                         index = round(index/fSrobot*fScalcium);
%                    index = Data{3,i_file}.pks{2};
%                    index = round(index*fScalcium);
                    
            end
            %create matrix with signal of i_cell
            M = NaN(length(index),round(MENO*fScalcium+LITTLEPIU*fScalcium+1));
            n_syn_cell = indexMap(i_cell,i_file);
            dist = NaN(length(index),1);
            %cycle for all the index
            for i_ind = 1:length(index)
                st = round(index(i_ind)-MENO*fScalcium);
                en = round(index(i_ind)+LITTLEPIU*fScalcium);
                if st>1 && en<size(Data{4,i_file}.cells_signal,1)
                    lineall = zscore(Data{4,i_file}.cells_signal(:,n_syn_cell));
                    linec = lineall(st:en);
                    % bring the activity to zero level to better calculate
                    % the integrate of the force
                    mm = min(linec);
                    %linec = linec - mm;
                    M(i_ind,:) = linec;
                    [~,dist(i_ind)] = max(linec);
                    dist(i_ind) = dist(i_ind)/fScalcium-MENO;
                end                
            end
        else
            M=[];
            index = [];
        end
        Matrix_allFiles = [Matrix_allFiles;M];
        line_day = ones(length(index),1)*i_file;
        info_days = [info_days;line_day];
        Distance_allFiles = [Distance_allFiles;dist];
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
    pact(i_cell) = kruskalwallis([Area(:,1);Area(:,2)],[ones(i_row,1);ones(i_row,1)*2],'off');
    if pact(i_cell) < 0.01 % statistical significant
        if nanmean(Area(:,1))<nanmean(Area(:,2))
            cell_type_act(i_cell) = 1; % excitatory
        else
            cell_type_act(i_cell) = 2; %inhibitory
        end
    else
        cell_type_act(i_cell) = 0; % not involved cell
    end
    else
        cell_type_act(i_cell) = 0;
    end
    
    Matrix_allCellsAct{i_cell} = Matrix_allFiles;
    Info_days_allCells{i_cell} = info_days;
    Distance_allCells{i_cell} = Distance_allFiles;
end

% plot the cells activity around index for type
title_plot = {'Active task-Excitatory','Active task-Inhibitory'};
for i_type = 1:2
    ind_ty = find(cell_type==i_type);
    PlotImagesc(Matrix_allCellsAct,ind_ty,title_plot{i_type},MENO,LITTLEPIU,cell_type_act)
end

%% Calculate if is there is significant differences in the time of
% activation

% NO SIGNIFICANT RESULTS (NOT useful)

% dist_stat = [];
% dist_group = [];
% for i_cell = 1:size(indexMap,1)
%     if cell_type(i_cell)==1
%         dist_stat = [dist_stat;Distance_allCells{i_cell}];
%         dist_group = [dist_group;ones(length(Distance_allCells{i_cell}),1)*i_cell];
%     end
% end
% 
% [p,~,stat] = kruskalwallis(dist_stat,dist_group);
% c = multcompare(stat);

end