%% Compare the amplitude of the single cells after injury respect to baseline
% Choose only units that are present in all recordings 
% Considerate the real fluorescence and check the evolution

function AmplitudeSingleCells_vDays(Data,indexMap,MENO,LITTLEPIU,in_case,LITTLE,fSrobot,fScalcium,fSkin,DAYS,cell_type)


col_DPI = contains(DAYS,'DPI');
ind_DPI = find(col_DPI==1);
n_days = sum(col_DPI)+1;

ind_allcell = ~any(isnan(cell_type), 2)& cell_type(:,1)==1;
indexMap = indexMap(ind_allcell,:);

Area_allCells = cell(size(indexMap,1),n_days);

% cycle for all the remaing cells
for i_cell = 1:size(indexMap,1)
    Area_allFiles = [];
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
                    lineall = Data{4,i_file}.cells_signal(:,n_syn_cell);
                    linec = lineall(st:en);%+abs(min(lineall(st:en)));
                    M(i_ind,:) = linec;                    
                end
            end
        else
            M=NaN(1,round(MENO*fScalcium+LITTLEPIU*fScalcium+1));
        end
        %integration of the signal before and after the index
        Area = NaN(size(M,1),2);
        for i_row = 1:size(M,1)
            middle = MENO*fScalcium;
            fine = (MENO+LITTLE)*fScalcium+1;
            Area(i_row,1) = trapz(M(i_row,1:middle))/MENO;
            %Area(i_row,2) = trapz(M(i_row,middle+1:fine))/LITTLE;
            Area(i_row,2) = max(M(i_row,middle+1:fine));
        end
        Area(isoutlier(Area(:,2)),:)=[];
        if i_file < ind_DPI(1)
            Area_allFiles = [Area_allFiles;Area];
            Area_allCells{i_cell,1} = Area_allFiles;
        else
            Area_allCells{i_cell,i_file-ind_DPI(1)+2} = Area;
        end
    end
end

%% Calculate the mean value of the area for every cell and the relation after DPI and Baseline

%cycle for every units
data_plot = NaN(size(Area_allCells));
std_plot = NaN(size(Area_allCells));
for i_cell = 1:size(Area_allCells,1)
    data_stat = [];
    group_stat = [];
    for i_gg = 1:size(Area_allCells,2)
        data_plot(i_cell,i_gg) = nanmedian(Area_allCells{i_cell,i_gg}(:,2));
        std_plot(i_cell,i_gg) = nanstd(Area_allCells{i_cell,i_gg}(:,2));
        data_stat = [data_stat;Area_allCells{i_cell,i_gg}(:,2)];
        group_stat = [group_stat;ones(size(Area_allCells{i_cell,i_gg},1),1)*i_gg];
    end
    [p,~,stat] = kruskalwallis(data_stat,group_stat);
    c = multcompare(stat);
end

for i_row = 1:size(data_plot,1)
    figure
    bar(data_plot(i_row,:))
end

end
