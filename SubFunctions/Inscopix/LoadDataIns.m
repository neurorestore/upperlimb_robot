% Extract Data from different days Matrix

function [Data_tot] = LoadDataIns(path,DAYS,Rat,Task)

n_days = size(DAYS,2);
Data_tot = cell(5,n_days); % matrix where to save all data

for i_gg = 1:n_days
    ListDAY = dir(path);
    ind = NaN;
    Data_task = [];
    for i_d = 3:length(ListDAY) % find the correct folder
        if contains(ListDAY(i_d).name,DAYS{i_gg}(1:5))
            ind = i_d;
        end
    end
    ListFiles = dir ([path,ListDAY(ind).name,'\',Rat]);
    for i_f = 3:length(ListFiles) % find the correct file
        if contains(ListFiles(i_f).name,'Peak_x04_Ins.mat')
            load ([path,ListDAY(ind).name,'\',Rat,'\',ListFiles(i_f).name]);
            if strcmp(Data.info.Status,Task)
                Data_task = Data;
            end
            clear Data
        end
    end
    if ~isempty(Data_task)
        Data_tot{1,i_gg} = Data_task.Recorded_Data;
        Data_tot{2,i_gg} = Data_task.VICON;
        Data_tot{3,i_gg} = Data_task.SIMI;
        Data_tot{4,i_gg} = Data_task.INSCOPIX;
        Data_tot{5,i_gg} = Data_task.good_trials;
    else
        disp(['Check ',DAYS{i_gg},', Data are missing'])
    end    
end
end