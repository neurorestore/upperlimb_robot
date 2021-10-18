% Extract Data from Baseline Matrix

function [Data_tot] = LoadDataBL(path,DAYS,Rat)

n_days = size(DAYS,2);
Data_tot = cell(5,n_days*2);

for i_gg = 1:n_days
    ListDAY = dir(path);
    ind = NaN;
    Half = [];
    Active = [];
    for i_d = 3:length(ListDAY) % find the correct folder
        if contains(ListDAY(i_d).name,DAYS{i_gg}(1:5))
            ind = i_d;
        end
    end
    ListFiles = dir ([path,ListDAY(ind).name,'\',Rat]);
    for i_f = 3:length(ListFiles) % find the correct file
        if contains(ListFiles(i_f).name,'Peak_x04_Ins.mat')
            load ([path,ListDAY(ind).name,'\',Rat,'\',ListFiles(i_f).name]);
            if contains(Data.info.Status,'Push_active')
                Half = Data;
            else
                Active = Data;
            end
            clear Data
        end
    end
    if ~isempty(Active) || ~isempty(Half)
        Data_tot{1,i_gg*2-1} = Active.Recorded_Data;
        Data_tot{2,i_gg*2-1} = Active.VICON;
        Data_tot{3,i_gg*2-1} = Active.SIMI;
        Data_tot{4,i_gg*2-1} = Active.INSCOPIX;
        Data_tot{5,i_gg*2-1} = Active.good_trials;
        Data_tot{1,i_gg*2} = Half.Recorded_Data;
        Data_tot{2,i_gg*2} = Half.VICON;
        Data_tot{3,i_gg*2} = Half.SIMI;
        Data_tot{4,i_gg*2} = Half.INSCOPIX;
        Data_tot{5,i_gg*2} = Half.good_trials;
    else
        disp(['Check ',DAYS{i_gg},', Data are missing'])
    end    
end
end