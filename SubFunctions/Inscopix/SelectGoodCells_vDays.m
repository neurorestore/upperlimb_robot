%2. Choose useful cell: delete cells that are not present enough in
%multiple recordings
function [indexMap] = SelectGoodCells_vDays(data_inscopix,Rat,add,DAYS,DAYSsel,DAYStot,Task)

%load the file output of .mat of CellReg.m
ListFiles = dir ([data_inscopix,Rat,add,'\',DAYS{1}]);
    
for i_f = 3: length(ListFiles)
    if contains(ListFiles(i_f).name,'cellRegistered')
        ind = i_f;
    end
end
    
if isempty(ind)
    disp('Error, the cellRegistered file does not exist')
    return
else
    load([data_inscopix,Rat,add,'\',DAYS{1},'\',ListFiles(ind).name]);
    % the mapping of each registered cell to the indices in each registered session
    indexMap = cell_registered_struct.cell_to_index_map;
    clear cell_registered_struct ListFiles ind i_f
end
% REMIND %
% in indexMap are saved indices of same units alongs days and session, in
% this case all data from baseline are put together(first active, then
% passive)
%%%%%%

% delete columns that correspond to days that are not present in DAYStot
% and that are not the selected task
for i_gg = size(DAYStot,2):-1:1
    if ~contains(DAYSsel,DAYStot{i_gg})
        indexMap(:,i_gg*2-1:i_gg*2) = [];
    end
end

if strcmp(Task,'Push_active')
    indexMap(:,1:2:end) = [];
else
    indexMap(:,2:2:end) = [];
end

% If a cell is not present in half of the recordings after injury at least,
% we don't consider it 

col_DPI = contains(DAYSsel,'DPI');

for i_row = size(indexMap,1):-1:1
    ind_zeros = find(indexMap(i_row,:)==0);
    if length(ind_zeros)> size(indexMap)/4 %(sum(col_DPI)/2)
        indexMap(i_row,:)=[];
    end
end

end
