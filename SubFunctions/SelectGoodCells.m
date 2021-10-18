%2. Choose useful cell: delete cells that are not present enough in
%multiple recordings
function [indexMap] = SelectGoodCells(data_inscopix,Rat,add,DAYS,DAYSBL,DAYSBLtot)

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
% passive also in this case
%%%%%%

% delete columns that correspond to days that are not present in DAYSBL
for i_gg = size(DAYSBLtot,2):-1:1
    if ~contains(DAYSBL,DAYSBLtot{i_gg})
        indexMap(:,i_gg*2-1:i_gg*2) = [];
    end
end

% If a cell is not present in at least a quarter of the recordings we don't
% consider the unit significative

rec_tot = size(indexMap,2);%n tot recordings

for i_row = size(indexMap,1):-1:1
    ind_zeros = find(indexMap(i_row,:)==0);
    if length(ind_zeros)> (rec_tot-rec_tot/4)
        indexMap(i_row,:)=[];
    end
end

end