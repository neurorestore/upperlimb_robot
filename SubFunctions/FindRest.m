% Find resting index for intervals of one seconds (no force peaks in the
% main direction
function Data = FindRest(Data,fSrobot)

time = 1; %seconds without force activity
INT = time*fSrobot;

for i_file = 1:size(Data,2)
    peaks = Data{1,i_file}.Analysis.Fzpeaks(:,1);
    force = Data{1,i_file}.Fz.data;
    onset = Data{1,i_file}.Analysis.Fzpeaks(:,2);
    endset = Data{1,i_file}.Analysis.Fzpeaks(:,3);
    diff = onset(2:end)-endset(1:end-1);
    index = find(diff>INT);
    rest_in = [];
    for j = 1:length(index)
        n_int = floor(diff(index(j))/INT);
        for i_in = 1:n_int
            st = endset(index(j))+1+(i_in-1)*INT;
            en = st+INT-1;
            % check if in an interval of one second there is not the area
            % of a previous peak or a consecutive
            distance = find(abs((peaks-st))<INT);
            % check if the standard deviation is too high or not
            if std(force(st:en))<0.06 && isempty(distance)
                rest_in = [rest_in,st];
            end
        end
    end
    % check if rest_in are in good trials
    for i_r = length(rest_in):-1:1
        cicle = Data{1,i_file}.cicles.data(rest_in(i_r));
        if ~ismember(cicle, Data{5,i_file})
            rest_in(i_r)= [];
        end
    end
    Data{1,i_file}.Analysis.OnsetRest = rest_in;
end
end