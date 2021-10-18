% Fill NaN value in the matrix to perform the PCA
% info: column 1 - number of day
%       column 2 - number of animals
function matrix_out = Reconstruct_matrix_PCA (matrix_in,info,REM)

days = unique(info(:,1));
matrix_out = NaN(size(matrix_in));
n_par = size(matrix_in,2);
n_cy = size(matrix_in,1);

if REM % remove outliers
    for i_par = 1:n_par
        line_par = matrix_in(:,i_par);
        TF = isoutlier(line_par,'quartiles');
        %sum(TF)
        matrix_in(TF==1,i_par)=NaN;
    end
end



for i_par = 1:n_par
    tot_nan = sum(isnan(matrix_in(:,i_par)));
    line_par = matrix_in(:,i_par);
    if tot_nan> n_cy/2
        disp(['Delay parameter number ',num2str(i_par)]);
    else
        for i_d = 1:length(days)
            index = find(info(:,1)==days(i_d));
            in_nan = find(isnan(line_par(index))==1);
            if ~isempty(in_nan)
                %out = isoutlier(line_par(index));
                %index(out==1)= [];
                M = max(line_par(index));
                m = min(line_par(index));
                add = randi([0 100],length(in_nan),1);
                add = add*(M-m)/100+m;
                line_par(index(in_nan))= add;
            end
        end
    end 
    sum(isnan(line_par))
    matrix_out(:,i_par) = line_par;
end

end