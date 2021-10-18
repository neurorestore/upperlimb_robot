% Code to plot
function PlotImagesc(Matrix_allCells,ind_ty,titolo,MENO,PIU,cell_type)

fScalcium = 20;
d = size(Matrix_allCells{1},2);

Matrix = NaN(length(ind_ty),d);
for i_cell = 1: length(ind_ty)
%     A = Matrix_allCells{ind_ty(i_cell)};
%     tutto = [];
%     for i = 1:size(A,1)
%         tutto = [tutto, A(i,:)];
%     end
% 
%     tutto = zscore(tutto);
%     B = NaN (size(A));
%     for i = 1:size(A,1)
%         B(i,:) = tutto(1+(i-1)*size(A,2):size(A,2)+(i-1)*size(A,2));
%     end
    Matrix (i_cell,:) = nanmean(Matrix_allCells{ind_ty(i_cell)});
end
x = [-MENO:1/fScalcium:PIU];
figure
imagesc(x,[1:1:size(Matrix,1)],Matrix);
colormap('jet')
colorbar
caxis([-1 1])
line([0 0], get(gca, 'YLim'), 'Color', [0 0 0], 'LineStyle', '--')
xlabel('Time (s)')
yyaxis left
ylabel('Neuron (#)')
yticks([1:length(ind_ty)])
yticklabels(num2str(ind_ty))
yyaxis right
yticks([0.5/length(ind_ty):1/length(ind_ty):1])
yticklabels(num2str(flip(cell_type(ind_ty))))
title(titolo)
end
