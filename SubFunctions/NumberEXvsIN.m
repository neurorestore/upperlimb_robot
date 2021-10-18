%% Function to compare the activity (in term of excitatory or inhibitory) of the units in different days
% 1-> excitatory; 2-> inhibitory; 0-> not linked
function NumberEXvsIN(Data,indexMap,MENO,LITTLEPIU,in_case,LITTLE,fSrobot,fScalcium)

p = NaN(size(indexMap,1),size(indexMap,2)); % save the p value for every unit for every day/task (column 1 == day1 active, column 2 = day1 half, column 3 = day2 active ecc)
cell_type = NaN(size(indexMap,1),size(indexMap,2));
mycolor = {[243 146 0]./255,[128 0 128]./255};
% cell_type_tot = NaN(size(indexMap,1),size(indexMap,2),2);
% % cycle for different indexes
% for j = 1:2
%     if j==1
%         in_case = 1;
%     else
%         in_case = 2;
%     end
%cycle for every file
for i_file = 1:size(indexMap,2)
    for i_cell = 1:size(indexMap,1)
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
%                     if floor(i_file/2)==i_file/2
%                         index = Data{3,i_file}.pks{3};
%                         index = round(index*fScalcium);
%                     else
% %                         pos_tot = Data{1,i_file}.pos_Spindle_drive.data;
% %                         speed = derivative(pos_tot',1/fSrobot);
% %                         acc = derivative(speed,1/fSrobot);
% %                         th = mean(speed)+std(speed);
% %                         [~,index] = findpeaks_GUI(-speed,'MINPEAKHEIGHT',th,'MINPEAKDISTANCE',10);
% %                         onset = NaN(1,length(index));
% %                         for n_p = 1:length(index)
% %                             punti_neg = find(acc(1:index(n_p))<0 & speed(1:index(n_p))<th);
% %                             if ~isempty(punti_neg)
% %                                 onset(n_p) = punti_neg(end);
% %                             end
% %                         end
%                         status = Data{1,i_file}.T_status.data;
%                         inin = find(status==2);
%                         index = inin(find (diff(inin)>1.5)+1);
%                         index = [inin(1),index];
%                         index = round(index/fSrobot*fScalcium);
%                     end
            end
            %create matrix with signal of i_cell
            M = NaN(length(index),round(MENO*fScalcium+LITTLEPIU*fScalcium+1));
            n_syn_cell = indexMap(i_cell,i_file);
            %cycle for all the index
            for i_ind = 1:length(index)
                st = round(index(i_ind)-MENO*fScalcium);
                en = round(index(i_ind)+LITTLEPIU*fScalcium);
                if st>1 && en<size(Data{4,i_file}.cells_signal,1)
                    lineall = zscore(Data{4,i_file}.cells_signal(:,n_syn_cell));
                    M(i_ind,:) = lineall(st:en);                 
                end
            end
            %integration of the signal before and after the index
            Area = NaN(size(M,1),2);
            for i_row = 1:size(M,1)
                middle = MENO*fScalcium;
                fine = (MENO+LITTLE)*fScalcium+1;
                Area(i_row,1) = trapz(M(i_row,1:middle))/MENO;
                Area(i_row,2) = trapz(M(i_row,middle+1:fine))/LITTLE;
            end
            if sum(~isnan(Area(:,1)))>0.5 && sum(~isnan(Area(:,2)))>0.5 % check for units that have something different from nan value
                p(i_cell,i_file) = kruskalwallis([Area(:,1);Area(:,2)],[ones(i_row,1);ones(i_row,1)*2],'off');
                if p(i_cell,i_file) < 0.05 % statistical significant
                    if nanmean(Area(:,1))<nanmean(Area(:,2))
                        cell_type(i_cell,i_file) = 1; % excitatory
                    else
                        cell_type(i_cell,i_file) = 2; %inhibitory
                    end
                else
                    cell_type(i_cell,i_file) = 0; % not involved cell
                end
            else
                cell_type(i_cell,i_file) = 0;
            end
        else
            M=NaN(1,round(MENO*fScalcium+LITTLEPIU*fScalcium+1));
            index = 1;
        end       

    end
end
% cell_type_tot(:,:,j) = cell_type;
% end

% % % count number of cell that are active excitatory and active inhibitory
% number_group = NaN (2,size(indexMap,2),2);
% for i_gr = 1:size(indexMap,2)
%     for j = 1:2
%         number_group(1,i_gr,j) = length(find(cell_type_tot(:,i_gr,j)==1))/(size(cell_type_tot,1)-sum(isnan(cell_type_tot(:,i_gr,j))));
%         number_group(2,i_gr,j) = length(find(cell_type_tot(:,i_gr,j)==2))/(size(cell_type_tot,1)-sum(isnan(cell_type_tot(:,i_gr,j))));
%     end
% end
% number_type = NaN(2,2,2,size(indexMap,2)/2);
% for i_gr = 1:2:size(indexMap,2)
%     for j=1:2
%         number_type(1,:,j,(i_gr+1)/2) = number_group(1,i_gr:i_gr+1,j);
%         number_type(2,:,j,(i_gr+1)/2) = number_group(2,i_gr:i_gr+1,j);
%     end
% end
% 
% % stat on different recordings
% data = [number_group(1,1:2:end,1),number_group(1,1:2:end,2),number_group(1,2:2:end,1),number_group(1,2:2:end,2),...
%     number_group(2,1:2:end,1),number_group(2,1:2:end,2),number_group(2,2:2:end,1),number_group(2,2:2:end,2)];
% group_stat = [ones(size(indexMap,2)/2,1);ones(size(indexMap,2)/2,1)*2;ones(size(indexMap,2)/2,1)*3;ones(size(indexMap,2)/2,1)*4;...
%     ones(size(indexMap,2)/2,1)*5;ones(size(indexMap,2)/2,1)*6;ones(size(indexMap,2)/2,1)*7;ones(size(indexMap,2)/2,1)*8];
% [p,~,stat] = kruskalwallis(data,group_stat);
% multcompare(stat)
% 
% figure
% hold on
% media = mean(number_type,4);
% deviazione = std(number_type,0,4)/(size(indexMap,2)/2);
% bar([1,4,8,11],[media(1,:,1),media(2,:,1)])
% errorbar([1,4,8,11],[media(1,:,1),media(2,:,1)],[deviazione(1,:,1),deviazione(2,:,1)],'k.','LineWidth',1)
% bar([2,5,9,12],[media(1,:,2),media(2,:,2)])
% errorbar([2,5,9,12],[media(1,:,2),media(2,:,2)],[deviazione(1,:,2),deviazione(2,:,2)],'k.','LineWidth',1)
% xticks([1.5,4.5,8.5,11.5])
% xticklabels({'Active Ex','Half Ex','Active In','Half In'})
% ylabel ('Percentuage of cells involved in the task')

for i_cell = 1:size(cell_type,1)
    if max(cell_type(i_cell,:))>1 && isempty(find(cell_type(i_cell,:)==1))
       index = find(cell_type(i_cell,:)==0);
        cell_type (i_cell,index) = 2;
    end
end

% % count number of cell that are active excitatory and active inhibitory
number_group = NaN (3,size(indexMap,2));
for i_gr = 1:size(indexMap,2)
    number_group(1,i_gr) = length(find(cell_type(:,i_gr)==1))/(size(cell_type,1)-sum(isnan(cell_type(:,i_gr))));
    number_group(2,i_gr) = length(find(cell_type(:,i_gr)==2))/(size(cell_type,1)-sum(isnan(cell_type(:,i_gr))));
    number_group(3,i_gr) = length(find(cell_type(:,i_gr)==0))/(size(cell_type,1)-sum(isnan(cell_type(:,i_gr))));
end
number_type = NaN(3,2,size(indexMap,2)/2);
for i_gr = 1:2:size(indexMap,2)
    number_type(1,:,(i_gr+1)/2) = number_group(1,i_gr:i_gr+1)*100;
    number_type(2,:,(i_gr+1)/2) = number_group(2,i_gr:i_gr+1)*100;
    number_type(3,:,(i_gr+1)/2) = number_group(3,i_gr:i_gr+1)*100;
end

% stat on different recordings
data = [number_group(1,1:2:end),number_group(1,2:2:end),number_group(2,1:2:end),number_group(2,2:2:end),number_group(3,1:2:end),number_group(3,2:2:end)];
group_stat = [ones(size(indexMap,2)/2,1);ones(size(indexMap,2)/2,1)*2;ones(size(indexMap,2)/2,1)*3;ones(size(indexMap,2)/2,1)*4;ones(size(indexMap,2)/2,1)*5;ones(size(indexMap,2)/2,1)*6];
[p,~,stat] = kruskalwallis(data,group_stat);
c = multcompare(stat);

figure
hold on
media = mean(number_type,3);
deviazione = std(number_type,0,3)/sqrt((size(indexMap,2)/2));
bb = bar([1,2,4,5,7,8],[media(1,:),media(2,:),media(3,:)])
errorbar([1,2,4,5,7,8],[media(1,:),media(2,:),media(3,:)],[deviazione(1,:),deviazione(2,:),deviazione(3,:)],'k.','LineWidth',1)
x = ones(length(number_type(1,1,:)),1);
scatter(x,number_type(1,1,:),10,'k','filled');
x = ones(length(number_type(1,2,:)),1)*2;
scatter(x,number_type(1,2,:),10,'k','filled');
x = ones(length(number_type(2,1,:)),1)*4;
scatter(x,number_type(2,1,:),10,'k','filled');
x = ones(length(number_type(2,2,:)),1)*5;
scatter(x,number_type(2,2,:),10,'k','filled');
x = ones(length(number_type(3,1,:)),1)*7;
scatter(x,number_type(3,1,:),10,'k','filled');
x = ones(length(number_type(3,2,:)),1)*8;
scatter(x,number_type(3,2,:),10,'k','filled');
xticks([1.5,4.5,7.5])
xticklabels({'ActN','QuiN','IndN'})%'Half Ex','Active In','Half In'})
ylabel ('Percentuage of cells involved in the task','FontSize',13,'FontName','Arial')
bb.EdgeColor = 'flat';
bb.FaceColor ='flat';
bb.CData(1,:) = mycolor{2};
bb.CData(2,:)= mycolor{1};
bb.CData(3,:) = mycolor{2};
bb.CData(4,:)= mycolor{1};
bb.CData(5,:) = mycolor{2};
bb.CData(6,:)= mycolor{1};
end