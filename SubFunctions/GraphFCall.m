% Function to compare the energy of the calcium signal compare to the force
% recorded by the load cell

function GraphFCall(Data, Matrix, indexMap, INT, fSrobot, fScalcium, fSemg, MENO, cell_type,in_case,act)

INTf = INT*fSrobot;
INTc = INT*fScalcium;
INTe = INT*fSemg;

n_int_rem = floor(1.5/INT)+1; % it is necessary to don't consider the last part of the recordings
                                  % because the interval for the calcium is
                                  % delayed compare to the one for the
                                  % force, so at the final part of the
                                  % recordings I don't have information
                                  % about the neural activity
X = [1/fScalcium:1/fScalcium:INT];
Xf = [1/fSrobot:1/fSrobot:INT];
Xe = [1/fSemg:1/fSemg:INT];

% cycle for all the remaing cells
% for i_cell = 1:size(indexMap,1)
%     if cell_type(i_cell) ~=0
%         %find distance between maximum force peak and maximum calcium activity
%         %in the mean signal to find the shift of integration
%         sign_medio = mean(Matrix{i_cell});
%         if cell_type(i_cell) == 1 %excitatory neuron
%             [~,i_m] = max(sign_medio);
%         elseif cell_type(i_cell) == 2 % inhibitory neuron
%             [~,i_m] = min(sign_medio);
%         end
%         dist = i_m/fScalcium -MENO;
%         Matrix_allFiles = [];
%         for i_file = 2:2:size(indexMap,2)
%             if indexMap(i_cell,i_file) ~= 0
%                 force = Data{1,i_file}.Fz.data;
%                 bicep = Data{2,i_file}.EMGfilt(:,1);
%                 tricep = Data{2,i_file}.EMGfilt(:,2);
%                 n_syn_cell = indexMap(i_cell,i_file);
%                 calcium = Data{4,i_file}.cells_signal(:,n_syn_cell);
%                 n_int = floor(length(force)/INTf)-n_int_rem;
%                 Area = NaN(n_int,4);
%                 for i_index = 1:n_int
%                     linef = -force(1+(i_index-1)*INTf:(i_index-1)*INTf+INTf);
%                     linec = calcium(1+dist+(i_index-1)*INTc:dist+(i_index-1)*INTc+INTc);
%                     lineb = bicep(1+(i_index-1)*INTe:(i_index-1)*INTe+INTe);
%                     linet = tricep(1+(i_index-1)*INTe:(i_index-1)*INTe+INTe);
%                     Area(i_index,1) = trapz(Xf,linef);
%                     Area(i_index,2) = trapz(X,linec);
%                     Area(i_index,3) = trapz(Xe,lineb);
%                     Area(i_index,4) = trapz(Xe,linet);
%                 end
%                 Area(:,2) = Area(:,2)-mean(Area(:,2));
%             else
%                 Area = [];
%             end
%         Matrix_allFiles = [Matrix_allFiles;Area];
%         end
% %         figure
% %         scatter (Matrix_allFiles(:,1),Matrix_allFiles(:,2))
% %         xlabel('Integration of force peaks')
% %         ylabel('Integration of calcium peaks')
% %         
% %         figure
% %         scatter (Matrix_allFiles(:,3),Matrix_allFiles(:,2))
% %         xlabel('Integration of Biceps')
% %         ylabel('Integration of calcium peaks')
% %         
% %         figure
% %         scatter (Matrix_allFiles(:,1),Matrix_allFiles(:,3))
% %         xlabel('Integration of force peaks')
% %         ylabel('Integration of Biceps')
%     end
% end

%%
% Graph F calcium only in the interesting part cell by cell
% cycle for all the remaing cells
% for i_cell = 1:size(indexMap,1)
%     if cell_type(i_cell) ~=0
%         Matrix_allFiles = [];
%         for i_file = 2:2:size(indexMap,2)
%             if indexMap(i_cell,i_file) ~= 0
%                 force = Data{1,i_file}.Fz.data;
%                 force = resample(force,20,100);
%                 bicep = Data{2,i_file}.EMGfilt(:,1);
%                 tricep = Data{2,i_file}.EMGfilt(:,2);
%                 n_syn_cell = indexMap(i_cell,i_file);
%                 calcium = Data{4,i_file}.cells_signal(:,n_syn_cell);
%                 switch in_case
%                     case 1 % index = force peaks
%                         index = Data{1,i_file}.Analysis.Fzpeaks(:,1); %index of the force peak
%                         if act==2
%                             phase = Data{1,i_file}.Analysis.Fzpeaks(:,15);
%                             index(phase<1.5)=[]; % select only index during pulling phase
%                                         % because we are working in the
%                                         % half tasks
%                         end
%                         index = round(index/fSrobot*fScalcium);
%                     case 2 % index = start movement
%                         if act == 1
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
%                             status = Data{1,i_file}.T_status.data;
%                             inin = find(status==2);
%                             index = inin(find (diff(inin)>1.5)+1);
%                             index = [inin(1),index];
%                             index = round(index/fSrobot*fScalcium);
%                         else
%                             index = Data{3,i_file}.pks{3};
%                             index = round(index*fScalcium);
%                         end
%                 end
%                 n_int = length(index);
%                 Area = NaN(n_int,4);
%                 for i_index = 1:n_int
%                     linef = -force(index(i_index)-INTf/2:index(i_index)+INTf/2-1);
%                     linec = calcium(index(i_index):index(i_index)+INTc-1);
%                     lineb = bicep(index(i_index)/fScalcium*fSemg-INTe/2:index(i_index)/fScalcium*fSemg+INTe/2-1);
%                     linet = tricep(index(i_index)/fScalcium*fSemg-INTe/2:index(i_index)/fScalcium*fSemg+INTe/2-1);
%                     Area(i_index,1) = trapz(Xf,linef);
%                     Area(i_index,2) = trapz(X,linec);
%                     Area(i_index,3) = trapz(Xe,lineb);
%                     Area(i_index,4) = trapz(Xe,linet);
%                 end
%                 Area(:,2) = Area(:,2)-mean(Area(:,2));
%             else
%                 Area = [];
%             end
%         Matrix_allFiles = [Matrix_allFiles;Area];
%         end
%         figure
%         scatter (Matrix_allFiles(:,1),Matrix_allFiles(:,2))
%         xlabel('Integration of force peaks')
%         ylabel('Integration of calcium peaks')
%         
%         figure
%         scatter (Matrix_allFiles(:,3),Matrix_allFiles(:,2))
%         xlabel('Integration of Biceps')
%         ylabel('Integration of calcium peaks')
%         
%         figure
%         scatter (Matrix_allFiles(:,1),Matrix_allFiles(:,3))
%         xlabel('Integration of force peaks')
%         ylabel('Integration of Biceps')
%     end
% end
%%
% Graph F calcium with all units (evaluation of the total activity)
% 
% Matrix_allFiles = [];
% for i_file = act:2:size(indexMap,2)
%     force = Data{1,i_file}.Fz.data;
%     bicep = Data{2,i_file}.EMGfilt(:,1);
%     tricep = Data{2,i_file}.EMGfilt(:,2);
%     n_int = floor(length(force)/INTf)-n_int_rem;
%     % consider only the part of the signal inside good trials
%     good_trials = Data{5,i_file}(:,1);
%     trial = Data{1,i_file}.cicles.data;
%     Area = NaN(n_int,4);
%     for i_index = 1:n_int
%         linef = -force(1+(i_index-1)*INTf:(i_index-1)*INTf+INTf);
%         Areac = 0;
%         for i_cell = 1:size(indexMap,1)
%             if indexMap(i_cell,i_file) ~= 0
%             if cell_type(i_cell) == 1 %excitatory neuron
%                 sign_medio = mean(Matrix{i_cell});
%                 [~,i_m] = max(sign_medio);
%                 dist = i_m/fScalcium -MENO;
%                 n_syn_cell = indexMap(i_cell,i_file);
%                 calciump = Data{4,i_file}.cells_signal(:,n_syn_cell);
%                 calcium = calciump-min(calciump);
%                 linec = calcium(1+dist+(i_index-1)*INTc:dist+(i_index-1)*INTc+INTc);
%                 Areac = Areac+trapz(X,linec);
%             end
%             end
%         end
%         lineb = bicep(1+(i_index-1)*INTe:(i_index-1)*INTe+INTe);
%         linet = tricep(1+(i_index-1)*INTe:(i_index-1)*INTe+INTe);
%         trial_int = unique(trial(1+(i_index-1)*INTf:(i_index-1)*INTf+INTf));
%         if ismember(trial_int,good_trials) % if it is not part of the good trials; it isn't considered
%             Area(i_index,1) = abs(trapz(Xf,linef));
%             Area(i_index,2) = Areac;
%             Area(i_index,3) = abs(trapz(Xe,lineb));
%             Area(i_index,4) = abs(trapz(Xe,linet));
%         else
%             Area(i_index,1) = NaN;
%             Area(i_index,2) = NaN;
%             Area(i_index,3) = NaN;
%             Area(i_index,4) = NaN;
%         end
%     end
%     Area(:,2) = Area(:,2)-nanmean(Area(:,2)); % units are not always with the same fluorescence during different days
%     Matrix_allFiles = [Matrix_allFiles;Area];
% end
% 
% % Remove outliers
% % for i_col = 1 : size(Matrix_allFiles,2)
% %     in_out = isoutlier(Matrix_allFiles(:,i_col));
% %     Matrix_allFiles(in_out==1,:) = NaN;
% % end
% 
% in_out = isoutlier(Matrix_allFiles(:,2));
% Matrix_allFiles(in_out==1,:) = [];
% 
% % This is only to visualize only positive values because it is more
% % intuitive
% Matrix_allFiles(:,2) = Matrix_allFiles(:,2)-min(Matrix_allFiles(:,2));
% 
% figure
% scatter (Matrix_allFiles(:,1),Matrix_allFiles(:,2))
% xlabel('Integration of force peaks')
% ylabel('Integration of calcium peaks')
%         
% figure
% scatter (Matrix_allFiles(:,3)+ Matrix_allFiles(:,4),Matrix_allFiles(:,2))
% xlabel('Integration of EMG')
% ylabel('Integration of calcium peaks')
% 
% figure
% scatter (Matrix_allFiles(:,1),Matrix_allFiles(:,3)+ Matrix_allFiles(:,4))
% xlabel('Integration of force peaks')
% ylabel('Integration of EMG')

%%
% Graph only of the area around onset peaks
Matrix_allFiles = [];
for i_file = act:2:size(indexMap,2)
    force = Data{1,i_file}.Fz.data;
    bicep = Data{2,i_file}.EMGfilt(:,1);
    tricep = Data{2,i_file}.EMGfilt(:,2);
    switch in_case
        case 1 % index = force peaks
           index = Data{1,i_file}.Analysis.Fzpeaks(:,1); %index of the force peak
           if act==2
               phase = Data{1,i_file}.Analysis.Fzpeaks(:,15);
               index(phase<1.5)=[]; % select only index during pulling phase
                                        % because we are working in the
                                        % half tasks
           end
           index = round(index/fSrobot*fScalcium);
        case 2 % index = start movement
            if act == 1
%                         pos_tot = Data{1,i_file}.pos_Spindle_drive.data;
%                         speed = derivative(pos_tot',1/fSrobot);
%                         acc = derivative(speed,1/fSrobot);
%                         th = mean(speed)+std(speed);
%                         [~,index] = findpeaks_GUI(-speed,'MINPEAKHEIGHT',th,'MINPEAKDISTANCE',10);
%                         onset = NaN(1,length(index));
%                         for n_p = 1:length(index)
%                             punti_neg = find(acc(1:index(n_p))<0 & speed(1:index(n_p))<th);
%                             if ~isempty(punti_neg)
%                                 onset(n_p) = punti_neg(end);
%                             end
%                         end
                status = Data{1,i_file}.T_status.data;
                inin = find(status==2);
                index = inin(find (diff(inin)>1.5)+1);
                index = [inin(1),index];
                index = round(index/fSrobot*fScalcium);
            else
                index = Data{3,i_file}.pks{3};
                index = round(index*fScalcium);
            end
    end
    n_int = length(index);
    Area = NaN(n_int,4);
    for i_index = 1:n_int
        linef = -force(index(i_index)/fScalcium*fSrobot-INTf/2:index(i_index)/fScalcium*fSrobot+INTf/2-1);
        Areac = 0;
        for i_cell = 1:size(indexMap,1)
            if indexMap(i_cell,i_file) ~= 0
            if cell_type(i_cell) == 1 %excitatory neuron
                n_syn_cell = indexMap(i_cell,i_file);
                calciump = Data{4,i_file}.cells_signal(:,n_syn_cell);
                calcium = calciump-min(calciump);
                linec = calcium(index(i_index):index(i_index)+INTc-1);
                Areac = Areac+trapz(X,linec);
            end
            end
        end
        lineb = bicep(index(i_index)/fScalcium*fSemg-INTe/2:index(i_index)/fScalcium*fSemg+INTe/2-1);
        linet = tricep(index(i_index)/fScalcium*fSemg-INTe/2:index(i_index)/fScalcium*fSemg+INTe/2-1);
        Area(i_index,1) = abs(trapz(Xf,linef));
        Area(i_index,2) = Areac;
        Area(i_index,3) = abs(trapz(Xe,lineb));
        Area(i_index,4) = abs(trapz(Xe,linet));
    end
    Area(:,2) = Area(:,2)-nanmean(Area(:,2)); % units are not always with the same fluorescence during different days
    Matrix_allFiles = [Matrix_allFiles;Area];
end
% Remove outliers
% for i_col = 1 : size(Matrix_allFiles,2)
%     in_out = isoutlier(Matrix_allFiles(:,i_col));
%     Matrix_allFiles(in_out==1,:) = NaN;
% end

in_out = isoutlier(Matrix_allFiles(:,2));
Matrix_allFiles(in_out==1,:) = [];

% This is only to visualize only positive values because it is more
% intuitive
Matrix_allFiles(:,2) = Matrix_allFiles(:,2)-min(Matrix_allFiles(:,2));

figure
scatter (Matrix_allFiles(:,1),Matrix_allFiles(:,2))
xlabel('Integration of force peaks')
ylabel('Integration of calcium peaks')
        
figure
scatter (Matrix_allFiles(:,3)+ Matrix_allFiles(:,4),Matrix_allFiles(:,2))
xlabel('Integration of EMG')
ylabel('Integration of calcium peaks')

figure
scatter (Matrix_allFiles(:,1),Matrix_allFiles(:,3)+ Matrix_allFiles(:,4))
xlabel('Integration of force peaks')
ylabel('Integration of EMG')
end