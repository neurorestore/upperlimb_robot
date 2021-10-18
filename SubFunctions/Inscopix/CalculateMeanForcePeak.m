% Function to create a matrix around the index of the force value for all
% the files
% act is the index that specify if we are evaluating half or active
% task(1->active, 2->half)
function [Mforce_allFiles,days_force] = CalculateMeanForcePeak(Data,in_case,MENO,PIU,fSrobot,fScalcium,act)
Mforce_allFiles= [];
days_force = [];
for i_file = act:2:size(Data,2)
        force = Data{1,i_file}.Fz.data;
        %find index to align force peaks
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
                    index = Data{3,i_file}.start;
                    index = round(index*fScalcium);
%                     if act == 1
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
%                     else
%                         index = Data{3,i_file}.pks{3};
%                         index = round(index*fScalcium);
%                     end

            end
            %create matrix with signal force
            Mforce = NaN(length(index),round(MENO*fSrobot+PIU*fSrobot+1));
            %cycle for all the index
            for i_ind = 1:length(index)
                st = round(index(i_ind)/fScalcium*fSrobot-MENO*fSrobot);
                en = round(index(i_ind)/fScalcium*fSrobot+PIU*fSrobot);
                if st>1 && en<length(force)
                    linec = force(st:en);
                    Mforce(i_ind,:) = linec;
                end
            end
            line_day = ones(length(index),1)*i_file;
            days_force = [days_force;line_day];
            Mforce_allFiles = [Mforce_allFiles;Mforce];
end

%plot force signal 
figure;hold on; set(gca, 'FontSize', 14)
hold on
x = [-MENO:1/fSrobot:PIU];
plot(x, nanmean(Mforce_allFiles),'LineWidth',2)
Mf_dev = nanstd(Mforce_allFiles)/sqrt(size(Mforce_allFiles,1));
X=[x, fliplr(x)];
Y=[nanmean(Mforce_allFiles) + Mf_dev,fliplr(nanmean(Mforce_allFiles) - Mf_dev)];
fill( X,Y,'b');
alpha(.10)
ylim([-1.7,0.1])
line([0 0],get(gca, 'YLim'), 'Color', [0 0 0], 'LineStyle', '--')
xlim([-MENO,PIU])
%ylim(get(gca, 'YLim'))
xlabel('Time (s)')
ylabel('Amplitude Force z')

end