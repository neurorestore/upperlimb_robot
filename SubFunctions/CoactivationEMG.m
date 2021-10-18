%% Study the coactivation of different muscles
% To evaluate coactivation of muscles (Probability density distribution)
% 1. find good trials
% 2. calculate average amplitudes every 40 ms for two EMGs
% 3. Use these values to target x and y coordinates of a bin in a grid
% 4. The relative number of counts in each bin defined the probability
% density distribution

function [EMGav_tot,EMGav_pull,EMGav_push] = CoactivationEMG (EMG,fS_EMG,GoodTrials,t_status,path,file,fS_rob)

nameEMG = {'Biceps','Triceps','DE'};
n_win = fS_EMG*0.010; %20 ms every window
n_area = 30;
GoodTrials(:,2:3) = GoodTrials(:,2:3)*fS_EMG;
ratio = fS_EMG/fS_rob;

EMGav_tot = [];
EMGav_pull = [];
EMGav_push = [];
for n_trials = 1: size(GoodTrials,1) %cicles for all the good trials
    % coactivation during complete trials
    d = floor((GoodTrials(n_trials,3)-GoodTrials(n_trials,2))/n_win);
    add = find(t_status(GoodTrials(n_trials,2)/ratio+10:GoodTrials(n_trials,3)/ratio)>1.5);
    add = add(1)*10;
    EMGav = NaN(d,size (EMG,2));
    for n_EMG = 1: size (EMG,2)
        EMGp = abs(EMG(GoodTrials(n_trials,2):GoodTrials(n_trials,3),n_EMG));
        change_status = NaN;
        for n_d = 1:d
            EMGav(n_d,n_EMG) = median(EMGp(1+(n_d-1)*n_win:n_win+(n_d-1)*n_win)); % consider if it is better median or mean
            %EMGav(n_d,n_EMG) = trapz(EMGp(1+(n_d-1)*n_win:n_win+(n_d-1)*n_win));
            if (1+(n_d-1)*n_win)<add % division of the trial between push and pull
                change_status = n_d;
            end
        end
    end
    EMGav_tot = [EMGav_tot;EMGav];
    EMGav_pull = [EMGav_pull;EMGav(change_status+1:end,:)];
    EMGav_push = [EMGav_push;EMGav(1:change_status,:)];

end
v_max = max(EMGav_tot);
% find maximum value for every EMG but avoiding outlier
for n_EMG = 1:size(EMGav_tot,2)
    values = sort(EMGav_tot(:,n_EMG),'descend');
    v_out = isoutlier(values(1:4));
    index = find(v_out==0);
    v_max(n_EMG) = values(index(1));
end
figure; scatter (EMGav_tot(:,1)/v_max(1),EMGav_tot(:,2)/v_max(2))
xlim([0 0.5])
ylim([0 0.5])
figure; scatter (EMGav_pull(:,1)/v_max(1),EMGav_pull(:,2)/v_max(2))
xlim([0 0.5])
ylim([0 0.5])
figure; scatter (EMGav_push(:,1)/v_max(1),EMGav_push(:,2)/v_max(2))
xlim([0 0.5])
ylim([0 0.5])

% % plot of the 90% of the data because otherwise is too little the active
% % part
% ordered = sort(EMGav_tot(:,1));
% enne = round(size(EMGav_tot,1)*90/100);
% xl_1 = round(ordered(enne)*100);
% xl_1 = xl_1/100;
% ordered = sort(EMGav_tot(:,2));
% xl_2 = round(ordered(enne)*100);
% xl_2 = xl_2/100;
% ordered = sort(EMGav_tot(:,3));
% xl_3 = round(ordered(enne)*100);
% xl_3 = xl_3/100;
% 
% %step = v_max/n_area;
% figure
% lim = max([xl_1,xl_2,xl_3]);
%h=
histogram2(EMGav_tot(:,1)/v_max(1),EMGav_tot(:,2)/v_max(2),n_area,'XBinLimits',[0 0.3],'yBinLimits',[0,0.3],'DisplayStyle','tile');
colormap ('Jet')
ylabel(nameEMG{2})
xlabel(nameEMG{1})
%imagesc(h.BinCounts)

figure
%h=
histogram2(EMGav_pull(:,1)/v_max(1),EMGav_pull(:,2)/v_max(2),n_area,'XBinLimits',[0 0.3],'yBinLimits',[0,0.3],'DisplayStyle','tile');
colormap ('Jet')
ylabel(nameEMG{2})
xlabel(nameEMG{1})
saveas(gcf,[path,file,'_CoactivationPull.fig'])
%saveas(gcf,[path,file,'_CoactivationPull.png'])
%imagesc(flip(h.BinCounts'))

figure
%h=
histogram2(EMGav_push(:,1)/v_max(1),EMGav_push(:,2)/v_max(2),n_area,'XBinLimits',[0 0.3],'yBinLimits',[0,0.3],'DisplayStyle','tile')%,'ShowEmptyBins','on');
colormap ('Jet')
ylabel(nameEMG{2})
xlabel(nameEMG{1})
saveas(gcf,[path,file,'_CoactivationPush.fig'])
%saveas(gcf,[path,file,'_CoactivationPush.png'])
%imagesc(flip(h.BinCounts'))
end

