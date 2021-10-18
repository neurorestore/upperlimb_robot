%% Function to create a graph of the force toward the EMG activity
% Select only parts where there is a good trial and then divide them in
% useful intervals

function [R,Rt,Rd,Area_cy] = PlotEMGvsF(fz,EMG,good_trials,fSrobot,fSEMG,INT,dist,filename)

INTf = INT*fSrobot;
INTe = INT*fSEMG;
name_EMG = {'Biceps','Triceps','EDC','ECR','FCU','Deltoid' ,'trigger'};

Xf = [1/fSrobot:1/fSrobot:INT];
Xe = [1/fSEMG:1/fSEMG:INT];

filename_new = filename(1:end-17);
dist = round(nanmean(dist)/fSEMG*fSrobot);
if isnan(dist)
    dist = zeros(size(EMG,2),1);
end

for i_emg = 1:size(EMG,2)
    if isnan(dist(i_emg)) % to fix the problem when there is not a distance because there are not enough peaks
        dist(i_emg) = 0;
    end
end
% cycle for different EMG channels
Area = [];
% create a variable where to save mean of the force and emg for cycle (as
% parameters)
Area_cy = NaN(size(good_trials,1),size(EMG,2)+1); % one column for the force and one for every EMG
for i_emg = 1:size(EMG,2)
    Area_ch = [];
    for i_cy = 1:size(good_trials,1)-1 % cycle for every good trials part
        fz_cy = fz(good_trials(i_cy,2)*fSrobot+dist(i_emg):good_trials(i_cy,3)*fSrobot+dist(i_emg));
        emg_cy = EMG(good_trials(i_cy,2)*fSEMG:good_trials(i_cy,3)*fSEMG,i_emg);
        n_int = floor((good_trials(i_cy,3)-good_trials(i_cy,2))/INT);
        energy_cy = NaN(n_int,2);
        for i_int = 1:n_int % cycle for division in intervals
            fz_int = fz_cy((1+(i_int-1)*INTf):((i_int-1)*INTf+INTf));
            emg_int = emg_cy((1+(i_int-1)*INTe):((i_int-1)*INTe+INTe));
            energy_cy(i_int,1) = abs(trapz(Xf,fz_int)); % calculate the integral area in each interval
            energy_cy(i_int,2) = abs(trapz(Xe,emg_int));
        end
        Area_ch = [Area_ch;energy_cy];
        Area_cy (i_cy,1) = nanmean(energy_cy(:,1))/(n_int*INT);
        Area_cy (i_cy,i_emg+1) = nanmean(energy_cy(:,2))/(n_int*INT);
    end
    if i_emg == 1
        Area = [Area,Area_ch];
    else
        Area = [Area,Area_ch(:,2)];
    end
end

%% Plot and linear relation 
Area(isoutlier(Area(:,2)),:) = [];
Area(isoutlier(Area(:,1)),:) = [];
Area(isoutlier(Area(:,3)),:) = [];
[R,P] = corrcoef(Area(:,1:2));
% Triceps
[Rt,Pt] = corrcoef(Area(:,1:2:3));
% Others EMG
Rd = [];
for i = 4:size(Area,2)
    [Rr,Pr] = corrcoef(Area(:,1:i-1:i));
    Rd = [Rd;Rr(1,2)];
end

% Plot
for i_fig = 2:size(Area,2)
    [p,S] = polyfit(Area(:,1),Area(:,i_fig),1);
    x = [0:0.01:max(Area(:,1))+0.05];
    [y_fit,delta] = polyval(p,x,S);
    figure; hold on
    set(gca,'FontSize',14)
    plot(Area(:,1),Area(:,i_fig),'ko','LineWidth',1)
    hold on
    plot(x,y_fit,'Color',[0 113 188]/255,'LineWidth',1.5)
    plot(x,y_fit+2*delta,'b--',x,y_fit-2*delta,'b--')
    xlim ([0,max(Area(:,1))+0.05])
    xlabel('Force Integration')
    ylabel(name_EMG{i_fig-1})
    if i_fig == 2
        title (['CorrCoeff ',num2str(R(1,2)),', p = ',num2str(P(1,2))])
        savefig([filename_new,'BicForceRelation.fig'])
    end
end
 
end







