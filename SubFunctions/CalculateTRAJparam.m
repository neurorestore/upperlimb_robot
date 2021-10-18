%% Function to extract parameters from the single trajectories
% every parameters is calculated for the whole trajectory, for the pushing
% phase and for the pulling phase

function Data = CalculateTRAJparam(Data,start_pulling)

trials = Data.good_trials(:,1);
trajectories = Data.SIMI.trajectory;
fS = Data.SIMI.fS_KIN;

%% initialize all parameters

height = NaN(length(trials),3);
amplitude = NaN(length(trials),3);
length_tr = NaN(length(trials),3);
length_x = NaN(length(trials),3);
length_y = NaN(length(trials),3);
AUC = NaN(length(trials),1);

speed = NaN(length(trials),3);
speed_x = NaN(length(trials),3);
speed_y = NaN(length(trials),3);
smoothness_th = NaN(length(trials),3);
smoothness_x = NaN(length(trials),3);
smoothness_y = NaN(length(trials),3);
speed_max = NaN(length(trials),3);
speed_max_x = NaN(length(trials),3);
speed_max_y = NaN(length(trials),3);

acceleration = NaN(length(trials),3);
acceleration_max = NaN(length(trials),3);

%% calculate length parameters
%(height,length, x length, y length, amplitude(x_max-x_min), area curve)
for i_cy = 1:length(trials)
   tr = trajectories(:,:,i_cy); 
   tr(isnan(tr(:,1)),:)= [];
   
   for i_ph = 1:3 % cycle to analyse all trajectory or only pull or push phase
       switch i_ph
           case 1 % parameter for all trajectory
               st = 1;
               en = length(tr);
           case 2 % parameter for push phase
               st = 1;
               en = start_pulling(i_cy);
           case 3 % parameter for pull phase
               st = start_pulling(i_cy);
               en = length(tr);
       end
       height(i_cy,i_ph) = max(tr(st:en,2))-min(tr(st:en,2));
       amplitude(i_cy,i_ph) = max(tr(st:en,1))-min(tr(st:en,1));
       % length
       length_tr (i_cy,i_ph) = 0.0;
       length_x (i_cy,i_ph) = 0.0;
       length_y (i_cy,i_ph) = 0.0;
       for i = st:en-1
           length_step = sqrt( (tr(i+1,1)-tr(i,1))^2 + (tr(i+1,2)-tr(i,2))^2 );
           length_tr (i_cy,i_ph) = length_tr(i_cy,i_ph) + length_step;
           length_x (i_cy,i_ph) = length_x(i_cy,i_ph) + abs(tr(i+1,1)-tr(i,1));
           length_y (i_cy,i_ph) = length_y(i_cy,i_ph) + abs(tr(i+1,2)-tr(i,2));
       end
       
   end
   AUC (i_cy) = polyarea(tr(:,1),tr(:,2));
end

%% calculate speed parameters
%(mean speed, smoothness, max speed -total, x and y direction-)

for i_cy = 1:length(trials)
   tr = trajectories(:,:,i_cy); 
   tr(isnan(tr(:,1)),:)= [];
   
   for i_ph = 1:3 % cycle to analyse all trajectory or only pull or push phase
       switch i_ph
           case 1 % parameter for all trajectory
               st = 1;
               en = length(tr);
           case 2 % parameter for push phase
               st = 1;
               en = start_pulling(i_cy);
           case 3 % parameter for pull phase
               st = start_pulling(i_cy);
               en = length(tr);
       end
       % speed
       vel_step = NaN(en-st,1);
       vel_step_x = NaN(en-st,1);
       vel_step_y = NaN(en-st,1);
       for i = st:en-1
           length_step = sqrt( (tr(i+1,1)-tr(i,1))^2 + (tr(i+1,2)-tr(i,2))^2 );
           vel_step(i-st+1) = length_step*fS;
           vel_step_x(i-st+1) = abs(tr(i+1,1)-tr(i,1))*fS; % it is the same if I use derivative
           vel_step_y(i-st+1) = abs(tr(i+1,2)-tr(i,2))*fS;
       end
       vel_step(isoutlier(vel_step,'mean','ThresholdFactor',5)) = NaN;
       vel_step_x(isoutlier(vel_step_x,'mean','ThresholdFactor',5)) = NaN;
       vel_step_y(isoutlier(vel_step_y,'mean','ThresholdFactor',5)) = NaN;
       speed (i_cy,i_ph)= nanmean(vel_step);
       speed_x (i_cy,i_ph)= nanmean(vel_step_x);
       speed_y (i_cy,i_ph)= nanmean(vel_step_y);
       speed_max (i_cy,i_ph) = max(vel_step);
       speed_max_x (i_cy,i_ph) = max(vel_step_x);
       speed_max_y (i_cy,i_ph) = max(vel_step_y);
       %th = nanstd(vel_step)/3;
       [smooth] = findpeaks_GUI(vel_step);%,'MINPEAKHEIGHT',th,'MINPEAKDISTANCE',10);
       smoothness_th (i_cy,i_ph) = length(smooth);
       %th = nanstd(vel_step_x)/3;
       [smooth] = findpeaks_GUI(vel_step_x);%,'MINPEAKHEIGHT',th,'MINPEAKDISTANCE',10);
       smoothness_x (i_cy,i_ph) = length(smooth);
       %th = nanstd(vel_step_y)/3;
       [smooth] = findpeaks_GUI(vel_step_y);%,'MINPEAKHEIGHT',th,'MINPEAKDISTANCE',10);
       smoothness_y (i_cy,i_ph) = length(smooth);
              
       %% calculate acceleration parameters
       % mean acceleration, max acceleration
       acc_step = NaN(en-st,1);
       
       for i = st:length(vel_step)-1
           acc_step(i-st+1) = (vel_step(i+1)-vel_step(i,1))*fS;
       end
       acc_step(isoutlier(acc_step,'mean','ThresholdFactor',5)) = NaN;
       acceleration(i_cy,i_ph) = nanmean(acc_step);
       acceleration_max (i_cy,i_ph) = max(acc_step);
   end
   
end

%% save all parameters in SIMI.Analysis struct

Data.SIMI.Analysis.height = height;
Data.SIMI.Analysis.amplitude = amplitude;
Data.SIMI.Analysis.length_tr = length_tr;
Data.SIMI.Analysis.length_x = length_x;
Data.SIMI.Analysis.length_y = length_y;
Data.SIMI.Analysis.AUC = AUC;

Data.SIMI.Analysis.speed = speed;
Data.SIMI.Analysis.speed_x = speed_x;
Data.SIMI.Analysis.speed_y = speed_y;
Data.SIMI.Analysis.smoothness_th = smoothness_th;
Data.SIMI.Analysis.smoothness_x = smoothness_x;
Data.SIMI.Analysis.smoothness_y = smoothness_y;
Data.SIMI.Analysis.speed_max = speed_max;
Data.SIMI.Analysis.speed_max_x = speed_max_x;
Data.SIMI.Analysis.speed_max_y = speed_max_y;

Data.SIMI.Analysis.acceleration = acceleration;
Data.SIMI.Analysis.acceleration_max = acceleration_max;

end