%% Function to calculate parameters extracted by the trajectory 
% 1. Calculate when there is a movement
% 2. Associate every movement a one force peaks (only for half task)
% 3. Count number of submovements, attempts, v media of peaks during each trial
% 4. Study the trajectory in x-y space (smoothness, vmedia, height, AUC,
%total length, y mov in the retraction phase...)

function Data = CalculateKINparam(Data)

% Find mov peaks in the main direction
x = sgolayfilt(Data.SIMI.x,3,21);
y = sgolayfilt(Data.SIMI.y,3,21);
fS = Data.SIMI.fS_KIN;
fSrobot = Data.Recorded_Data.fS_robot;
pks = cell(size(x,2),1);
for n_points = 1:1%size(x,2)
    if sum(isnan(x))<1
        if strcmp (Data.info.Status,'Push_active')
            speed = -derivative(x(:,n_points),1/fS);
        elseif strcmp (Data.info.Status,'Pull_active')
            speed = derivative(x(:,n_points),1/fS);
        elseif strcmp (Data.info.Status,'active')
            speed = -derivative(x(:,n_points),1/fS);
        end

        acc = derivative(speed,1/fS);
    else
        index = find(isnan(x)==0);
        if strcmp (Data.info.Status,'Push_active')
            speed = -derivative(x(index(1):end,n_points),1/fS);
        elseif strcmp (Data.info.Status,'Pull_active')
            speed = derivative(x(index(1):end,n_points),1/fS);
        elseif strcmp (Data.info.Status,'active')
            speed = -derivative(x(index(1):end,n_points),1/fS);
        end
        acc = derivative(speed,1/fS);
        add = NaN(index(1)-1,1);
        speed = [add;speed];
        acc = [add;acc];
    end
    th = nanmean(speed)+1/2*nanstd(speed);
    [pks{n_points,1},pks{n_points,2}] = findpeaks_GUI(speed,'MINPEAKHEIGHT',th,'MINPEAKDISTANCE',10);
    % find onset pks speed -> points where the acceleration change
    % direction before the threshold
    onset = NaN(1,length(pks{n_points,2}));
    for n_p = 1:length(pks{n_points,2})
        punti_neg = find(acc(1:pks{n_points,2}(n_p))<0 & speed(1:pks{n_points,2}(n_p))<th);
        if ~isempty(punti_neg)
            onset(n_p) = punti_neg(end)+1;
        end
    end
    pks{n_points,3} = onset/fS;       
    pks{n_points,2} = pks{n_points,2}/fS;
    
    figure
    hold on
    plot(speed)
    plot(x(:,n_points)-nanmean(x(:,n_points)))
    scatter(pks{n_points,2}*fS,pks{n_points,1})
    onset(isnan(onset))=[];
    scatter(onset,x(onset,n_points)-nanmean(x(:,n_points)))
end

Data.SIMI.pks = pks;
Data.SIMI.name_pks = {'Speed Amplitude','Index Peak','Index onset'};

% check if pks speed are in good trials
[SIMI] = InsideGoodTrialsSpeed (Data.SIMI, Data.Recorded_Data.cicles.data, Data.good_trials, Data.Recorded_Data.T_status.data,fSrobot,fS);
Data.SIMI = SIMI;

%% Associate movement to force peaks to divide attempts and submov
%it is possible only if it is a half-half task
if ~strcmp (Data.info.Status,'active')
    if strcmp (Data.info.Status,'Push_active')
        force_pks = Data.Recorded_Data.Analysis.Fzpeaks(:,1)/Data.Recorded_Data.fS_robot;
    elseif strcmp (Data.info.Status,'Pull_active')
        force_pks = Data.Recorded_Data.Analysis.Fzpeaks_contra(:,1)/Data.Recorded_Data.fS_robot;
    end
    speed_pks = Data.SIMI.pks{1,2};
    ind_pks = ones(length(speed_pks),1); %vector where to save if it creates movements or not 
    ind_out_force_pks = NaN(length(speed_pks),1);
    dist = NaN(length(speed_pks),1); % vector where to save the distance between index force peak and index speed peak
    for i_pks = 1:length(speed_pks)
        in_for = find (force_pks<speed_pks(i_pks)+0.1);
        if ~isempty(in_for)
            %ind_pks(in_for(end)) = 1;
            dist(i_pks) = speed_pks(i_pks)-force_pks(in_for(end));
            ind_out_force_pks(i_pks) = in_for(end);
        else
            ind_pks(i_pks) = 0;
        end
    end
    ind_pks(dist>0.5)= 0; % check if some peaks has been associated with a wrong and far force peak
    ind_out_force_pks (dist>0.5)= NaN;
    dist (dist>0.5) = NaN;
    ind_out_force_pks(isnan(ind_out_force_pks))=[];
    ind_force_pks = zeros(length(force_pks),1);
    ind_force_pks(ind_out_force_pks) = 1;
    if strcmp (Data.info.Status,'Push_active')
        Data.Recorded_Data.Analysis.Fzpeaks = [Data.Recorded_Data.Analysis.Fzpeaks,ind_force_pks];
    elseif strcmp (Data.info.Status,'Pull_active')
        Data.Recorded_Data.Analysis.Fzpeaks_contra = [Data.Recorded_Data.Analysis.Fzpeaks_contra,ind_force_pks];
    end
    Data.Recorded_Data.Analysis.name_peaks_matrix{1,size(Data.Recorded_Data.Analysis.Fzpeaks,2)}= 'create_mov';
    % add column to specify which speed pks are associated with a force
    % peak
    % This is important especially after injury when animals are not able
    % to finish the movement so the slide is bring back to the home
    % position manually
    Data.SIMI.pks{1,size(Data.SIMI.pks,2)+1} = ind_pks';
    Data.SIMI.name_pks{1,size(Data.SIMI.pks,2)} = 'associate_force_pks';
end

%% count submovements and attempts and calculate vmedia
if strcmp (Data.info.Status,'Push_active')
    fz_pks = Data.Recorded_Data.Analysis.Fzpeaks;
    submov = NaN(size(Data.good_trials,1),1);
    attem = NaN(size(Data.good_trials,1),1);
    vm_pks = NaN(size(Data.good_trials,1),1);
    for i_t = 1:size(Data.good_trials,1)
        tot_pks = find(fz_pks(:,14)==Data.good_trials(i_t,1) & fz_pks(:,15)==2);% only peaks in the pulling phase
        submov(i_t) = sum(fz_pks(tot_pks,16)); % this coloumn is one for mov peak and zero for no mov
        attem(i_t) = length(tot_pks)-submov(i_t);
        tot_pks_speed = find (Data.SIMI.pks{1,4} == Data.good_trials(i_t,1));
        vm_pks(i_t) = mean(Data.SIMI.pks{1,1}(tot_pks_speed)& Data.SIMI.pks{1,6}(tot_pks_speed)==1);
    end
elseif strcmp (Data.info.Status,'Pull_active')
    fz_pks = Data.Recorded_Data.Analysis.Fzpeaks_contra;
    submov = NaN(size(Data.good_trials,1),1);
    attem = NaN(size(Data.good_trials,1),1);
    vm_pks = NaN(size(Data.good_trials,1),1);
    for i_t = 1:size(Data.good_trials,1)
        tot_pks = find(fz_pks(:,14)==Data.good_trials(i_t,1) & fz_pks(:,15)==1);% only peaks in the pushing phase
        submov(i_t) = sum(fz_pks(tot_pks,16)); % this coloumn is one for mov peak and zero for no mov
        attem(i_t) = length(tot_pks)-submov(i_t);
        tot_pks_speed = find (Data.SIMI.pks{1,4} == Data.good_trials(i_t,1));
        vm_pks(i_t) = mean(Data.SIMI.pks{1,1}(tot_pks_speed)& Data.SIMI.pks{1,6}(tot_pks_speed)==1);
    end
elseif strcmp (Data.info.Status,'active')
    submov = NaN; %motors move the robot, animals can't do submovements
    attem = NaN;
    vm_pks = NaN(size(Data.good_trials,1),1);
    for i_t = 1:size(Data.good_trials,1)
        tot_pks_speed = find (Data.SIMI.pks{1,4} == Data.good_trials(i_t,1));
        vm_pks(i_t) = mean(Data.SIMI.pks{1,1}(tot_pks_speed));
    end     
end
Data.SIMI.Analysis.submov = submov;
Data.SIMI.Analysis.attempts = attem;
Data.SIMI.Analysis.vm_pks = vm_pks;

%% Define the trajectory cycle by cycle
% find the longest trajectory to create a matrix

max_l = floor(max(Data.good_trials(:,3)-Data.good_trials(:,2))*fS)+1;
trajectory = NaN(max_l,2,size(Data.good_trials,1));
xx = sgolayfilt(Data.SIMI.x,3,51);
yy = sgolayfilt(Data.SIMI.y,3,51);
start_pulling = NaN(size(Data.good_trials,1),1);
for i_cy = 1:size(Data.good_trials,1)
    st = round(Data.good_trials(i_cy,2)*fS);
    % in the active task the end of the cycle is the end of the cycle
    if strcmp (Data.info.Status,'active')
        en = round(Data.good_trials(i_cy,3)*fS);
    elseif strcmp (Data.info.Status,'Push_active')
        % in the push active task the cycle ends with the last speed peaks
        % we don't consider the part of the movement manually push by the
        % operator when the animal is not able to finish the task
        index_speed_pks = find (Data.SIMI.pks{1,4}==Data.good_trials(i_cy,1) & Data.SIMI.pks{1,6}==1);
        if ~isempty(index_speed_pks)
            add = Data.SIMI.pks{1,2}(index_speed_pks(end))-Data.SIMI.pks{1,3}(index_speed_pks(end));
            en = round((Data.SIMI.pks{1,2}(index_speed_pks(end))+ add)*fS);
        else
            en = round(Data.good_trials(i_cy,3)*fS);
        end
    elseif strcmp (Data.info.Status,'Pull_active')
        en = round(Data.good_trials(i_cy,3)*fS);
    end
    trajectory(1:en-st+1,1,i_cy) = xx(st:en)-xx(st);
    trajectory(1:en-st+1,2,i_cy) = -(yy(st:en)-yy(st));
    % plot to check trajectories
    figure
    plot(trajectory(:,1,i_cy),trajectory(:,2,i_cy))
    title (num2str(i_cy))
    %ylim([0 15])
    % ADD FUNCTION MEAN PLOT
    % create vector info to divide pushing and pulling phase
    phase_cycle = Data.Recorded_Data.T_status.data(st*fSrobot/fS:en*fSrobot/fS);
    pulling_phase = find(phase_cycle==2);
    continuo = find(diff(pulling_phase)>1.5);
    if isempty(continuo)
        start_pulling(i_cy) = round(pulling_phase(1)/fSrobot*fS);
    else
        start_pulling(i_cy) = round(pulling_phase(continuo(1)+1)/fSrobot*fS);
    end
end

Data.SIMI.trajectory = trajectory;
%% Define the trajectory cycle by cycle for top view
% find the longest trajectory to create a matrix
if 0
    max_l = floor(max(Data.good_trials(:,3)-Data.good_trials(:,2))*fS)+1;
    trajectory = NaN(max_l,2,size(Data.good_trials,1));
    xx = sgolayfilt(Data.SIMI.x2,3,51);
    yy = sgolayfilt(Data.SIMI.z,3,51);
    start_pulling = NaN(size(Data.good_trials,1),1);
    for i_cy = 1:size(Data.good_trials,1)
        st = round(Data.good_trials(i_cy,2)*fS);
        % in the active task the end of the cycle is the end of the cycle
        if strcmp (Data.info.Status,'active')
            en = round(Data.good_trials(i_cy,3)*fS);
        elseif strcmp (Data.info.Status,'Push_active')
            % in the push active task the cycle ends with the last speed peaks
            % we don't consider the part of the movement manually push by the
            % operator when the animal is not able to finish the task
            index_speed_pks = find (Data.SIMI.pks{1,4}==Data.good_trials(i_cy,1) & Data.SIMI.pks{1,6}==1);
            if ~isempty(index_speed_pks)
                add = Data.SIMI.pks{1,2}(index_speed_pks(end))-Data.SIMI.pks{1,3}(index_speed_pks(end));
                en = round((Data.SIMI.pks{1,2}(index_speed_pks(end))+ add)*fS);
            else
                en = round(Data.good_trials(i_cy,3)*fS);
            end
        elseif strcmp (Data.info.Status,'Pull_active')
            en = round(Data.good_trials(i_cy,3)*fS);
        end
        trajectory(1:en-st+1,1,i_cy) = xx(st:en)-xx(st);
        trajectory(1:en-st+1,2,i_cy) = -(yy(st:en)-yy(st));
        % plot to check trajectories
        figure
        plot(trajectory(:,1,i_cy),trajectory(:,2,i_cy))
        title (num2str(i_cy))
        %ylim([0 15])
        % ADD FUNCTION MEAN PLOT
        % create vector info to divide pushing and pulling phase
        phase_cycle = Data.Recorded_Data.T_status.data(st*fSrobot/fS:en*fSrobot/fS);
        pulling_phase = find(phase_cycle==2);
        continuo = find(diff(pulling_phase)>1.5);
        if isempty(continuo)
            start_pulling(i_cy) = round(pulling_phase(1)/fSrobot*fS);
        else
            start_pulling(i_cy) = round(pulling_phase(continuo(1)+1)/fSrobot*fS);
        end
    end

    Data.SIMI.trajectory2 = trajectory;
end
%% Extract parameters for every trajectory
Data = CalculateTRAJparam(Data,start_pulling);

end








