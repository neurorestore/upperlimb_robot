% Function to select only the good trials before continuing the analysis of
% the recordings, not all the recordings are good for all trials for
% multiple reason, this function create a vector where it is specified
% which trials have to be excluded in the future analysis

% !!!! RECORDINGS ARE NOT CUT in the middle to eliminate wrong trials
% to avoid filtering of cutted signals!!!!!

function Data = Cleaning_Trials(Data)

Data.Recorded_Data.cicles.data(1)=0; % the first value of the matrix recorded by cRIO is ugual at the last one of the previous file

% eliminate last trials, cut when one of the data finishes
timeEMG = size(Data.VICON.EMG,1)/Data.VICON.fS_EMG;
timeRobot = length(Data.Recorded_Data.t.data)/Data.Recorded_Data.fS_robot;
timeKIN = Data.SIMI.duration;

[tTot,i_min] = min([timeEMG,timeRobot,timeKIN]);
Data.tTot = tTot;
fields = fieldnames (Data.Recorded_Data);
d = floor (tTot*Data.Recorded_Data.fS_robot);
Data.VICON.EMG = Data.VICON.EMG(1:d*Data.VICON.fS_EMG/Data.Recorded_Data.fS_robot,:);
for i_field =1:length(fields)
    if isstruct(Data.Recorded_Data.(fields{i_field})) && ~isempty(Data.Recorded_Data.(fields{i_field}))
        Data.Recorded_Data.(fields{i_field}).data=Data.Recorded_Data.(fields{i_field}).data(1:d);        
    end
end

trials = unique(Data.Recorded_Data.cicles.data);
del_buffer =0;

% Correction VICON error: when VICON blocks, I add a vector of zero at the
% beginning of the EMG to synchronize but data are missing, so these trials
% are not good
zero_add = find(Data.VICON.EMG(:,1)==0);
if ~isempty(zero_add) && zero_add(1) ==1 && zero_add(2) ==1
    zero_add_cons = zero_add(diff(zero_add==1));
    if zero_add_cons(1)==1 && length(zero_add_cons)>100
        delete_buffer = (1:Data.Recorded_Data.cicles.data(zero_add_cons(end)));
        del_buffer = [del_buffer, delete_buffer];
    end
end


% eliminate trials that have been indicated as wrong during recordings
index_wrong = find(Data.Recorded_Data.Wrong_trials.data==1);
delete_buffer = unique(Data.Recorded_Data.cicles.data(index_wrong));
del_buffer = [del_buffer; delete_buffer'];

% final check visivo
figure
hold on
x = [0:1/Data.VICON.fS_EMG:size(Data.VICON.EMG,1)/Data.VICON.fS_EMG];
plot(x(1:end-1),Data.VICON.EMG(:,1)*10,'b')
x = [0:1/Data.Recorded_Data.fS_robot:length(Data.Recorded_Data.Fz.data)/Data.Recorded_Data.fS_robot];
plot(x(1:end-1),Data.Recorded_Data.Fz.data,'r')
plot(x(1:end-1),Data.Recorded_Data.cicles.data,'g')
for i_del = 1:length(del_buffer)
    index = find(Data.Recorded_Data.cicles.data==del_buffer(i_del));
    if ~isempty(index)
    scatter(index(1)/Data.Recorded_Data.fS_robot, del_buffer(i_del)/10,'k')
    end
end
ylim([-20,25])
hold on

prompt = {'Do you want to delete other trials? (y=1/n=0)','Number of trial to delete','Number of trial to delete',...
    'Number of trial to delete','Number of trial to delete','Number of trial to delete','Number of trial to delete',...
    'Number of trial to delete','Number of trial to delete','Number of trial to delete','Number of trial to delete',...
    'Number of trial to delete','Number of trial to delete','Number of trial to delete','Number of trial to delete'};
x = inputdlg(prompt)

if str2num(x{1})==1
    for i=2:length(x)
        del_buffer =[del_buffer;str2num(x{i})];
    end
end

%% eliminate in trials the del_buffer numbers and update Data
del_buffer = unique (del_buffer);
trials_good = [];
for j_trials = 1:length(trials)
    check = find(del_buffer==trials(j_trials));
    if isempty(check)
        trials_good = [trials_good,trials(j_trials)];
    end
end
% save also the start and the end of the good trials
beg_trials = NaN (length(trials_good),1);
fine_trials = NaN (length(trials_good),1);
for i_good = 1 : length(trials_good)
    index_single_trial = find (Data.Recorded_Data.cicles.data==trials_good(i_good) & Data.Recorded_Data.T_status.data~=0);
    beg_trials(i_good) = index_single_trial(1)/Data.Recorded_Data.fS_robot;
    fine_trials(i_good) = index_single_trial(end)/Data.Recorded_Data.fS_robot;
end
Data.good_trials = [trials_good' beg_trials fine_trials];

end