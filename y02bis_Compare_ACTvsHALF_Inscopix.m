%% Comparison Active-Half in Baseline Inscopix
% 1.upload all data
% 2.choose useful cells (units presents in more than 1/4 files)
% 3.study a window of 2s around index for every cell to discriminate 
%   eccitatory units from inibitory from nothing
% ATTENTION the index can be force peak or movement onset
% 4.graph Force (integration) vs Calcium signal
% 5.find resting intervals
% 6.frequency rate and overall fluorescence
% 7.check push vs pull activity

close all
clear all

% Select Rat
Rat = 'R04';
% Select Index (case 1-> force peak, case 2-> start movement, case 3 -> isolated and high peaks)
in_case = 1; %possibility to add other switches
% Size window
MENO = 2; %seconds
LITTLEPIU = 0.7; %seconds to calculate if it is excitatory or inhibitory cell
PIU = 1.2; %seconds window plot
INT = 0.5; %seconds for the integration for the graph between calcium activity and force
% Select folder for file mat analysis
data_robot = 'C:\R-Platform\DATA\2019_02_Group1\Recordings_robot\';
% Select folder for inscopix
data_inscopix = 'C:\R-Platform\DATA\2019_02_Group1\SUs-software inscopix\';
add = 'RegCells'; %folder with the matrix output of GUI CellReg
DAYS = {'BL'};
%DAYS = {'ALLdays'};
DAYSBL = {'DAY03','DAY04','DAY05','DAY07'};

% DAYS tot depends on the RegCell matrix, that is how many days have been
% calculated together to find the same cells between different days
DAYSBLtot = {'DAY01','DAY02','DAY03','DAY04','DAY05','DAY07'}; %R04
%DAYSBLtot = {'DAY01','DAY02','DAY03','DAY04','DAY05','DAY07','DAY09','DAY10','DAY11','DAY12','DAY13'}; % Rat 04
%DAYSBLtot = {'DAY03','DAY04','DAY05','DAY08'}; %R08

% Sample frequency
fSrobot = 100; %robot
fSemg = 1000; %VICON
fSkin = 50; %SIMI
fScalcium = 20; %inscopix

CodePath = cd;

%% 1. LOAD DATA
% Data-> cell structure; columns different days(first active, then passive)
%        rows (1.robot, 2. EMG, 3.kin, 4.inscopix, 5.goodTrials)

cd([CodePath,'\SubFunctions'])
[Data] = LoadDataBL(data_robot,DAYSBL,Rat);


%% 2. Choose useful cells
%cells that are present in at least half of the data
[indexMap] = SelectGoodCells(data_inscopix,Rat,add,DAYS,DAYSBL,DAYSBLtot);


%% 3. Discriminate different cell type in the half tasks
% cell type: 1->excitatory, 2->inhibitory, 0-> do not partecipate

% Compare number of cells that are excitatory or inhibitory in different
% days
NumberEXvsIN(Data,indexMap,MENO,PIU,in_case,LITTLEPIU,fSrobot,fScalcium);

cd([CodePath,'\SubFunctions\Inscopix'])
[Matrix_allCells,Matrix_allCellsAct,cell_type,cell_type_act] = DivideCellType(Data,indexMap,MENO,PIU,in_case,LITTLEPIU,fSrobot,fScalcium,fSkin);

cd([CodePath,'\SubFunctions'])
% plot to check the number of cells that are involved on the task
% figure
% bar([1,2,4,5],[length(find(cell_type==1))/length(cell_type),length(find(cell_type_act==1))/length(cell_type_act)...
%     length(find(cell_type==2))/length(cell_type),length(find(cell_type_act==2))/length(cell_type_act)])
% xticks([1,2,4,5])
% xticklabels({'Half Ex','Active Ex','Half In','Active In'})
% ylabel ('Percentuage of cells involved in the task')

% plot mean of all cells together and then plot mean inhibitory and mean
% excitatory
PlotInvsEx(Data,indexMap,cell_type,cell_type_act,in_case,fScalcium,fSrobot);

%% 4. Graph Force-Calcium
% I have tried:
% - all units together / units separated
% - along all signal (only good trial) / only in the intervals selected by
% Matrix_allCells
% NO INTERESTING RESULTS!!!!! 

% GraphFCall(Data, Matrix_allCells, indexMap, INT, fSrobot, fScalcium, fSemg,MENO, cell_type,in_case,2);
% GraphFCall(Data, Matrix_allCellsAct, indexMap, INT, fSrobot, fScalcium, fSemg,MENO, cell_type,in_case,1);

%% 5. Find resting intervals of one seconds (intervals without force peaks)
Data = FindRest(Data,fSrobot);

%% 6. Find frequency rate and overall fluorescence
% Plot firing rate in intervals of activity and resting

%[p,f_rate] = FreqRate(indexMap,Data,cell_type,fScalcium,fSrobot);
RasterCalciumEvents(indexMap,Data,cell_type,fScalcium,fSrobot,in_case);
cd(CodePath)
%R = Calculate_correlationCells(Half,Active,indexMap);

%% 7. Check activity neurons between pks and pks contra in half task
% when they push they could activate different neurons
% - To obtain results data from a real pushing activity are necessary, so task pull_active
% I only obtain less cell in peaks contra that are active as the one in the
% peaks, because as peaks contra the end of the real peak is selected for
% the majority of the time
% - The code can be used also to compare different peaks in different
% direction (Fx and Fy)
% NO INTERESTING RESULTS!!!!!

% cd([CodePath,'\SubFunctions\Inscopix'])
% Neural_pushVSpull(Data,indexMap,MENO,PIU,in_case,LITTLEPIU,fSrobot,fScalcium);
% 
% cd(CodePath)
