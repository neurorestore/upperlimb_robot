%% Comparison Half along days Inscopix
% 1.upload all data
% 2.choose useful cells (units presents in at least half of the files after injury)
% 3.study a window of 2s around index for every cell to discriminate 
%   eccitatory units from inibitory from nothing
% ATTENTION the index can be force peak or movement onset
% 4.
% 5.find resting intervals
% 6.frequency rate and and raster calcium events
% 7.create a pie chart for the evolution of SUs

close all
clear all

% Select Rat
Rat = 'R04';
% DAYS tot depends on the RegCell matrix, that is how many days have been
% calculated together to find the same cells between different days
DAYStot = {'DAY01','DAY02','DAY03','DAY04','DAY05','DAY07','DAY09','DAY10','DAY11','DAY12','DAY13'}; % Rat 04
%DAYStot = {'DAY03','DAY04','DAY05','DAY08','DAY09','DAY10','DAY11','DAY12','DAY13'}; % Rat 08

% Select Index (case 1-> force peak, case 2-> start movement)
in_case = 1; %possibility to add other switches
% Select Task ('Push_active' or 'active')
Task = {'Push_active'};

% Size window
MENO = 2; %seconds
LITTLEPIU = 0.7; %seconds to calculate if it is excitatory or inhibitory cell
PIU = 1.2; %seconds window plot

% Select folder for file mat analysis
data_robot = 'C:\R-Platform\DATA\2019_02_Group1\Recordings_robot\';
% Select folder for inscopix
data_inscopix = 'C:\R-Platform\DATA\2019_02_Group1\SUs-software inscopix\';
add = 'RegCells'; %folder with the matrix output of GUI CellReg
DAYS = {'ALLdays'};
DAYSsel = {'DAY03','DAY04','DAY05','DAY07','DAY09_4DPI','DAY10_8DPI','DAY11_14DPI','DAY12_21DPI','DAY13_35DPI'};

% Sample frequency
fSrobot = 100; %robot
fSemg = 1000; %VICON
fSkin = 50; %SIMI
fScalcium = 20; %inscopix

CodePath = cd;

%% 1. LOAD DATA
% Data-> cell structure; columns different days
%        rows (1.robot, 2. EMG, 3.kin, 4.inscopix, 5.goodTrials)

cd([CodePath,'\SubFunctions\Inscopix'])
[Data] = LoadDataIns(data_robot,DAYSsel,Rat,Task);


%% 2. Choose useful cells
%cells that are present in at least half of the data after injury
[indexMap] = SelectGoodCells_vDays(data_inscopix,Rat,add,DAYS,DAYSsel,DAYStot,Task);


%% 3. Discriminate different cell type in the half tasks
% cell type: 1->mov active, 2->mov silent, 0-> do not partecipate

% Compare number of cells that are excitatory or inhibitory in different
% days
cell_type = NumberEXvsIN_vDays(Data,indexMap,MENO,PIU,in_case,LITTLEPIU,fSrobot,fScalcium,DAYSsel);

[Matrix_allCells] = DivideCellType_vDays(Data,indexMap,MENO,PIU,in_case,LITTLEPIU,fSrobot,fScalcium,fSkin,DAYSsel,cell_type);

% plot mean of all cells together and then plot mean inhibitory and mean
% excitatory
PlotInvsEx_vDays(Data,indexMap,cell_type,in_case,fScalcium,fSrobot,Task,DAYSsel);

%% 4. Graph Force-Calcium

%% 5. Find resting intervals of one seconds (intervals without force peaks)
cd([CodePath,'\SubFunctions'])
Data = FindRest(Data,fSrobot);

%% 6. Find frequency rate and rasters calcium events
% Plot firing rate in intervals of activity and resting
cd([CodePath,'\SubFunctions\Inscopix'])
%[p,f_rate] = FreqRate(indexMap,Data,cell_type,fScalcium,fSrobot);
[Ins_results.activity_in_ToT, Ins_results.RP_tot] = RasterCalciumEvents_vDays(indexMap,Data,cell_type,fScalcium,fSrobot,fSkin,in_case,DAYSsel);

%R = Calculate_correlationCells(Half,Active,indexMap);
%% 7. Create a pie chart to study evolution of single units
%if ~strcmp(Task,'active')
    [Ins_results.n_cell_type,Ins_results.n_evolution_BL,Ins_results.n_evolution_DPI] = PieChartActive_vDays(cell_type,DAYSsel)
    save([data_inscopix,Rat,add,'\Inscopix_results_pass'],'Ins_results');
%end

% checking amplitude single cells after injury (no results)
% AmplitudeSingleCells_vDays(Data,indexMap,MENO,PIU,in_case,LITTLEPIU,fSrobot,fScalcium,fSkin,DAYSsel,cell_type);

cd(CodePath)
