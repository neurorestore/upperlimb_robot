%% Create the matrix to use CellReg to perform a registration of the neurons between 
% different days. The spatial footprints of cellular activity (ROIs) must be provided for each session separately. 
%The matrix of the spatial footprints is of size NxMxK, where N is the number of neurons,
%M is the number of pixels in the y axis and K is the number of pixels in the x axis. 
%Each entry in the matrix is equal to the corresponding pixel's value which 
%represents its contribution to the overall cell's fluorescence.

% Select folder for inscopix
data_inscopix = 'C:\R-Platform\DATA\Group1_2019_02\SUs-software inscopix\';

% Select animal
Rat = 'R08';

% Select days to analyse
DAYS = {'DAY04'};%,'DAY03','DAY04','DAY05','DAY08','DAY09_4DPI','DAY10_8DPI','DAY11_14DPI','DAY12_21DPI','DAY13_35DPI'};%'DAY10_8DPI','DAY13_35DPI'};

% Select type of task (half, active)
Task = {'half'};

%
FOV = [900,800]; %IMPROVE THIS DATA!!! Calculate from different images

%% Load data
% initialize a varibale where all the cells info are collected:
% for each day: 1. name of the unit (by inscopix software)
% 2. activity of units along time
% 3. position (x,y) of the centroid
% 4. size of the unit

allCells = cell(length(DAYS),4);

for i_gg = 1:length(DAYS)
    pathname = [data_inscopix,Rat,'\'];
    filename = [DAYS{i_gg},'_',Task{1},'.csv'];
    inscopix = readtable([pathname,filename],'ReadVariableNames',false);
    name_cells = table2cell(inscopix(1,2:end));
    status_cells = table2cell(inscopix(2,2:end));
    cells = csvread([pathname,filename],2,1);
    dim = size(cells);
    time_cells = csvread([pathname,filename],2,0,[2,0,dim(1),0]);
    clear inscopix
    filename = [DAYS{i_gg},'_',Task{1},'-props.csv'];
    pos = csvread([pathname,filename],1,5,[1,5,dim(2),6])';
    size_cells = csvread([pathname,filename],1,8,[1,8,dim(2),8])';
    for i_c = dim(2):-1:1
        if contains (status_cells{i_c},'rejected')
            name_cells(i_c)=[];
            cells(:,i_c)=[];
            pos(:,i_c) = [];
            size_cells(:,i_c) = [];
        end
    end
    allCells{i_gg,1} = name_cells;
    allCells{i_gg,2} = cells;
    allCells{i_gg,3} = pos;
    allCells{i_gg,4} = size_cells;    
end

%% Create matrix for CellReg (with homogeneous signal for every cell)

for i_gg = 1:length(DAYS)
    totCells = size(allCells{i_gg,1},2);
    matCell = NaN(totCells,FOV(1),FOV(2));
    for i_c = 1:totCells
        center = allCells{i_gg,3}(:,i_c);
        diam = allCells{i_gg,4}(i_c)/2;
        x=(1:FOV(1)).';
        y=1:FOV(2);
        binaryMap=(x-center(1)).^2+(y-center(2)).^2 <=diam^2;
        plane = zeros(FOV(1),FOV(2));
        plane(binaryMap==1) = max(allCells{i_gg,2}(:,i_c));
        matCell(i_c,:,:)= plane;
    end
    save([pathname,DAYS{i_gg},Task{1},'Matrix_CellReg.mat'],'matCell');
end
disp('A matrix for each day has been saved. To continue use the GUI CellReg to register cell along different recordings')






