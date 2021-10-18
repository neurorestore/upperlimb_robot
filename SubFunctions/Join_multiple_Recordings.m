% Function to join multiple recordings (same day, same type, same animal)

function Data_tot = Join_multiple_Recordings(Data_tot,selection,pathname)
% unisci tutti i dati
% robot (modifica tempo e numero ciclo), ma anche emg e cinematica e 
% inscopix se ci sono
cd(pathname)
forcename = {'Fx','Fy','Fz'};
%trova la registrazione con tempo minore per tagliare tutti i segnali a
%quel tempo
timeEMG = size(Data_tot.VICON.EMG,1)/Data_tot.VICON.fS_EMG;
timeRobot = length(Data_tot.Recorded_Data.t.data)/Data_tot.Recorded_Data.fS_robot;
timeKIN = Data_tot.SIMI.duration;
    
[tTot,i_min] = min([timeEMG,timeRobot,timeKIN]);

% Troncamento dei segnali EMG e dei parametri del robot al tempo minore
Data_tot.Props.SingleTime = tTot;
fields = fieldnames (Data_tot.Recorded_Data);
d = floor (tTot*Data_tot.Recorded_Data.fS_robot);
Data_tot.VICON.EMG = Data_tot.VICON.EMG(1:d*Data_tot.VICON.fS_EMG/Data_tot.Recorded_Data.fS_robot,:);
for i_field =1:length(fields)
    if isstruct(Data_tot.Recorded_Data.(fields{i_field})) && ~isempty(Data_tot.Recorded_Data.(fields{i_field}))
        Data_tot.Recorded_Data.(fields{i_field}).data=Data_tot.Recorded_Data.(fields{i_field}).data(1:d);
    end
end
% change direction force for right forelimb
    if contains(Data_tot.info.Paw,'right') % check if I have to change also moment!
        Data_tot.Recorded_Data.Fz.data = -Data_tot.Recorded_Data.Fz.data;
        Data_tot.Recorded_Data.Mx.data = -Data_tot.Recorded_Data.Mx.data;
    end
    
% remove offset to force signal
for i_force = 1:size(forcename,2)
    fx = Data_tot.Recorded_Data.(forcename{i_force}).data;
    mfx  = median(fx(Data_tot.Recorded_Data.T_status.data==0));
    Data_tot.Recorded_Data.(forcename{i_force}).data = fx-mfx;
end

Data_tot.SIMI.times = Data_tot.SIMI.duration;
Data_tot.SIMI.duration = tTot;

%% Ciclo while per aggiungere un numero di file variabile
while contains(selection,'y')
    [filename, pathname] = uigetfile({'*.mat','Matlab extracted data (*.mat)'}, 'Choose Files:','MultiSelect', 'on'); % Select File
    load([pathname,filename]);
    % change direction force for right forelimb
    if contains(Data.info.Paw,'right') % check if I have to change also moment!
        Data.Recorded_Data.Fz.data = -Data.Recorded_Data.Fz.data;
        Data.Recorded_Data.Mx.data = -Data.Recorded_Data.Mx.data;
    end
    
    % remove offset to force signal
    for i_force = 1:size(forcename,2)
        fx = Data.Recorded_Data.(forcename{i_force}).data;
        mfx  = median(fx(Data.Recorded_Data.T_status.data==0));
        Data.Recorded_Data.(forcename{i_force}).data = fx-mfx;
    end
    dim = size(Data_tot.Props.name,1);
    Data_tot.Props.name(dim+1,:)= Data.Props.name;
    Data_tot.VICON.name(dim+1,:)= Data.VICON.name;
    Data_tot.VICON.trig(dim+1)= Data.VICON.trig;
    
    % trova il tempo minore del nuovo file 
    timeEMG = size(Data.VICON.EMG,1)/Data.VICON.fS_EMG;
    timeRobot = length(Data.Recorded_Data.t.data)/Data.Recorded_Data.fS_robot;
    timeKIN = Data.SIMI.duration;
    
    % troncamento dei dati EMG e Robot secondo la registrazione più corta
    [tTot,i_min] = min([timeEMG,timeRobot,timeKIN]);
    Data_tot.Props.SingleTime(dim+1) = tTot;
    
    Data_tot.SIMI.trig(dim+1) = Data.SIMI.trig;
    Data_tot.SIMI.times(dim+1) = Data.SIMI.duration;
    Data_tot.SIMI.duration = Data_tot.SIMI.duration+tTot;
    Data_tot.SIMI.Px2cm(dim+1) = Data.SIMI.Px2cm;
    
    fields = fieldnames (Data.Recorded_Data);
    d = floor (tTot*Data.Recorded_Data.fS_robot);
    Data.VICON.EMG = Data.VICON.EMG(1:d*Data.VICON.fS_EMG/Data.Recorded_Data.fS_robot,:);
    % e unione dei file
    Data_tot.VICON.EMG = [Data_tot.VICON.EMG;Data.VICON.EMG];
    for i_field =1:length(fields)
        if isstruct(Data.Recorded_Data.(fields{i_field})) && ~isempty(Data.Recorded_Data.(fields{i_field}))
            if contains(fields{i_field},'cicles') % update of the count of cicles as the sum 
                lastcicles = unique(Data_tot.Recorded_Data.(fields{i_field}).data);
                if lastcicles(end)-lastcicles(end-1)==1
                    lastc = lastcicles(end);
                else
                    lastc = lastcicles(end-1);
                    stop_index=find(Data_tot.Recorded_Data.(fields{i_field}).data==lastcicles(end));
                    Data_tot.Recorded_Data.(fields{i_field}).data(stop_index)= lastc;
                end
                Data.Recorded_Data.(fields{i_field}).data(1:2)=0; % it starts from the last one of the previous recording
                Data.Recorded_Data.(fields{i_field}).data=Data.Recorded_Data.(fields{i_field}).data(1:d);
                Data_tot.Recorded_Data.(fields{i_field}).data = [Data_tot.Recorded_Data.(fields{i_field}).data,Data.Recorded_Data.(fields{i_field}).data+lastc];
            else
                Data.Recorded_Data.(fields{i_field}).data=Data.Recorded_Data.(fields{i_field}).data(1:d);
                Data_tot.Recorded_Data.(fields{i_field}).data = [Data_tot.Recorded_Data.(fields{i_field}).data,Data.Recorded_Data.(fields{i_field}).data];
            end
        end
    end
    tempo = length(Data_tot.Recorded_Data.t.data)/Data.Recorded_Data.fS_robot;
    Data_tot.Recorded_Data.t.data = (0.00:0.01:tempo); % aggiornamento del tempo che non deve ripartire da zero ma deve proseguire come se fosse un unico file   
    
    prompt = 'Do you want to upload other files? (y/n)';
    selection = input(prompt,'s');
end
    
end