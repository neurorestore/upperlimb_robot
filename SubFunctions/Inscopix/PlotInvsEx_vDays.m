% Function to create a plot for every file where all the units are divided
% between excitatory and inibhitory
function PlotInvsEx_vDays(Data,indexMap,cell_type,in_case,fScalcium,fSrobot,Task,DAYS)

INT = 1; % 1 second of interval area

col_DPI = contains(DAYS,'DPI');
ind_DPI = find(col_DPI==1);
n_days = sum(col_DPI)+1; %total number of useful days (that is all post injuries days and one day of baseline)
n_BL = length(find(col_DPI==0)); %number of days of baseline

%reconstruct cell_type for every file
if n_BL>1
    for i_type = 2: n_BL
        cell_type = [cell_type(:,1),cell_type];
    end
end

for i_file = 1:size(Data,2)
    switch in_case
        case 1 % index = force peaks
            index = Data{1,i_file}.Analysis.Fzpeaks(:,1); %index of the force peak
            if strcmp(Task,'Push_active')
                phase = Data{1,i_file}.Analysis.Fzpeaks(:,15);
                index(phase<1.5)=[]; % select only index during pulling phase
                                        % because we are working in the
                                        % half tasks
            end
            index = round(index/fSrobot*fScalcium);
            
        case 2 % index = start movement
            index = Data{3,i_file}.start;
            index = round(index*fScalcium);   
    end
    index = index/fScalcium;
    type = cell_type(:,i_file);
    if strcmp(Task,'Push_active')
        mov = resample(Data{3,i_file}.x,fScalcium,Data{3,i_file}.fS_KIN);               
    else % half file
        mov = resample(Data{1,i_file}.pos_Spindle_drive.data,fScalcium,fSrobot); 
    end
    
    figure; hold on
    set(gca,'FontSize',14)
    for i_ind = 1:length(index)
        xx=[index(i_ind) index(i_ind)+INT index(i_ind)+INT index(i_ind)];
        yy=[-12 -12 3 3];
        p = patch( xx,yy,[76,189,237]/255,'EdgeColor','none');
        set(p,'FaceAlpha',0.3)
    end
    datamedia = zscore(Data{4,i_file}.cells_signal);
    X= [1/fScalcium:1/fScalcium:size(datamedia,1)/fScalcium];
    plot(X,mean(datamedia,2),'Color',[125,46,141]/255,'LineWidth',1)
    inin = indexMap(type==1,i_file);
    inin(inin==0)=[];
    datamedia = zscore(Data{4,i_file}.cells_signal(:,inin));
    plot(X,mean(datamedia,2)-2,'Color',[118,171,47]/255,'LineWidth',1)
    inin = indexMap(type==2,i_file);
    inin(inin==0)=[];
    datamedia = zscore(Data{4,i_file}.cells_signal(:,inin));
    plot(X,mean(datamedia,2)-4,'Color',[236,176,31]/255,'LineWidth',1)
    
    X = [1/fScalcium:1/fScalcium:length(mov)/fScalcium];
    plot(X,(mov-1138)/100-7,'k','LineWidth',1)
    force = resample(Data{1,i_file}.Fz.data,fScalcium,fSrobot);
    X = [1/fScalcium:1/fScalcium:size(force,2)/fScalcium];
    plot(X,force-9,'k','LineWidth',1)
    xlabel ('Time (s)')

end

end