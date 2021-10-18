% Function to create a plot for every file where all the units are divided
% between excitatory and inibhitory
function PlotInvsEx(Data,indexMap,cell_type,cell_type_act,in_case,fScalcium,fSrobot)

INT = 1; % 1 second of interval area

for i_file = 1:size(Data,2)
    switch in_case
        case 1 % index = force peaks
            index = Data{1,i_file}.Analysis.Fzpeaks(:,1); %index of the force peak
            if mod(i_file,2)<0.5
                phase = Data{1,i_file}.Analysis.Fzpeaks(:,15);
                index(phase<1.5)=[]; % select only index during pulling phase
                                        % because we are working in the
                                        % half tasks
            end
            index = round(index/fSrobot*fScalcium);
            
        case 2 % index = start movement
            index = Data{3,i_file}.start;
            index = round(index*fScalcium);
%             if mod(i_file,2)>0.5
%                 status = Data{1,i_file}.T_status.data;
%                 inin = find(status==2);
%                 index = inin(find (diff(inin)>1.5)+1);
%                 index = [inin(1),index];
%                 index = round(index/fSrobot*fScalcium);
%             else
%                 index = Data{3,i_file}.pks{2};
%                 index = round(index*fScalcium);
%             end              
    end
    index = index/fScalcium;
    if mod(i_file,2)>0.5 % active file
        type = cell_type_act;
        mov = resample(Data{1,i_file}.pos_Spindle_drive.data,fScalcium,fSrobot);        
    else % half file
        type = cell_type;
        mov = resample(Data{3,i_file}.x,fScalcium,Data{3,i_file}.fS_KIN);
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