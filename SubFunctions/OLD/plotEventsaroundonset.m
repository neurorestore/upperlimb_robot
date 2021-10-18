% Plot Neurons'events around the onset of the force peaks in zeta direction

function plotEventsaroundonset(cells,Data,day,indexMap,fS,nameday,datapath)
%%TO BE CONTINUED! FIND THE EVENTS AND CALCULATE THE HISTOGRAM
path = cd;
fS_force = Data.Recorded_Data.fS_robot;
peaks = Data.Recorded_Data.Analysis.Fzpeaks;
onset = round(peaks(:,1)/fS_force*fS);
phase = find(peaks(:,15)>1.5);
onset = onset(phase);
MENO = 0.5*fS; % seconds before onset
PIU = 1*fS; % seconds after onset
x = [-MENO/fS:1/fS:PIU/fS];

% eliminate onset too near to the end
delet =[];
for x_time = 1:length(onset)
    if onset(x_time)+PIU>length(cells{5})
        delet = x_time;
    end
end
onset(delet) = [];

%% Find Neurons common to all the selected recordings
for i_gg = 1: size(indexMap,2)
    indexMap(indexMap(:,i_gg)==0,:)=[];
end

%% Normalize each cell on its maximum
cells_raw = cells{2};
for i_r = 1:size(cells{2},2)
    cells{2}(:,i_r)= cells{2}(:,i_r)/max(cells{2}(:,i_r));
end

%% Fix the threshold to 2 std
th = NaN(size(cells{2},2),1);
activeMap = zeros(size(cells{2},1),size(cells{2},1));
for i_n = 1:length(th)
    th = mean(cells{2}(:,i_n))+2*std(cells{2}(:,i_n));
    ind_up = find(cells{2}(:,i_n)>th);
    activeMap(ind_up,i_n) = 1;
end
    
%% Analysis of the events around the force onset

NeuronMatrix_tot = NaN (size(indexMap,1),MENO+PIU+1);

for n_neuro = 1: size(indexMap,1) %cicles for the different selected neurons
    Nmatrix = NaN(length(onset),MENO+PIU+1);
    valueC = indexMap(n_neuro,day);
    for n_peaks = 1 : length (onset)
        data = activeMap(onset(n_peaks)-MENO:onset(n_peaks)+PIU,valueC);
        Nmatrix(n_peaks,:)= data;
    end 
     y = sum(Nmatrix)/length(onset);
%     y_dev = std(Nmatrix)/sqrt(length(onset));
%     figure;hold on; set(gca, 'FontSize', 14)
%     plot(x,y,'LineWidth',2)
%     X=[x, fliplr(x)];
%     Y=[y + 2 * y_dev,fliplr(y - 2 * y_dev)];
%     fill( X,Y,'b');
%     alpha(.10)
%     xlim([-MENO/fS,PIU/fS])
%     ylim(get(gca, 'YLim'))
%     line([0 0], get(gca, 'YLim'), 'Color', [0 0 0], 'LineStyle', '--')
%     xlabel('Time (s)')
%     ylabel(['Fluorescence neuron ',num2str(n_neuro)])
    
    %cd(pathfold)
    %saveas(gca,strcat(['MeanEnvelope ' , nameEMG{n_EMG},'_',filename]),'fig')
    NeuronMatrix_tot(n_neuro,:) = y;
end

%% Plot of the sum of all neurons 
y = sum(NeuronMatrix_tot,1);
figure;hold on; set(gca, 'FontSize', 14)
plot(x,y)
ylim(get(gca, 'YLim'))
line([0 0], get(gca, 'YLim'), 'Color', [0 0 0], 'LineStyle', '--')
xlabel('Time (s)')
ylabel('Events cells')
xlim([-MENO/fS,PIU/fS])
if contains(nameday,'DPI')
    titolo = nameday(7:end);
else
    titolo = 'Baseline';
end
title (titolo);
cd(datapath)
saveas(gca,strcat(['Events Cells ' , titolo]),'fig')


%% Plot surface
figure
imagesc(x,[1:1:size(NeuronMatrix_tot,1)],NeuronMatrix_tot);
colormap('jet')
line([0 0], get(gca, 'YLim'), 'Color', [0 0 0], 'LineStyle', '--')
xlabel('Time (s)')
ylabel('Neuron events(#)')
title (titolo);
saveas(gca,strcat(['Events all neurons ', titolo]),'fig')

% figure
% hold on
% for i_n = 1:size(NeuronMatrix_tot,1)
%     plot3(x,ones(length(x))*i_n,NeuronMatrix_tot(i_n,:))
%     grid
% end
% xlabel('Time (s)')
% ylabel('Neuron (#)')
% title (titolo);
cd(path)
end