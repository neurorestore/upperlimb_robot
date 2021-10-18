% Plot Neurons around the onset of the force peaks in zeta direction

function [Rel_max_tot] = plotNeuronsaroundonset(cells,Data,day,indexMap,fS,nameday,datapath)

path = cd;

% if 1 remove outliers peaks
OUT = 1;

fS_force = Data.Recorded_Data.fS_robot;
peaks = Data.Recorded_Data.Analysis.Fzpeaks;
onset = round(peaks(:,1)/fS_force*fS);
amp = peaks(:,4);
phase = find(peaks(:,15)>1.5);
onset = onset(phase);
amp = amp(phase);
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

% eliminate onset too near to each other -> dynamic of calcium is low


%% Find Neurons common to all the selected recordings
for i_gg = 1: size(indexMap,2)
    indexMap(indexMap(:,i_gg)==0,:)=[];
end

%% Normalize on the maximum activation of that cell
cells_raw = cells{2};
for i_r = 1:size(cells{2},2)
    cells{2}(:,i_r)= cells{2}(:,i_r)/max(cells{2}(:,i_r));
end

%% Analysis of the fluorescence around the force onset

NeuronMatrix_tot = NaN (size(indexMap,1),MENO+PIU+1);
Rel_max_tot = NaN(length(onset),size(indexMap,1));
figure;hold on; set(gca, 'FontSize', 14)

for n_neuro = 1: size(indexMap,1) %cicles for the different selected neurons
    Nmatrix = NaN(length(onset),MENO+PIU+1);
    Rel_max = NaN (length(onset),1);
    valueC = indexMap(n_neuro,day);
    for n_peaks = 1 : length (onset)
        data = cells{2}(onset(n_peaks)-MENO:onset(n_peaks)+PIU,valueC);
        Nmatrix(n_peaks,:)= data;
        %figure;plot(Nmatrix(n_peaks,:))
        [~,Rel_max(n_peaks)] = max(data);
        Rel_max(n_peaks) = (Rel_max(n_peaks)-MENO)/fS;
    end
    if OUT % try to remove peaks with the same activation of fluorescence
        in_out = isoutlier(Rel_max);
        Rel_max(in_out) = NaN;
        Nmatrix(in_out,:) = NaN;
    end
     y = nanmean(Nmatrix);
    y_dev = nanstd(Nmatrix)/sqrt(length(onset));
    subplot(5,7,n_neuro)
    hold on
    plot(x,y,'LineWidth',2)
    X=[x, fliplr(x)];
    Y=[y + 2 * y_dev,fliplr(y - 2 * y_dev)];
    fill( X,Y,'b');
    alpha(.10)
    xlim([-MENO/fS,PIU/fS])
    ylim(get(gca, 'YLim'))
    line([0 0], get(gca, 'YLim'), 'Color', [0 0 0], 'LineStyle', '--')
    xlabel('Time (s)')
    ylabel(['N° ',num2str(n_neuro)])    
    NeuronMatrix_tot(n_neuro,:) = y;
    Rel_max_tot(:,n_neuro) = Rel_max;
end
cd(datapath)
saveas(gca,strcat(['Mean Fluorescence single Neurons' ,'_',nameday]),'fig')

%% Plot of the mean of all neurons 
y = mean(NeuronMatrix_tot,1);
y_dev = std(NeuronMatrix_tot,1)/sqrt(size(indexMap,1));
figure;hold on; set(gca, 'FontSize', 14)
plot(x,y,'LineWidth',2)
X=[x, fliplr(x)];
Y=[y + 2 * y_dev,fliplr(y - 2 * y_dev)];
fill( X,Y,'b');
alpha(.10)
xlim([-MENO/fS,PIU/fS])
ylim(get(gca, 'YLim'))
line([0 0], get(gca, 'YLim'), 'Color', [0 0 0], 'LineStyle', '--')
xlabel('Time (s)')
ylabel('Mean of the Fluorescence')
if contains(nameday,'DPI')
    titolo = nameday(7:end);
else
    titolo = 'Baseline';
end
title (titolo);
cd(datapath)
saveas(gca,strcat(['MeanFluorescence ' , titolo]),'fig')

%% Plot of the most active neuron

[M,iM]= max(sum(NeuronMatrix_tot,2));
[m,im]=max(sum(cells{2},1)); %select the same neuron in the files where I
%have tried to use both
y = NeuronMatrix_tot(iM,:);
Nmatrix = NaN(length(onset),MENO+PIU+1);
valueC = indexMap(iM,day);

for n_peaks = 1 : length (onset)
    data = cells{2}(onset(n_peaks)-MENO:onset(n_peaks)+PIU,valueC);
    Nmatrix(n_peaks,:)= data;
end 
if OUT % try to remove peaks with the same activation of fluorescence
    in_out = isnan(Rel_max_tot(:,iM));
    Nmatrix(in_out,:) = NaN;
end
y_dev = nanstd(Nmatrix)/sqrt(length(onset));
figure;hold on; set(gca, 'FontSize', 14)
plot(x,y,'LineWidth',2)
X=[x, fliplr(x)];
Y=[y + 2 * y_dev,fliplr(y - 2 * y_dev)];
fill( X,Y,'b');
alpha(.10)
xlim([-MENO/fS,PIU/fS])
ylim(get(gca, 'YLim'))
line([0 0], get(gca, 'YLim'), 'Color', [0 0 0], 'LineStyle', '--')
xlabel('Time (s)')
ylabel(['Fluorescence of the most active neuron',num2str(iM)])
title (titolo);
saveas(gca,strcat(['Fluorescence of the most active neuron',num2str(iM), titolo]),'fig')


%% Plot surface
figure
imagesc(x,[1:1:size(NeuronMatrix_tot,1)],NeuronMatrix_tot);
colormap('jet')
colorbar
line([0 0], get(gca, 'YLim'), 'Color', [0 0 0], 'LineStyle', '--')
xlabel('Time (s)')
ylabel('Neuron (#)')
title (titolo);
saveas(gca,strcat(['Fluorescence all neurons ', titolo]),'fig')

% figure
% hold on
% for i_n = 1:size(NeuronMatrix_tot,1)
%     plot3(x,ones(length(x))*i_n,NeuronMatrix_tot(i_n,:))
%     grid
% end

cd(path)
end