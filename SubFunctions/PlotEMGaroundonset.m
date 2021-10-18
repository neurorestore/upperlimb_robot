% Plot EMG around onset

function [P_ch,Rel_max_tot,EMGmatrix_tot] = PlotEMGaroundonset (peaks,start_mov,fS_force,EMG,EMGenv,fS,Fz, pathfold,filename, in_case,type)

nameEMG = {'Biceps','Triceps','EDC','ECR','FCU','Deltoid'};
SPECT = 0;
switch in_case
    case 1 % index = force peaks
        onsetF = peaks(:,1); %index of the force peak
        phase = peaks(:,15);
        if strcmp(type,'Pull_active')
            onsetF(phase>1.5)=[]; % select only index during pushing phase
            %onsetF(peaks(:,16)==0)=[]; % select only index associated with movement
        else
            onsetF(phase<1.5)=[]; % select only index during pulling phase
        end
        onset = round(onsetF/fS_force*fS);
    case 2 % index = start movement
        onsetF = round(start_mov*fS_force);
        onset = round(start_mov*fS);
end
% onsetF = peaks(:,1);
% % onsetF = peaks(:,1)*100;
% % onset = round(peaks(:,1)*fS);
% onset = round(peaks(:,1)/fS_force*fS);
% phase = find(peaks(:,15)>1.5); % select peaks only during pulling phase
% onset = onset(phase);
% onsetF = onsetF(phase);
MENO = 0.3*fS; % seconds before onset 0.3
PIU = 0.5*fS; % seconds after onset 0.5
x = [-MENO/fS:1/fS:PIU/fS];
xFz = [-MENO/fS:1/fS_force:PIU/fS];
Path = cd;

% eliminate onset too near to the end
delet=[];
for x_time = 1:length(onset)
    if onset(x_time)+PIU>=length(EMGenv)
        delet = [delet, x_time];
    end
    if onsetF(x_time)+PIU/fS*fS_force>=length(Fz)
        delet = [delet, x_time];
    end
end
onset(unique(delet)) = [];
onsetF(unique(delet)) = [];

Rel_max_tot = NaN(length(onset),2,size(EMG,2));
EMGmatrix_tot = NaN (size(EMG,2),MENO+PIU+1);
StartEMG_min_tot = NaN (length(onset),2,size(EMG,2));

% filter of the force
Fz = sgolayfilt(Fz,3,21);

%% Analysis of the Envelope of the EMGs around the force onset

for n_EMG = 1: size(EMG,2) %cicles for the different EMG target
    EMGmatrix = NaN(length(onset),MENO+PIU+1);
    Fzmatrix = NaN(length(onset),MENO/fS*fS_force+PIU/fS*fS_force+1);
    Rel_max = NaN (length(onset),2); % save the maximum point of the activation of the EMG
    StartEMG_min = NaN (length(onset),2); % save the minimum of the EMG before every onset
%     B = 1/25*ones(25,1);
%     EMGfilt = filtfilt(B,1,EMG(:,n_EMG));
%     EMGfilt = EMG(:,n_EMG);
figure;hold on
    for n_peaks = 1 : length (onset)
        data = EMGenv(onset(n_peaks)-MENO:onset(n_peaks)+PIU,n_EMG);
        EMGmatrix(n_peaks,:)= data;
        Fzmatrix(n_peaks,:) = Fz(onsetF(n_peaks)-MENO/fS*fS_force:onsetF(n_peaks)+PIU/fS*fS_force);
        plot(EMGmatrix(n_peaks,:))
        [Rel_max(n_peaks,1),Rel_max(n_peaks,2)] = max(data);
        Rel_max(n_peaks,2) = Rel_max(n_peaks,2)-MENO;
        % find the previous minimum to see when it start to increase the
        % activity on the EMG
        %[min_EMG,locs_min] = findpeaks_GUI(-EMGenv(onset(n_peaks)-MENO:onset(n_peaks),n_EMG));
        [min_EMG,locs_min] = min(data(1:MENO));
        StartEMG_min(n_peaks,1)= -min_EMG(end);
        StartEMG_min(n_peaks,2)= locs_min(end)-MENO;
        if Rel_max(n_peaks,2)< StartEMG_min(n_peaks,2)
            Rel_max(n_peaks,:) = NaN;
            StartEMG_min(n_peaks,:) = NaN;
            EMGmatrix(n_peaks,:) = NaN;
            Fzmatrix(n_peaks,:) = NaN;
        end
    end 
    if size(EMGmatrix,1)<2
        if isempty(n_peaks)
            y = NaN(1,MENO+PIU+1);
            yFz = NaN(1,MENO/fS*fS_force+PIU/fS*fS_force+1);
        else
            y = EMGmatrix;
            yFz = Fzmatrix;
        end
        y_dev = NaN (1,length(y));
        yFz_dev = NaN (1,length(yFz));
        Rel_max = NaN;
    else
        y = nanmean(EMGmatrix);
        yFz = nanmean(Fzmatrix);
        y_dev = nanstd(EMGmatrix)/sqrt(length(onset));
        yFz_dev = nanstd(Fzmatrix)/sqrt(length(onset));
    end
    figure;hold on; set(gca, 'FontSize', 14)
    subplot(2,1,1)
    hold on
    plot(x,y,'LineWidth',2)
    X=[x, fliplr(x)];
    Y=[y + y_dev,fliplr(y - y_dev)];
    fill( X,Y,'b');
    alpha(.10)
    xlim([-MENO/fS,PIU/fS])
    ylim(get(gca, 'YLim'))
    line([0 0], get(gca, 'YLim'), 'Color', [0 0 0], 'LineStyle', '--')
    xlabel('Time (s)')
    ylabel(['Amplitude ',nameEMG{n_EMG},' envelope'])
    
    subplot(2,1,2)
    hold on
    plot(xFz,yFz,'r','LineWidth',2)
    X=[xFz, fliplr(xFz)];
    Y=[yFz + 2 * yFz_dev,fliplr(yFz - 2 * yFz_dev)];
    fill( X,Y,'r');
    alpha(.10)
    
    xlim([-MENO/fS,PIU/fS])
    ylim(get(gca, 'YLim'))
    line([0 0], get(gca, 'YLim'), 'Color', [0 0 0], 'LineStyle', '--')
    xlabel('Time (s)')
    ylabel(['Amplitude Fz'])
    
%     cd(pathfold)
%     saveas(gca,strcat(['MeanEnvelope ' , nameEMG{n_EMG},'_',filename]),'fig')
%     cd(Path)
    
    Rel_max_tot(:,:,n_EMG)=Rel_max;
    EMGmatrix_tot(n_EMG,:) = y;
    StartEMG_min_tot(:,:,n_EMG) = StartEMG_min;
end

%% plot position maximum of the EMGenv for ch 1 and ch 2 and statistic
%remove outliers
out = isoutlier (Rel_max_tot(:,2,1));
Rel_max_tot(out,2,1)= NaN;

out = isoutlier (Rel_max_tot(:,2,2));
Rel_max_tot(out,2,2)= NaN;

%out = isoutlier (Rel_max_tot(:,2,3));
%Rel_max_tot(out,2,3)= NaN;
group = [ones(length(onset),1);ones(length(onset),1)*2];
[p,tbl,stats] = kruskalwallis([Rel_max_tot(:,2,1);Rel_max_tot(:,2,2)],group);
if length(stats)>1
    c = multcompare(stats);
end
% boxplot
data = [nanmean((Rel_max_tot(:,2,1))/fS),nanmean((Rel_max_tot(:,2,2))/fS)];
xbar = [1,2];
figure;hold on; set(gca, 'FontSize', 14)
bar(xbar,data)
errhigh = [nanstd((Rel_max_tot(:,2,1))/fS)/sqrt(length(onset)),nanstd((Rel_max_tot(:,2,2))/fS)/sqrt(length(onset))];
er = errorbar(xbar,data,-errhigh,errhigh,'HandleVisibility','off');
er.Color = [0 0 0];                            
er.LineStyle = 'none';
er.LineWidth = 1;
ylabel ('Temporal distance onsetF-EMGenvelope(s)')
xticks(xbar)
xticklabels(nameEMG)
title (['p=',num2str(p)])
 
cd(pathfold)
saveas(gca,strcat(['MaxEMG-vs-OnsetF' , filename]),'fig')
%% plot position start EMGenv for ch 1 and ch 2 and statistic
% %remove outliers
% out = isoutlier (StartEMG_min_tot(:,2,1));
% StartEMG_min_tot(out,:,:)= NaN;
% 
% out = isoutlier (StartEMG_min_tot(:,2,2));
% StartEMG_min_tot(out,:,:)= NaN;
% 
% out = isoutlier (StartEMG_min_tot(:,2,3));
% StartEMG_min_tot(out,:,:)= NaN;
% group = [ones(length(onset),1);ones(length(onset),1)*2;ones(length(onset),1)*3];
% [p,tbl,stats] = kruskalwallis([StartEMG_min_tot(:,2,1);StartEMG_min_tot(:,2,2);StartEMG_min_tot(:,2,3)],group);
% c = multcompare(stats);
% % boxplot
% data = [nanmean((StartEMG_min_tot(:,2,1))/fS),nanmean((StartEMG_min_tot(:,2,2))/fS),nanmean((StartEMG_min_tot(:,2,3))/fS)];
% xbar = [1,2,3];
% figure;hold on; set(gca, 'FontSize', 14)
% bar(xbar,data)
% errhigh = [nanstd((StartEMG_min_tot(:,2,1))/fS)/sqrt(length(onset)),nanstd((StartEMG_min_tot(:,2,2))/fS)/sqrt(length(onset)),nanstd((StartEMG_min_tot(:,2,3))/fS)/sqrt(length(onset))];
% er = errorbar(xbar,data,-errhigh,errhigh,'HandleVisibility','off');
% er.Color = [0 0 0];                            
% er.LineStyle = 'none';
% er.LineWidth = 1;
% ylabel ('Temporal distance onsetF-startEMG')
% xticks(xbar)
% xticklabels(nameEMG)
% title (['p=',num2str(p)])

% cd(pathfold)
% saveas(gca,strcat(['StartEMG-vs-OnsetF' , filename]),'fig')


%% Power Spectrum 
%%Per ora NON utile!!! Non porta nessun risultato!!

% for n_EMG = 1: size(EMG,2) %cicles for the different EMG target
%     pxxmatrix = NaN(129,length(onset));
%     fmatrix = NaN(129,length(onset));
%     %B = 1/25*ones(25,1);
%     %EMGfilt = filtfilt(B,1,EMG(:,n_EMG));
%     EMGfilt = EMG(:,n_EMG);
%     for n_peaks = 1 : length (onset)
%         data = EMGfilt(onset(n_peaks)-MENO:onset(n_peaks)+PIU);
%         [pxx,f] = pwelch(data,200,100,[],fS);
%         %figure
%         %M = medfreq(pxx,f);
%         pxxmatrix(:,n_peaks)= pxx;
%         fmatrix(:,n_peaks) = f;
%         %figure
%         %plot(f,10*log10(pxx))
%         %xlabel('Frequency (Hz)')
%         %ylabel('Magnitude (dB)')
%     end  
% end
%% calculate Spectrogram around force onset

if SPECT
window = 64;
noverlap = 60;
maxCol = 500;
minCol = 30;
L = length(onset);
nfft =1024; %size of the fft
n_freq= nfft/2+1; num_bin = floor((length(EMG)-noverlap)/(window-noverlap));
n_ch = size(EMG,2);
Spectr = NaN(n_freq,num_bin,n_ch);
F = NaN(n_freq,n_ch);
T = NaN(num_bin, n_ch);
P = NaN(n_freq,num_bin,n_ch);
% Calculate the Spectrogram for every EMG for the whole signal
for nEMG=1:n_ch
    [Spectr(:,:,nEMG),F(:,nEMG),T(:,nEMG),P(:,:,nEMG)] = spectrogram(EMG(:,nEMG)',window,noverlap,nfft,fS,'yaxis');
end

nbin=floor(((MENO+PIU)/(window-noverlap)))+1;

binMENO = round(MENO/(window-noverlap));
binPIU =nbin-binMENO;
deltaF = fS/2/(nfft/2+1);
i_maxCol = round(maxCol/deltaF);
i_minCol = round(minCol/deltaF);
P_ch = NaN (n_freq,nbin+1,n_ch);

% Calculate the mean of the spectrogram around the onset of the force peaks
for nEMG = 1:2%n_ch
    P_3D=NaN(n_freq,nbin+1,L-1);
    for ii = 1:length(onset)
        i = round(onset(ii)/(window-noverlap));
        P_3D(:,:, ii) = P(:, i-binMENO:i+binPIU,nEMG);
    end
    P_ch(:,:,nEMG)=nanmean(P_3D,3);    
end

% Plot the Spectrogram of every channel
figure, set(gca, 'FontSize', 14),
for i_plot = 1 : 2%n_ch
    subplot(2,1,i_plot)
    surf([-binMENO*(window-noverlap)/fS:(window-noverlap)/fS:binPIU*(window-noverlap)/fS],F(i_minCol:i_maxCol,i_plot),10*log10(P_ch(i_minCol:i_maxCol,:,i_plot)),'edgecolor','none'); view(0,90);
    hold on, line([0 0], get(gca, 'YLim'), 'Color', [0 0 0], 'LineStyle', '--')
    xlim([-MENO/fS,PIU/fS]); ylim([minCol,maxCol]);
    colorbar
    %caxis([-80 -30])
    xlabel('Time (s)'), ylabel(['Fr (Hz) ',nameEMG{i_plot}])
end
cd(pathfold)
saveas(gca,strcat(['Spectrogram' , filename]),'fig')
else
    P_ch = NaN;
end

end