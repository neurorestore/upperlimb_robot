%% Function to analyse peaks in  moment
function [peaks,peaks_contra,name] = CalculatePEAKS_moments (signal,fS)

Int_SearchMaxMin=40;
     
% Filter
fx  = sgolayfilt(signal,3,31);
%fx  = cheb2LPfilt(signal,30,3,fS);
        
% derivata forza 
df   = derivative(fx',1/fS);

%trova un valore medio (nel caso ci fosse un offset) e la deviazione
  %standard
       
mfx  = median(fx);
fstd = std(fx);
thd = mfx-2*fstd;
th = mfx+2*fstd;

[pks,locs] = findpeaks_GUI(-fx,'MINPEAKHEIGHT',-thd,'MINPEAKDISTANCE',30);
pks=-pks;

[pks_contra,locs_contra] = findpeaks_GUI(fx,'MINPEAKHEIGHT',th,'MINPEAKDISTANCE',30);
% elimino i picchi di forza che si hanno prima di 10 samples e 10 samples prima della fine
for i=1:length(pks)
    if locs(i) <Int_SearchMaxMin+1 || locs(i)> length(fx)-Int_SearchMaxMin-1
        pks(i)=0;
    end
end
locs(pks==0)=[];
pks(pks==0)=[];

for i=1:length(pks_contra)
    if locs_contra(i) <Int_SearchMaxMin+1 || locs_contra(i)> length(fx)-Int_SearchMaxMin-1
        pks_contra(i)=0;
    end
end
locs_contra(pks_contra==0)=[];
pks_contra(pks_contra==0)=[];
 %elimino i picchi di forza nel cui intorno la deviazione standard è troppo bassa
 % perchè non sono picchi dell'animale ma effetti del rumore
for i=1:length(pks)
    sD(i) =  std(fx(locs(i)-10 : locs(i)+10));
    if sD(i) <0.025
        pks(i)=0;
    end
end
locs(pks==0)=[];
pks(pks==0)=[];

for i=1:length(pks_contra)
    sD(i) =  std(fx(locs_contra(i)-10 : locs_contra(i)+10));
    if sD(i) <0.025
        pks_contra(i)=0;
    end
end
locs_contra(pks_contra==0)=[];
pks_contra(pks_contra==0)=[];

%% Find onset and end of every peak (pks)

ThFP = 0.5;
fwhm_a = NaN(length(pks),1); % width of the peak at half value
width_peak = NaN(length(pks),1); % distance between onset and end of the peak
AUC_F_fwhm = NaN(length(pks),1); % area under the curve calculate on the fwhm
AUC_F_min =  NaN(length(pks),1); % area under the curve calculate on the onset
onset_F =  NaN(length(pks),1); % start of the force peak
end_F =  NaN(length(pks),1); % end of the force peak
ratio_ampl = NaN(length(pks),1); % ratio between the amplitude of the peak and the amplitude of the onset
ratio_pk_fwhm = NaN(length(pks),1); % ratio between the amplitue of the peak and the width
smoothness = NaN(length(pks),1); % number of peaks in the derivative of the peaks
df_up = NaN(length(pks),1); % mean value of the derivative during the rise of peaks
df_down = NaN(length(pks),1); % mean value of the derivative during the decrease of peaks
            
for i=1:length(pks)
    %%% trova onset picco forza (a partire dal primo SubPeak) -> uso
    %%% la derivata e cerco il primo valore sotto una certa soglia
    IntervalDerF_O = df( locs(i)-Int_SearchMaxMin : locs(i));
            
    %trovo il massimo della derivata
    [Max_intDF, i_Max_intDF] = min(IntervalDerF_O);
            
    i_IntervalDerF_SubTh_O = [];
    ThFP_N = ThFP;
    while isempty(i_IntervalDerF_SubTh_O)
         [i_IntervalDerF_SubTh_O] = find( IntervalDerF_O(1:i_Max_intDF) > - ThFP_N);
         ThFP_N = ThFP_N+0.1;
    end
    
    %%onset Force Peak %%%
    onset_force_peak =  i_IntervalDerF_SubTh_O(end) + locs(i) - Int_SearchMaxMin;
            
    %%% trova END picco forza (a partire dall'ultimo SubPeak) -> uso
    %%% la derivata e cerco il primo valore sotto una certa soglia
    IntervalDerF_E = df( locs(i) : locs(i)+Int_SearchMaxMin);
            
     %trovo il minimo della derivata
     [Min_intDF, i_Min_intDF] = max(IntervalDerF_E);
            
     i_IntervalDerF_SubTh_E = [];
     ThFP_N = ThFP;
     while isempty(i_IntervalDerF_SubTh_E)
         [i_IntervalDerF_SubTh_E] = find(IntervalDerF_E(i_Min_intDF:end) < ThFP_N);
          ThFP_N = ThFP_N+0.1;
     end
      
     %%end Force Peak %%%
     end_force_peak = i_IntervalDerF_SubTh_E(1) + locs(i) + i_Min_intDF;
            
            
     start_fwhm=find(fx(1:locs(i))>=pks(i)/2);
     stop_fwhm=find(fx(locs(i):end)>=pks(i)/2)+locs(i);
     rel_pre = find((islocalmin(fx(onset_force_peak:locs(i))))==1)+onset_force_peak;
     rel_post = find((islocalmin(fx(locs(i):end_force_peak)))==1)+locs(i);
     if isempty(rel_pre)
         rel_pre = onset_force_peak;
     end
     if isempty(rel_post)
         rel_post = end_force_peak;
     end
    if isempty(stop_fwhm)
        stop_fwhm=end_force_peak;
    end
    if isempty(start_fwhm)
        start_fwhm=onset_force_peak;
    end
    if rel_pre(end)>start_fwhm(end)
        st_fwhm = rel_pre(end);
    else
        st_fwhm = start_fwhm(end);
    end
    if rel_post(1)>stop_fwhm(1)
        en_fwhm = stop_fwhm(1);
    else
        en_fwhm = rel_post(1);
    end
    fwhm_a(i)=en_fwhm-st_fwhm;
    
    % primo AUC: prende la forma della curva fino a metà altezza e poi scende diritto quindi area di un rettangolo
    AUC_F_fwhm(i) = trapz( abs (  fx( st_fwhm : en_fwhm)))+fwhm_a(i)*abs(pks(i))/2;
    % secondo AUC:area sotto la curva fino ai punti di minimo relativo precedenti e successivi e poi scende dritto fino a zero
    A_trapezio = (fx(rel_pre(end))+fx(rel_post(1)))*(rel_post(1)-rel_pre(end))/2;
    AUC_F_min(i) = trapz(abs(fx(rel_pre(end):rel_post(1))))+ A_trapezio;
    onset_F(i)  = onset_force_peak;
    end_F(i)    = end_force_peak;
    ratio_ampl(i) = abs(pks(i)/fx(onset_force_peak));
    ratio_pk_fwhm(i) = abs(pks(i)/fwhm_a(i));
    df_up(i) = mean(df(onset_force_peak:locs(i)));
    df_down(i) = mean(df(locs(i):end_force_peak));
    width_peak(i) = end_force_peak-onset_force_peak;
    
    
    % find the smoothness of the peaks (pks)
    if end_force_peak-onset_force_peak>=3
       d = findpeaks(-df(onset_force_peak:end_force_peak));
    else
       d = 1;
    end

    if ~isempty(d)
       smoothness(i) = length(d);
    else
       smoothness(i) = 1;
    end
end 
peaks = [locs',onset_F,end_F,pks',fwhm_a,width_peak,ratio_ampl,...
    ratio_pk_fwhm,AUC_F_fwhm, AUC_F_min,df_up,df_down,smoothness];

figure
hold on 
plot (signal)
plot(fx)
scatter(locs,pks)
scatter(onset_F,fx(onset_F))
scatter(end_F,fx(end_F))
hold off

% media = mean(peaks);
% media(4)=media(4)*20;
% media(8)=media(8)*100;
% dev= std(peaks);
% dev(4)=dev(4)*20;
% dev(8)=dev(8)*100;
% figure
% hold on
% bar(media(4:end))
% errorbar(media(4:end),dev(4:end),'.k')
%         

%% Find onset and end of every peak CONTRA (pks_contra)

ThFP = 0.5;
fwhm_a = NaN(length(pks_contra),1); % width of the peak at half value
width_peak = NaN(length(pks_contra),1); % distance between onset and end of the peak
AUC_F_fwhm = NaN(length(pks_contra),1); % area under the curve calculate on the fwhm
AUC_F_min =  NaN(length(pks_contra),1); % area under the curve calculate on the onset
onset_F =  NaN(length(pks_contra),1); % start of the force peak
end_F =  NaN(length(pks_contra),1); % end of the force peak
ratio_ampl = NaN(length(pks_contra),1); % ratio between the amplitude of the peak and the amplitude of the onset
ratio_pk_fwhm = NaN(length(pks_contra),1); % ratio between the amplitue of the peak and the width
smoothness = NaN(length(pks_contra),1); % number of peaks in the derivative of the peaks
df_up = NaN(length(pks_contra),1); % mean value of the derivative during the rise of peaks
df_down = NaN(length(pks_contra),1); % mean value of the derivative during the decrease of peaks
            
for i=1:length(pks_contra)
    %%% trova onset picco forza (a partire dal primo SubPeak) -> uso
    %%% la derivata e cerco il primo valore sotto una certa soglia
    IntervalDerF_O = df( locs_contra(i)-Int_SearchMaxMin : locs_contra(i));
            
    %trovo il massimo della derivata
    [Max_intDF, i_Max_intDF] = max(IntervalDerF_O);
            
    i_IntervalDerF_SubTh_O = [];
    ThFP_N = ThFP;
    while isempty(i_IntervalDerF_SubTh_O)
         [i_IntervalDerF_SubTh_O] = find( IntervalDerF_O(1:i_Max_intDF) <  ThFP_N);
         ThFP_N = ThFP_N+0.1;
    end
    
    %%onset Force Peak %%%
    onset_force_peak =  i_IntervalDerF_SubTh_O(end) + locs_contra(i) - Int_SearchMaxMin;
            
    %%% trova END picco forza (a partire dall'ultimo SubPeak) -> uso
    %%% la derivata e cerco il primo valore sotto una certa soglia
    IntervalDerF_E = df( locs_contra(i) : locs_contra(i)+Int_SearchMaxMin);
            
     %trovo il minimo della derivata
     [Min_intDF, i_Min_intDF] = min(IntervalDerF_E);
            
     i_IntervalDerF_SubTh_E = [];
     ThFP_N = ThFP;
     while isempty(i_IntervalDerF_SubTh_E)
         [i_IntervalDerF_SubTh_E] = find(IntervalDerF_E(i_Min_intDF:end) >- ThFP_N);
          ThFP_N = ThFP_N+0.1;
     end
      
     %%end Force Peak %%%
     end_force_peak = i_IntervalDerF_SubTh_E(1) + locs_contra(i) + i_Min_intDF;
            
            
     start_fwhm=find(fx(1:locs_contra(i))>=pks_contra(i)/2);
     stop_fwhm=find(fx(locs_contra(i):end)>=pks_contra(i)/2)+locs_contra(i);
     rel_pre = find((islocalmin(fx(onset_force_peak:locs_contra(i))))==1)+onset_force_peak;
     rel_post = find((islocalmin(fx(locs_contra(i):end_force_peak)))==1)+locs_contra(i);
     if isempty(rel_pre)
         rel_pre = onset_force_peak;
     end
     if isempty(rel_post)
         rel_post = end_force_peak;
     end
    if isempty(stop_fwhm)
        stop_fwhm=length(fx);
    end
    if isempty(start_fwhm)
         start_fwhm = onset_force_peak;
     end
    if rel_pre(end)>start_fwhm(end)
        st_fwhm = rel_pre(end);
    else
        st_fwhm = start_fwhm(end);
    end
    if rel_post(1)>stop_fwhm(1)
        en_fwhm = stop_fwhm(1);
    else
        en_fwhm = rel_post(1);
    end
    fwhm_a(i)=en_fwhm-st_fwhm;
    
    % primo AUC: prende la forma della curva fino a metà altezza e poi scende diritto quindi area di un rettangolo
    AUC_F_fwhm(i) = trapz( abs (  fx( st_fwhm : en_fwhm)))+fwhm_a(i)*abs(pks_contra(i))/2;
    % secondo AUC:area sotto la curva fino ai punti di minimo relativo precedenti e successivi e poi scende dritto fino a zero
    A_trapezio = (fx(rel_pre(end))+fx(rel_post(1)))*(rel_post(1)-rel_pre(end))/2;
    AUC_F_min(i) = trapz(abs(fx(rel_pre(end):rel_post(1))))+ abs(A_trapezio);
    onset_F(i)  = onset_force_peak;
    end_F(i)    = end_force_peak;
    ratio_ampl(i) = abs(pks_contra(i)/fx(onset_force_peak));
    ratio_pk_fwhm(i) = abs(pks_contra(i)/fwhm_a(i));
    df_up(i) = mean(df(onset_force_peak:locs_contra(i)));
    df_down(i) = mean(df(locs_contra(i):end_force_peak));
    width_peak(i) = end_force_peak-onset_force_peak;
    
    
    % find the smoothness of the peaks (pks)
    if end_force_peak-onset_force_peak>=3
       d = findpeaks(-df(onset_force_peak:end_force_peak));
    else
       d = 1;
    end

    if ~isempty(d)
       smoothness(i) = length(d);
    else
       smoothness(i) = 1;
    end
end 
peaks_contra = [locs_contra',onset_F,end_F,pks_contra',fwhm_a,width_peak,ratio_ampl,...
    ratio_pk_fwhm,AUC_F_fwhm, AUC_F_min,df_up,df_down,smoothness];

name = {'ind_peak','onset','endset','peak_amp','fwhm','width','ratio_amp','ratio_pk_fwhm'...
    'AUC_fwhm','AUC_min','df_up','df_down','smoothness'};

figure
hold on 
plot (signal)
plot(fx)
scatter(locs_contra,pks_contra)
scatter(onset_F,fx(onset_F))
scatter(end_F,fx(end_F))
hold off
end