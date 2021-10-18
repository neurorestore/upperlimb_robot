%% EMG cose da fare

% Find the burst:
% Calculate amplitude (mean/maximum) duration and starting time
% 
% 1. envelope della curva= rettifica+low pass filter
% 2. find threshold
% 3. define every burst
%
% Spectral analysis?? Power spectrum of EMG and median power frequency

function [EMGenvtot, burst] = CalculateEMGparam(EMG, fS)
% I am not sure that these parameters are the best for the data! Check very
% carefully if you find better analysis

burst = cell(size(EMG,2),1);
EMGenvtot = NaN(size(EMG,1),size(EMG,2));
for nEMG =1:size(EMG,2)
    % calculate envelope of the curve
%     N_window = fS*0.002; %0.008
%     bb = ones(1,N_window)/N_window;
%     %EMGfilt =filtfilt(bb,1,EMG(:,nEMG));
%     EMGfilt = sgolayfilt(EMG(:,nEMG),3,13);
%     %EMGfilt = EMG(:,nEMG);
%     EMGenv = envelope(EMGfilt,400,'peak');
    
    %% Calculate envelope of the curve v.1.0
    EMGrect = abs(EMG(:,nEMG));
    B = 1/4*ones(4,1); % 4-point moving average filter, i.e., a 250-Hz low-pass filter
    EMGfilt = filtfilt(B,1,EMGrect);
    % to correct files where VICON starts later and EMG is 0 at the
    % beginning
    zerov = find(EMGfilt == 0);
    if ~isempty(zerov) && zerov(1) ==1
        d = length(zerov)+100;
        EMGenv = envelope(EMGfilt(d:end),300,'peak');
        add = zeros(d-1,1);
        EMGenv = [add;EMGenv];
    else
        EMGenv = envelope(EMGfilt,300,'peak');
    end
    EMGenvtot(:,nEMG)=EMGenv;
    
    %%
    figure;plot(EMGrect);%(EMG(:,nEMG))
    hold on
    plot(EMGfilt)
    plot(EMGenv)
    legend ('EMG','EMGfilt','EMGenv')
    
    % find threshold for burst
    EMGmed = mean(EMGfilt);
    EMGsd = std(EMGfilt);
    
    th = EMGmed+2*EMGsd;
    
    % detect start and end of every burst
    upth = find(EMGenv>th);
    index = find (diff(upth)>1);
    st_burst = [upth(1); upth(index+1)];
    en_burst = [upth(index); upth(end)];
    
    figure;plot(EMG(:,nEMG))
    hold on
    plot(EMGfilt)
    plot(EMGenv)
    plot([1,length(EMGenv)],[th,th],'r')
    scatter ( st_burst, ones(1,length(st_burst))*th)
    scatter ( en_burst, ones(1,length(en_burst))*th)
    
    %% Calculate mean/max/duration/start of the burst 
    % I don't want to normalize on the maximum here because it is
    % interesting to verify the reduction of the EMG after injury compared
    % to the healthy, maybe to evaluate all together we have to normalize
    % because the amplitude depends a lot on the implant
    
    nburst = length(st_burst);
    burst_mean = NaN(1,nburst);
    burst_max = NaN(1,nburst);
    burst_duration = NaN(1,nburst);
    burst_AUC = NaN(1,nburst);
    
    for i_burst = 1: nburst
        burst_mean(i_burst) = mean(EMGenv(st_burst(i_burst):en_burst(i_burst)));
        burst_max(i_burst) = max(EMGenv(st_burst(i_burst):en_burst(i_burst)));
        burst_duration(i_burst) = (en_burst(i_burst)-st_burst(i_burst))/fS;
        burst_AUC(i_burst) = trapz(EMGenv(st_burst(i_burst):en_burst(i_burst)))/(burst_duration(i_burst));
    end
    
    burst{nEMG} = [st_burst/fS,burst_mean',burst_max',burst_duration',burst_AUC'];
    
%     %% Frequency analysis (Spectrum-> variation of the median frequency during healthy and injured rats)
%     d=1;
%     for n = 1:20%nburst
%         if burst_duration(n)>0.5
%             figure
%             plot(EMG(st_burst(n):en_burst(n),nEMG))
%             [pxx,f] = pwelch(EMG(st_burst(n):en_burst(n),nEMG),100,50,[],fS);
%             %plot(f,10*log10(pxx))
%             %xlabel('Frequency (Hz)')
%             %ylabel('Magnitude (dB)')
%             figure
%             medfreq(pxx,f);
%             p(d) = medfreq(pxx,f);
%             d = d+1;
%         end
%     end
%    mean(p)
    
end 
end