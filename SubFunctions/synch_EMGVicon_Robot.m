% function to extract the start point with Vicon

function [start_EMG] = synch_EMGVicon_Robot(trigger)
%moving average filter
bb = 1/25*ones(25,1);
trig_filt = filtfilt(bb,1,double(trigger));
up = find (trig_filt>1);

for i = 2:length(up)
    diff(i) = up(i)-up(i-1);
end

doubleON = find(diff>1); % error double start(or double stop) robot during recordings 

if isempty(doubleON) % no error signal up to 5 only once
    start_EMG = up(1);
else % error during recordings, look for the correct start point
    d = length(up);
    doubleON = [1, doubleON, d];
    for j = 2:length(doubleON)
        int(j-1) = doubleON(j)-doubleON(j-1);
    end
    [~,i_int] = max(int);
    start_EMG = up(doubleON(i_int));
    disp ('Double-click during recordings')
end

figure;plot(trig_filt)
hold on 
scatter (start_EMG,trig_filt(start_EMG))
end