%% Find Led in the video of SIMI to synchronize the kinematic data %%
function [start_LED_s,time,Pix2cm]= Extraction_LED_called

% Select file video
[videofile_new videopath] = ...
        uigetfile({'*.avi','Video Files (*.avi)'; ...
        '*.*', 'All Files (*.*)'}, 'Select the Video File');

%[videofile_new] = select_video(videofile, videopath);

% Elaborate the chosen video

obj = VideoReader([videopath videofile_new]);
obj.CurrentTime=600/obj.FrameRate;
frame = readFrame(obj);
grayImm = rgb2gray(frame); 
filtImm = medfilt2(grayImm,[6 6]);  

% Select the lenght of 2 cm

figure('Name','Select a lenght of 2 cm: two ginput')
imshow(frame)
[Lx,Ly] = ginput(2);
Pix2cm = sqrt(diff(Lx)*diff(Lx)+diff(Ly)*diff(Ly))/2;
close



% Select the rectangular area of the LED
rect_ROI_default_pos_LED =   [  300   305];
ROI_size_LED = [  46   48];
rect_LED = [rect_ROI_default_pos_LED ROI_size_LED];

figure ('Name','Select the area of the LED','NumberTitle','off')
imshow(frame);
h_rect = imrect(gca, [rect_ROI_default_pos_LED(1) rect_ROI_default_pos_LED(2) ROI_size_LED(1) ROI_size_LED(2)]); %-> finestra di 150x150 o più grande px
pause
rect_LED             = round(getPosition(h_rect));

Grgb_LED = frame(rect_LED(2):rect_LED(2)+rect_LED(4), rect_LED(1):rect_LED(1)+rect_LED(3),:);
close

% Select the pixel of the LED

figure('Name','Indicate the LED position','NumberTitle','off')
imshow(Grgb_LED);
[PxLED, PyLED] = ginput(1)
PLED = [PxLED(1,1) PyLED(1,1)];
close

% Scroll all frames and calculate the LED value

totalframes = obj.Duration*obj.FrameRate;

for i=1:600 %totalframes-2
   
    display(strcat(num2str((i/totalframes)*100),'_% completed'))
    framenumber = i;       
    obj.CurrentTime=framenumber/obj.FrameRate;
    video = readFrame(obj);
    
    grayImmLED = rgb2gray(video(...
                    rect_LED(2):rect_LED(2)+rect_LED(4),...
                    rect_LED(1):rect_LED(1)+rect_LED(3),...
                    :));            
    filtImmLED = medfilt2(grayImmLED,[6 6]);
    
    LEDgrayscale(i) = filtImmLED(round(PLED(2)),round(PLED(1)));
end

%filtering and find the threshold
LEDgrayscale = double(LEDgrayscale);
LEDgrayscale_filt = sgolayfilt(LEDgrayscale,3,41);

d = length(LEDgrayscale_filt);
der = derivative(LEDgrayscale_filt',0.02);
[peaks, locs] = findpeaks(der);
[B,I] = sort(peaks);
variance = B(end)-B(1);
n_OnOff = length(find(B>B(end)-variance/4));

if n_OnOff ==1
    start_LED = locs(I(end));
else
    display('Double-click during recordings')
    I_up = locs(I(end-n_OnOff+1:end));
    doubleON = [1, sort(I_up'), d];
    for j = 2:length(doubleON)
        int(j-1) = doubleON(j)-doubleON(j-1);
    end
    [~,i_int] = max(int);
    start_LED = doubleON(i_int);
end  

figure
plot (LEDgrayscale)
hold on
plot (LEDgrayscale_filt)
plot(der)
scatter(start_LED,der(start_LED))

%%Qui si può aggiungere una scelta a mano nel caso in cui la
%%sincronizzazione sia sbagliata
f = uifigure;
selection = uiconfirm(f,'Is start_LED correct?','Confirm Close');

if ~contains(selection,'OK')
    figure ('Name','Select the point of the start_LED','NumberTitle','off')
    plot(der)
    [Px, Py] = ginput;
    dist = locs-Px;
    dist = abs(dist);
    [~,i_locs] = min(dist);
    start_LED = locs(i_locs);
end
start_LED_s = start_LED/obj.FrameRate;
time = obj.Duration-start_LED_s;

end
