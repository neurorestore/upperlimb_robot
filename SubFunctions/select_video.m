% Select video where the LED is more visible and check if video are all the
% same number of frames

function [videofile_new] = select_video(videofile, videopath)

name = videofile(1:end-5);
FigH = figure('Name','Choose the video for the led extraction','UserData',[]);

hold on
for i_v = 1:4% we have four cameras
    name_ok = [name,int2str(i_v),'.avi'];
    name_obj = VideoReader([videopath name_ok]);
    tot_frame(i_v) = name_obj.Duration*name_obj.FrameRate;
    fr = readFrame(name_obj);
    gray_fr = rgb2gray(fr);
    subplot(2,2,i_v)
    imshow(gray_fr)
end

f =gcf;
pause
axesClicked = gca;
allAxes = findobj(f.Children,'Type','axes');
numClicked = find(axesClicked==allAxes);
number = 5-numClicked;
videofile_new = [name,int2str(number),'.avi'];
close

if length(unique(tot_frame))>1
    display('Different length of the video')
    tot_frame
end
end