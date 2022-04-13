function [TIMER,himage] = StartLiveHist(VIEW,HISTOGRAM,app)
global mmc
%Stop any existing sequence acquisition and start a new acquisition every
%100 ms
mmc.stopSequenceAcquisition();
mmc.startContinuousSequenceAcquisition(100)
%Set the view of the figure to himage and update using CData
himage = imshow(zeros(2048,'uint16'),'Parent',VIEW);
set(VIEW,'CLimMode','auto')
set(himage,'CDataMapping','scaled')
%Create a timer that counts every 100ms
TIMER = timer;
TIMER.ExecutionMode = 'fixedSpacing';
TIMER.Period = 0.1;

%Use timer to update histogram
TIMER.TimerFcn = {@livehistogram,HISTOGRAM,himage,app};
set(HISTOGRAM,'Visible','on')
start(TIMER)


%%%%%%%% How the image is set %%%%%%%%%%%%

%A figure is created (vid), then an axes (VIEW) is created, then an image
%(himage) is displayed
%The image is updated in livehistogram by pulling images form the camera
%and setting the CData to the new image --> refreshes the CData rather than
%creating a new image with imshow each time