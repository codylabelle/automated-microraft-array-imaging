function [] = StopLiveHist(mTimer,~,HISTOGRAM,himage)
global mmc TIMER

% delete(timerfindall)
mmc.stopSequenceAcquisition();
% clearvars -global TIMER
cla(HISTOGRAM)
set(HISTOGRAM,'Visible','off')
set(himage,'CData',zeros(2048,'uint16'))
% imshow(zeros(2048,'uint8')),'Parent',himage)
% delete(timerfindall)