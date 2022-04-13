function livehistogram(obj,event,ax,himage,app)
global mmc vid LED

%If the external figure (image) is closed then stop the timer and sequence
%acquisition and if the live override is off then turn off the light
if ~ishandle(vid)
    app.LiveButton.Value=0;
    app.LiveButton.Text='Live';
    mmc.stopSequenceAcquisition();
    stop(obj)
    delete(timerfindall)
    clear -global TIMER
    mmc.stopSequenceAcquisition();
    % clearvars -global TIMER
    cla(ax)
    if app.LiveOverrideCheckBox.Value==0
        powervalue=num2str(0);
        fwrite(LED,['SET:LEVEL:CHANNEL1,' powervalue ';SET:LEVEL:CHANNEL2,' powervalue ';'])
        pause(0.5)
        while LED.BytesAvailable>0
            LEDMessage = fscanf(LED);
        end
    end
%If the external figure (image) is open then create the histogram using
%the image CData
else
    % pull last image from camera
    img = mmc.getLastImage();
    img = typecast(img,'uint16');
    width = mmc.getImageWidth();
    height = mmc.getImageHeight();
    img = imrotate(reshape(img,[width,height]),-90);
    %Set image (himage) to the CData of the last image taken from camera
    %This displays on the figure (vid) and axes (VIEW)
    set(himage,'CData',img);
    [counts,~] = imhist(img, 2^8);
    bar(ax,counts);
    xlim(ax,[0 255])
    set(ax,'yscale','log')
end


