function fm = FScoreIsolation(Z_Pos)
global mmc Div

mmc.setPosition('ZStage',Z_Pos);
mmc.waitForSystem()
%Uses Snap2 function to take an image with camera
% % img = Snap2(mmc);

for counterSnap = 1:500
    SnapFlag = 1;
    try
        img = mmc.getLastImage();
    catch
        SnapFlag = 0;
    end
    if SnapFlag == 1
        break
    end
end

img = typecast(img,'uint16');
width = mmc.getImageWidth();
height = mmc.getImageHeight();
img = imrotate(reshape(img,[width,height]),-90);

%BB = bounding box; if brightfield makes boundaries of 3x3 box
fm_temp = zeros(1,Div^2);
count = 0;
for ii = 0:(Div-1)
    for jj = 0:(Div-1)
        count = count+1;
        %Takes standard deviation of image (in vector form)
        %could replace fmeasure with std2
        fm_temp(count) = fmeasure(img,'GLVA',[100+(ii*(1848/Div)),100+(jj*(1848/Div)),(1848/Div),(1848/Div)]);
    end
end
%Finds the negative median focus score from the different 3x3 microraft areas
fm = -median(fm_temp);
