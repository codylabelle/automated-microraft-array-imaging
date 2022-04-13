function [img,camera_error_count] = Snap2(mmc)
%Take a picture with camera
mmc.waitForSystem()
camera_error_count = 0;
for counterSnap = 1:500
    SnapFlag = 1;
    try
        mmc.snapImage();
    catch
        SnapFlag = 0;
        camera_error_count = camera_error_count + 1;
        disp(camera_error_count)
    end
    if SnapFlag == 1
        break
    end
end
%if no errors occurred, set camera_error to an empty matrix
camera_error_count(find(camera_error_count == 0)) = [];

img = mmc.getImage();
img = typecast(img,'uint16');
width = mmc.getImageWidth();
height = mmc.getImageHeight();
img = imrotate(reshape(img,[width,height]),-90);

%%%  Previous code
% % mmc.waitForSystem()
% % mmc.snapImage();
% % img = mmc.getImage();
% % img = typecast(img,'uint16');
% % width = mmc.getImageWidth();
% % height = mmc.getImageHeight();
% % img = imrotate(reshape(img,[width,height]),-90);