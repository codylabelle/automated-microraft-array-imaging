%Subtract image from image
imgS = medfilt2(double(NeedleImg{1}) - double(NeedleImg{2}));
imtool(imgS,[])
%Binarize image; Otsu's method
%thresh = multithresh(imgS,3);
BW = imgS>max(imgS(:)-(1/7)*max(imgS(:)));
imtool(BW)
%BW = imbinarize(imgS);
%Remove small objects
se2 = strel('disk',2);
BW2 = imopen(BW,se2);
%Fill in needle object -- if any holes exist
se = strel('disk',15);
BW3 = imclose(BW2,se);
%Dilate point to create a disk shape --> improve finding centroid
BW3 = imdilate(BW3,se);

%label needle tip
LabelNeedle = bwlabel(BW3);
%Find the centroid of the needle
NeedleCent = regionprops(LabelNeedle,'centroid');