function [L_bw] = uraft(img,Width,Gap,pixel_size)
%segments individual microrafts in a brightfield image

%flat field correction  --> makes the original image have uniform
%brightness
sigma = (Gap)/(1.5*pixel_size);
BF =  imflatfield(img,sigma);

%New thresholding method
bw = imbinarize(BF);
%invert the logical image --> Borders white (and numbers) everything else
%black
bw = imcomplement(bw);

%Removes dark debris
%Fills in all of the rafts with white (except on border of image)
bw = imclose(bw,strel('disk',floor(0.25*Gap/pixel_size)));
bw = imfill(bw,'holes');
%Attemps to separate microrafts that were segmented together
bw = imopen(bw,strel('square',floor(0.5*Width/pixel_size)));

%%%%%%%%%%%%%%%%%% Try to remove rafts near border
bw = imdilate(bw,strel('square',floor((Width/pixel_size)/25)));
%%%%%%%%%%%%%%%%%
bw = imclearborder(bw);
bw = imerode(bw,strel('square',4));

% %Remove microraft that is smaller than 60% of a nominal microraft
% %Removes any microraft larger than 150% are larger than nominal microraft
%Filter out round objects/debris
ZL2 = regionprops(bw,'MajorAxisLength','Eccentricity','PixelIdxList','MinorAxisLength','Area');
removeObjects = [];
keepObjects = false(2048,2048);
for ii = 1:length(ZL2)
    if ((ZL2(ii).MajorAxisLength/ZL2(ii).MinorAxisLength) > 1.35) && ZL2(ii).Eccentricity > 0.5 || ...
            ZL2(ii).Area < (0.6*(Width/pixel_size)^2) || ZL2(ii).Area > (1.5*(Width/pixel_size)^2)
        removeObjects = [removeObjects; ii];
    else
        keepObjects(ZL2(ii).PixelIdxList) = 1;
    end
end
L_bw = keepObjects;