function [L_bw] = uraftIsolation(img,Width,Gap,pixel_size)
BF = img;
%%%% Old Thresholding method
% %Otsu's method (two level) threshold
% thresh = multithresh(BF,2);
% bw = BF>(thresh(1));

sigma = (Gap)/(1.5*pixel_size);
BF =  imflatfield(img,sigma);
bw = imbinarize(BF);

bw = imcomplement(bw);
%Removes dark debris
%Fills in all of the rafts with white (except on border of image)
bw = imclose(bw,strel('disk',floor(0.3*Gap/pixel_size)));
bw = imfill(bw,'holes');
%Attemps to separate microrafts that were segmented together
bw = imopen(bw,strel('square',floor(0.8*Width/pixel_size)));
bw = imclearborder(bw);
%Remove microraft that is smaller than 2% of a nominal microraft
bw = bwareaopen(bw,round(0.02*(Width/pixel_size)^2));
%Removes any microraft larger than 150% are larger than nominal microraft
L_bw = bw&~bwareaopen(bw,round(1.5*(Width/pixel_size)^2));
%imtool(L_bw)