global mmc
%%
mmc.waitForSystem();
img = mmc.getLastImage();
img = typecast(img,'uint16');
width = mmc.getImageWidth();
height = mmc.getImageHeight();
img = imrotate(reshape(img,[width,height]),-90);

pixel_size = 0.5247;
Width= 100;
Gap= 50;
%%

BF = img;
sigma = (Gap)/(1.5*pixel_size);
BF =  imflatfield(img,sigma);
% bw = imbinarize(BF);
% 
thresh = multithresh(BF,2);
bw = BF>(thresh(1));

imtool(bw)
%%
bw = imcomplement(bw);
%imtool(bw)
bw = imclose(bw,strel('disk',floor(0.3*Gap/pixel_size)));
bw = imfill(bw,'holes');
%Attemps to separate microrafts that were segmented together
bw = imopen(bw,strel('square',floor(0.8*Width/pixel_size)));
%%%%%%%%%%%%%%%%%
bw = imclearborder(bw);
%Remove microraft that is smaller than 2% of a nominal microraft
bw = bwareaopen(bw,round(0.02*(Width/pixel_size)^2));
%Removes any microraft larger than 150% are larger than nominal microraft
L_bw = bw&~bwareaopen(bw,round(1.5*(Width/pixel_size)^2));

imtool(L_bw)