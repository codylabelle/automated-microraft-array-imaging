function [j,k,CENTROID] = ANALYZE_IMG_rescan(img,j,k,pixel_size,Width,Gap)

%segment microrafts and finds centroids 
L_bw = uraft(img{6},Width,Gap,pixel_size);
ZL = regionprops(L_bw,'Centroid');
CENTROID = cell2mat({ZL.Centroid}');