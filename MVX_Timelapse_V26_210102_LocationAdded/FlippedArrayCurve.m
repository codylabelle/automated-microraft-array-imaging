function [needle_Travel_Distances] = FlippedArrayCurve(x,y,z,ZCorrection)
% x, y, and z are matched vectors or matrices containing 3-dimensional
% point-cloud data of the position of the microraft arrays surface 
% x =Positions{1}(:,:,1);
% y =Positions{1}(:,:,2);
% z =Positions{1}(:,:,3);

G=[x(:),y(:),z(:),ones(length(x(:)),1)];
[~, ~, v] = svd(G, 0);
P = v(:,4);
%Equation of the xyz plane; z = ax + by + cz + d; solving for a, b, c, d
%Flip over axis
z_slant=z-(-P(1).*x-P(2).*y-P(4))./P(3);
%Correct for offset created by finding the equation of the plane
%Output is in microns
needle_Travel_Distances=(z_slant+abs(prctile(z_slant(:),5)))*ZCorrection;

