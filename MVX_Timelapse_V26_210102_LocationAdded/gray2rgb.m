%Pass in the image, the color channel you want (red, green, or blue), and
%the contrast adjustment values
function [rgbImage] = gray2rgb(Image,channel,low,high)
Image = imadjust(Image,[low high],[]);
[m n] = size(Image);
Gray = Image;
X = zeros(m,n);
if channel == 1
    RGB1 = cat(3, Gray, X, X);
elseif channel == 2
    RGB1 = cat(3, X, Gray, X);
else
    RGB1 = cat(3, X, X, Gray);
end

rgbImage = RGB1;