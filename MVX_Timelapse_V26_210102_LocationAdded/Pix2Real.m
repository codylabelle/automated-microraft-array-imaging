function [x_real,y_real] = Pix2Real(x_pix,y_pix,pixel_size,Center_X,Center_Y,img_size)

% if size(pixel_size,2)==1
%     pixel_size = [pixel_size,pixel_size];
% end
% 
% % pix_move = [x_pix,y_pix]-[fliplr(img_size)/2+0.5];
% pix_move = [x_pix,y_pix]-[fliplr(img_size)/2+0.5];
% 
% real_move = -pix_move.*pixel_size;
% 
% Real_pos = real_move+[Center_X,Center_Y];
% x_real = Real_pos(1);
% y_real = Real_pos(2);

pix_move = [fliplr(img_size)/2+0.5]-[x_pix,y_pix];

x_real = (-pix_move(1)*pixel_size)+Center_X;
y_real = (-pix_move(2)*pixel_size)+Center_Y;