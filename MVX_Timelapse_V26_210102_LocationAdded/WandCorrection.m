function [xModified,yModified] = WandCorrection(currentX, currentY, Raft_locations_cell)
%Define the wand radius in microns: 2.5mm =  2500 microns
wandRadius = 1500/2;  %4760 for large wand
% wandRadius = 4760/2; 
%Load in the x and y coordinates of entire array
xCords = Raft_locations_cell{1,1}(:,1);
yCords = Raft_locations_cell{1,1}(:,2);
%Find any outlier points and set x,y to zero
TFX = isoutlier(xCords);
TFY = isoutlier(yCords);
xCords(TFX) = 0;
yCords(TFX) = 0;
xCords(TFY) = 0;
yCords(TFY) = 0;
%Find the minimum and maximum location of the rafts -- don't include zero
%points
minX = min(xCords(xCords~=0));
maxX = max(xCords(xCords~=0));
minY = min(yCords(yCords~=0));
maxY = max(yCords(yCords~=0));
raftIndex = [currentX, currentY];
%Find distances between the target microraft and edges of array
deltaX1 = abs(raftIndex(:,1) - minX);
deltaX2 = abs(raftIndex(:,1) - maxX);
deltaY1 = abs(raftIndex(:,2) - minY);
deltaY2 = abs(raftIndex(:,2) - maxY);
%Shift the x coordinate if too close to x min or x max
newX = [];
newY = [];

if deltaX1 < wandRadius
    newX = raftIndex(:,1) + (wandRadius - deltaX1);
elseif deltaX2 < wandRadius
    newX = raftIndex(:,1) - (wandRadius - deltaX2);
end
if ~isempty(newX)
    raftIndex(:,1) = newX;
end

%Shift the y coordinate if too close to y min or y max
if deltaY1 < wandRadius
    newY = raftIndex(:,2) + (wandRadius - deltaY1);
elseif deltaY2 < wandRadius
    newY = raftIndex(:,2) - (wandRadius - deltaY2);
end
if ~isempty(newY)
    raftIndex(:,2) = newY;
end
%Return modified coordinates for wand movement: output from function
xModified = raftIndex(:,1);
yModified = raftIndex(:,2);
end