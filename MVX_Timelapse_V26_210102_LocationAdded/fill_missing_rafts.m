function [Raft_locations_fixed] = fill_missing_rafts(Raft_locations,Array_Dim)
%Use this code to estimate locations of rafts that were not segmented
% Raft Matrix = RM
RM_temp = Raft_locations;
%Find values that have coordinates of (0,0) and set to NaN for fillmissing
RM_temp((RM_temp(:,1)==0),:) = NaN;
%Reshape X and Y data into grid the size of the array
NN = Array_Dim(3);
RMX1 = reshape(RM_temp(:,1),NN,NN);
RMY1 = reshape(RM_temp(:,2),NN,NN);
%Fill in missing raft locations using linear interpolation
RMX_fixed =  fillmissing(RMX1,'linear',2,'EndValues','extrap');
RMY_fixed =  fillmissing(RMY1,'linear',1,'EndValues','extrap');
%Reshape the matrix data into a column array
RMX = reshape(RMX_fixed, [NN^2,1]);
RMY = reshape(RMY_fixed,[NN^2,1]);
%Add the X and Y column data into a single array
RM_fixed(:,1) = RMX;
RM_fixed(:,2) = RMY;
%Rename for data output
Raft_locations_fixed = RM_fixed;

%scatter(Raft_locations_fixed(:,1),Raft_locations_fixed(:,2))
%scatter(Raft_locations(:,1),Raft_locations(:,2))
end
