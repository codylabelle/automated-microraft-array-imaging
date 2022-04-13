function [Center_Loc2] = ConditionRafts(Center_Loc,Array_Dim)

%Distances used for range finding --> may need to be tweaked
GW = sum(Array_Dim(1:2)); % Used to add in neighbor points
GW2 = sum(Array_Dim(1:2))+Array_Dim(2)/2; %Used to look for all neighbor points
W = Array_Dim(1)/1.25; %Used for averaging raft locations and to check for existing raft locations compared to added in neighbor points

%Averaging points within 'W' distance (about 2/3 microraft)
Center_Loc2 = zeros(size(Center_Loc));
IND = rangesearch(Center_Loc,Center_Loc,W);
for ii = 1:length(IND)
    IND2 = IND{ii}';
    if length(IND2)>1
        IND2 = sort(IND2,'ascend');
        Center_Loc2(ii,:) = mean(Center_Loc(IND2,:));
    else
        Center_Loc2(ii,:) = Center_Loc(ii,:);
    end
end
Center_Loc2 = unique(Center_Loc2,'rows','stable');

%Fills in missing holes
Center_Loc3 = Center_Loc2;
k = boundary(Center_Loc3,1); %find boundary points and delete
Center_Loc3(k,:) = [];
D = rangesearch(Center_Loc2,Center_Loc3,GW2); %looking for points within 175 microns for each point
IND = find(cellfun(@(x)size(x,2),D)<5)'; %Find points with less than 5 neighbors --> includes self as neighbor

%While there are points that have less than 5 neighbors and less than 50
%attempts have been made to fill holes
fillCounter = 1;
while ~isempty(IND) && fillCounter < 25
    pts = [];
    for i = IND  %Adding in neighbor points
        pts = [pts;Center_Loc3(i,1)-GW,Center_Loc3(i,2);Center_Loc3(i,1)+GW,Center_Loc3(i,2);Center_Loc3(i,1),Center_Loc3(i,2)-GW;Center_Loc3(i,1),Center_Loc3(i,2)+GW];
    end
    %Check if new points are within 'W' distance to existing rafts
    D = rangesearch(Center_Loc2,pts,W);
    Center_Loc2 = [Center_Loc2;pts(cellfun(@(x)size(x,2),D)==0,:)]; %Add points that have no exisiting point in that location
    
    %Average new points that were added
    Center_Loc = Center_Loc2;
    Center_Loc2 = zeros(size(Center_Loc));
    IND = rangesearch(Center_Loc,Center_Loc,W);
    for ii = 1:length(IND)
        IND2 = IND{ii}';
        if length(IND2)>1
            IND2 = sort(IND2,'ascend');
            Center_Loc2(ii,:) = mean(Center_Loc(IND2,:));
        else
            Center_Loc2(ii,:) = Center_Loc(ii,:);
        end
    end
    Center_Loc2 = unique(Center_Loc2,'rows','stable');
    %Re-identify points with less than five neighbors
    Center_Loc3 = Center_Loc2;
    k = boundary(Center_Loc3,1);
    Center_Loc3(k,:) = [];
    D = rangesearch(Center_Loc2,Center_Loc3,GW2);
    IND = find(cellfun(@(x)size(x,2),D)<5)';
    fillCounter = fillCounter + 1;
end
%Remove any points that have less than two neighbors --> random points
%outside of the array
D = rangesearch(Center_Loc2,Center_Loc2,175);
IND = find(cellfun(@(x)size(x,2),D)<2)';
Center_Loc2(IND,:) = [];