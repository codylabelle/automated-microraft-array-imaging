function [Loc,Raft_locations,Raft_locations_real] = Label_uraft(m,n,pixel_size,POS,Centroid_Cell,q,Array_Dim,DIR)
%Convert the raft locations to real world coordinates
%Creates a matrix containing the array number, imaging location row and
%column, and raft number
W = Array_Dim(1)/2;   % Y or X range 
waitb = waitbar(0);
waitbar(1/4,waitb,['Determining Raft Locations...'])
Center_Loc = [];
Loc = [];
count = 0;
for ii = 1:m
    for j = 1:n
        count = count+1;
        %Centroid_cell pixel value centers of each raft in each image
        Cent = Centroid_Cell{ii,j};
        for k = 1:size(Cent,1)
            [x_real,y_real] = Pix2Real(Cent(k,1),Cent(k,2),pixel_size,POS(ii,j,1),POS(ii,j,2),[2048 2048]);
            Center_Loc = [Center_Loc;x_real y_real];
            Loc = [Loc;q ii j k];
        end
    end
end

if exist([DIR,'\data'],'dir')
    save([DIR,'\data\label_uraft_troubleshoot.mat'],'Center_Loc','POS','Centroid_Cell','q','Array_Dim')
end

%Averages all points --> fills in gaps, removes points outside of main
%array
disp('Starting raft conditioning')
[Center_Loc2] = ConditionRafts(Center_Loc,Array_Dim);
Center_Loc_temp = Center_Loc2;
disp('Finished raft conditioning')

waitbar(2/4,waitb,['Identifying Rows...'])
Row_pts = {};
count = 0;

while ~isempty(Center_Loc2)
    count = count+1;
    %find top left corner
    TL = [min(Center_Loc2(:,1)),max(Center_Loc2(:,2))]; %min x, max y
    IND = knnsearch(Center_Loc2,TL);  %Find closest point to TL
    PT = Center_Loc2(IND,:); %Finds coordinates of point
    %first row point
    Row = PT;
    Center_Loc2(IND,:) = []; %Remove point from list of points
    %Finding the 2nd point in the row
    IND = 1:size(Center_Loc2,1);
    IND((Center_Loc2(:,1)<=PT(1))|(Center_Loc2(:,2)<PT(2)-W)|(Center_Loc2(:,2)>PT(2)+W)) = []; %Looks for all rafts within +/- 'W' distance in the Y direction and along the row
    if ~isempty(IND)
        [M,IND2] = min(pdist2(PT,Center_Loc2(IND,:))); %Finds the closest point in the row
        IND = IND(IND2); %Finds indice of that point
        PT = Center_Loc2(IND,:); %Finds the coordinatest of that point
    end
    %Find rest of points in the row
    while Row(end,1)<PT(1)
        Row = [Row;PT];
        Center_Loc2(IND,:) = [];
        IND = 1:size(Center_Loc2,1);
        IND((Center_Loc2(:,1)<=PT(1))|(Center_Loc2(:,2)<PT(2)-W)|(Center_Loc2(:,2)>PT(2)+W)) = []; %Looks for all rafts within +/- 'W' distance in the Y direction and along the row
        if ~isempty(IND)
            [~,IND2] = min(pdist2(PT,Center_Loc2(IND,:))); %Uses the closest raft in the row
            IND = IND(IND2);
            PT = Center_Loc2(IND,:);
        end
    end
    
    Row_pts{count} = Row;
end

%Find the order of each row
Row_avg = cellfun(@(x)mean(x(:,2)),Row_pts); %Finds the average y coord of each row
[~,IND] = sort(Row_avg,'descend'); %Sort the rows based on average y coord
Row_pts = Row_pts(IND); %Adds all row points to each row cell
%Making a rows array
S = cellfun(@(x)size(x,1),Row_pts); %Finds the size of each row
Rows = zeros(sum(S),3); %Creates array the size of all the points counted in all rows
SS = cell(size(Row_pts')); %Creates array the size of the rows
for ii = 1:length(S)
    SS{ii} = repmat(ii,S(ii),1); %Assigns row number to each array indice
end
SS = cell2mat(SS); %Converts from cell array to matrix
Rows(:,1:2) = cell2mat(Row_pts'); %Add row coords to array (columns 1 and 2)
Rows(:,3) = SS;  %Adds row number to Rows array in third column

%Finding the columns
waitbar(3/4,waitb,['Identifying Columns...'])
Center_Loc2 = Center_Loc_temp;
Col_pts = {};
count = 0;

while ~isempty(Center_Loc2)
    count = count+1;
    TL = [min(Center_Loc2(:,1)),max(Center_Loc2(:,2))];
    IND = knnsearch(Center_Loc2,TL);
    PT = Center_Loc2(IND,:);
    Col = PT;
    Center_Loc2(IND,:) = [];
    IND = 1:size(Center_Loc2,1);
    IND((Center_Loc2(:,2)>=PT(2))|(Center_Loc2(:,1)<PT(1)-W)|(Center_Loc2(:,1)>PT(1)+W)) = [];
    if ~isempty(IND)
        [M,IND2] = min(pdist2(PT,Center_Loc2(IND,:)));
        IND = IND(IND2);
        PT = Center_Loc2(IND,:);
    end
    
    while Col(end,2)>PT(2)
        Col = [Col;PT];
        Center_Loc2(IND,:) = [];
        IND = 1:size(Center_Loc2,1);
        IND((Center_Loc2(:,2)>=PT(2))|(Center_Loc2(:,1)<PT(1)-W)|(Center_Loc2(:,1)>PT(1)+W)) = [];
        if ~isempty(IND)
            [~,IND2] = min(pdist2(PT,Center_Loc2(IND,:)));
            IND = IND(IND2);
            PT = Center_Loc2(IND,:);
        end
    end
    
    Col_pts{count} = Col;
end

%Averaging the columns
Col_avg = cellfun(@(x)mean(x(:,1)),Col_pts);
[~,IND] = sort(Col_avg,'descend');
Col_pts = Col_pts(IND);
%Making column array
S = cellfun(@(x)size(x,1),Col_pts);
Cols = zeros(sum(S),3);
SS = cell(size(Col_pts'));
for ii = 1:length(S)
    SS{ii} = repmat(ii,S(ii),1);
end
SS = cell2mat(SS);
Cols(:,1:2) = cell2mat(Col_pts');
Cols(:,3) = SS;

waitbar(4/4,waitb,['Labeling Microrafts...'])
Center_Loc2 = Center_Loc_temp;
Max_Cols = length(Col_pts);

%Determine which rows and columns each raft belongs to
[~,~,RIDX] = intersect(Center_Loc2,Rows(:,1:2),'rows','stable');
[~,~,CIDX] = intersect(Center_Loc2,Cols(:,1:2),'rows','stable');
%Assign an index to each raft and sort from low to high
ind = ((Rows(RIDX,3)-1).*Max_Cols)+Cols(CIDX,3);
% Loc = [Loc,ind];
[ind,IDX] = sort(ind);
%Create a matrix containing the raft locations and their assigned index
Raft_locations_real = [Center_Loc2(IDX,:),ind];
%Create an empty matrix based on the number of rows and columns identified
Raft_locations = zeros(max(Rows(:,3))*max(Cols(:,3)),2);

%Averages the redundant microrafts
for ii = ind'
    %Raft_locations(ii,:) = mean(Raft_locations_real((Raft_locations_real(:,3)==ii),1:2));
    Raft_locations(ii,:) = mean(Raft_locations_real((Raft_locations_real(:,3)==ii),1:2),1);
end
disp(size(Raft_locations,1))
% Loc = [Loc,knnsearch(Raft_locations,Center_Loc)];
Loc = [Loc,knnsearch(Raft_locations,Center_Loc)];

close(waitb)