function displacements = findRaftInCenter(img,th,calibxy,width,gap,pixel_size)
%% segment image
% Output lm =  [displacementX,displacementY]
%Creaet binary image of full image (full FOV)
[bw] = uraftIsolation(img,width,gap,pixel_size);
%Find centroid of all objects
st=regionprops((bw),'Centroid');
cs=vertcat(st(:).Centroid);

%Create indices for expected location of raft
expect=[calibxy(1),calibxy(2)];
if ~isempty(cs)
    ds=sqrt( (expect(1)-cs(:,1)).^2 + (expect(2)-cs(:,2)).^2);
    %Find raft with centroid within 30 pixels distance from expected
    %location
    id=find(ds<(40/pixel_size));
    if ~isempty(id)
    displacements=cs(id(1),:)-expect;
    %If no raft is within 20 distance then no shift
    else
        'no raft within 40 microns of needle'
        displacements=[0 0];
    end
else
    %If no centroids are found then no shift
    'no centroids found'
    displacements=[0 0];
end
