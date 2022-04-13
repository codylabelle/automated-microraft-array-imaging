function [j,k,integ_temp,area_temp,cell_count_temp,idx_temp,intensity_temp,integ_inblue_temp,integ_total_temp,L,CENTROID,cell_loc_temp] = ANALYZE_IMG_t0(img,j,k,i,DIR,Filters,pixel_size,Size_Filt,Width,Gap,Mask_dye)

%Save into direction j by k
DIR2 = [DIR,'\',num2str(j),'x',num2str(k)];
%Save brightfield image
imwrite(img{6},[DIR2,'\',num2str(i,'%2.2d'),'_BF.tif'])
%Create binary mask
L_bw = uraft(img{6},Width,Gap,pixel_size);
%Labeling each raft with a different integer value (left to right)
L = bwlabel(L_bw);
%Find the centroid and bounding box of each raft
ZL = regionprops(L_bw,'Centroid','BoundingBox');
CENTROID = round(cell2mat({ZL.Centroid}'));
%Max number of rafts in each image
L_max = size(CENTROID,1);

%Create max and min pixel quantities based on upper and lower cell limits
%and clump limits
cell_min = floor(pi*(Size_Filt(1)/pixel_size)^2);
cell_max = floor(pi*(Size_Filt(2)/pixel_size)^2);
clump_max = floor(pi*(Size_Filt(3)/pixel_size)^2);
cell_rad = round(Size_Filt(2)/pixel_size);

%Create arrays to hold parameters, for the five channels
integ_temp = cell(5,1);         %Measure #1: integrated intensity of thresholded cells
area_temp = cell(5,1);          %Measure #2: total area of thresholded cells
cell_count_temp = cell(5,1);    %Measure #3: cell count
idx_temp = cell(5,1);           %Measure #4: individual thresholded pixel locations
intensity_temp = cell(5,1);     %Measure #5: individual thresholded pixel intensities
integ_total_temp = cell(5,1);   %Measure #6: integrated intensity of entire raft and contents
integ_inblue_temp = cell(5,1);  %Measure #7: Integrated intensity of fluorescence within chosen mask
cell_loc_temp = cell(5,1);      %Measure #8: Cell location from fastpeakfind


%Starts for a for loop for all channels and puts mask dye first
for z = unique([Mask_dye,find(Filters)],'stable')
    %Save with particular name rather than number
    switch z
        case 1
            STRING = 'Blue';
        case 2
            STRING = 'Violet';
        case 3
            STRING = 'Green';
        case 4
            STRING = 'Red';
        case 5
            STRING = 'FarRed';
    end
    %Pick an image
    img2 = img{z};
    %Save the image
    imwrite(img2,[DIR2,'\',num2str(i,'%2.2d'),'_',STRING,'.tif'])
    %Apply weiner filter to reduce noise
    %%img2 = wiener2(img2,[round((Size_Filt(2)/2)/pixel_size),round((Size_Filt(2)/2)/pixel_size)]);
    %img2 = wiener2(img2,[round((Size_Filt(1)*3)/pixel_size),round((Size_Filt(1)*3)/pixel_size)]);
    img2 = wiener2(img2,[round((Width/10)/pixel_size),round((Width/10)/pixel_size)]);
    %Find max value of image
    M = double(max(img2(:)));
    %Find median of image (may use for threshold)
    thresh1 = median(img2(:));
    %Apply top hat filter --> Removes objects larger than SE
    %The image minus the morphological opening of the image
    %Remove any brightspots that are larger than the structuring element
    %se = strel('disk',ceil((Size_Filt(1)*15)/pixel_size));
    se = strel('disk',ceil((Width/6.5)/pixel_size));
    img2 = imtophat(img2,se);
    %Find multiplier
    mult = M/max(img2(:));
    %Apply multiplier to get back to original max
    img2 = uint16(img2*mult);
    %Set all pixels not within raft boundaries to zero
    img2 = (img2.*uint16(L_bw));

    %Segment cells -- two level otsu's method
    thresh2 = multithresh(img2,2);
    %Comparing median of denoised image vs otsu's two level threshold
    THRESH = [thresh1,thresh2];
    %If Otsu has a level less than X% of the median, then make infinite
    %so that it is never chosen
    THRESH(THRESH<=(0.1*thresh1)) = Inf;
    %Pick the lowest threshold for best sensitivity
    bw = img2>min(THRESH);    
    
    %Remove any bright objects smaller than the lower cell size
    if cell_min>1
        bw = bwareaopen(bw,cell_min);
        %Always remove single pixel objects
    else
        bw = bwareaopen(bw,1);
    end
    if cell_max>0
        bw = bw&~bwareaopen(bw,clump_max);
        %Filter out edges of pdms from segmented objects
        ZL2 = regionprops(bw,'MajorAxisLength','Eccentricity','PixelIdxList','Solidity');
        removeObjects = [];
        keepObjects = false(2048,2048);
        for ii = 1:length(ZL2)
            if ZL2(ii).MajorAxisLength >(Width/2/pixel_size) && ZL2(ii).Eccentricity > 0.95 || ...
                    ZL2(ii).MajorAxisLength >(Width/2/pixel_size) && ZL2(ii).Solidity < 0.5 && ZL2(ii).Eccentricity > 0.8 || ...
                    ZL2(ii).Eccentricity > 0.98
                removeObjects = [removeObjects; ii];
            else
                keepObjects(ZL2(ii).PixelIdxList) = 1;
            end
        end
        bw = keepObjects;
    end
    %If you chose a mask dye
    if ~isempty(Mask_dye)
        %And is the mask dye
        if z == Mask_dye
            img_temp = img2.*uint16(bw);
            %Find peaks and indices
            p = FastPeakFind(img_temp);
            %p = FastPeakFind(imgaussfilt(img_temp,(round((Size_Filt(1)*2)/pixel_size)));
            R = p(2:2:end);
            C = p(1:2:end);
            ind = sub2ind(size(img2),R,C);
            bw2 = false(size(img2));
            %Use indices to make discs the size of a cell for the mask
            bw2(ind) = true;
            bw2 = imdilate(bw2,strel('disk',cell_rad));
        end
    else
        bw2 = false(size(bw));
    end
    
    %Within each parameter create a matrix for the number of rafts
    %Initilize variables
    %Integrates/sums pixel values for all cells on raft -->
    integ_temp{z} = zeros(L_max,1);
    %Find the number of cells on each raft
    cell_count_temp{z} = zeros(L_max,1);
    %Find the total area of fluorescence on each raft
    area_temp{z} = zeros(L_max,1);
    %Find indices of all pixels on each raft
    idx_temp{z} = cell(L_max,1);
    %Find the total fluorescence intensity on each raft
    intensity_temp{z} = cell(L_max,1);
    %Find the total fluorescence within the boundaries of the Hoechst
    %nuclear stain mask
    integ_inblue_temp{z} = zeros(L_max,1);
    %Find the total fluoresence on each raft --> all fluorescence,
    %including autofluoresence
    integ_total_temp{z} = zeros(L_max,1);
    %Find cell location based on peak fluorescence
    cell_loc_temp{z} = cell(L_max,1);

    %For every raft -->
    for kk = 1:L_max
        %crop images to only look at single rafts and contents
        %Create a cropped mask for each  raft for brightfield image -->
        %black and white image of single microraft
        mask = imcrop(L_bw,ZL(kk).BoundingBox);
        %Create a cropped image of each raft for fluoresence image -->
        %grayscale image of microraft
        img_temp = imcrop(img2,ZL(kk).BoundingBox);
        %Create a mask of each cell on each raft --> black and white
        %image of cells on microraft
        bw_temp = imcrop(bw,ZL(kk).BoundingBox);
        %Create a mask of each cell on each raft for the given
        %selected overall mask dye
        bw_temp2 = imcrop(bw2,ZL(kk).BoundingBox);
        %Use regionprops function to find fluorescence area, pixel
        %values, and and index list --> uses raft regions (bw_temp)
        %within image (img_temp)
        Z = regionprops(bw_temp,img_temp,'Area','PixelValues','PixelIdxList');
        %Integrates/sums pixel values for all cells on raft -->
        integ_temp{z}(kk) = sum(cat(1,Z.PixelValues));
        %Finds total area of bright pixels on raft
        area_temp{z}(kk) = sum(cat(1,Z.Area));
        %Creates indices of pixels on each raft
        idx_temp{z}{kk} = cat(1,Z.PixelIdxList);
        %Find the individual intensity of each pixel that is
        %fluorescing
        intensity_temp{z}{kk} = cat(1,Z.PixelValues);
        %Find totel pixel intensity within the selected mask
        integ_inblue_temp{z}(kk) = sum(img_temp(bw_temp2));
        %Find the total fluoresence on each raft --> all fluorescence;
        %particles, autofluoresence, etc.; no thresholding
        %Useful for raft packed with cells or dim cells
        integ_total_temp{z}(kk) = sum(img_temp(mask));
        %Create image that has a black background with cells in
        %grayscale --> condition image for FastPeakFind
        img_temp2 = img_temp.*uint16(bw_temp);
        %Find the number of cells on each raft using FastPeakFind
        if z == 3 || z == 4
            p = FastPeakFind(imgaussfilt(img_temp2,(round((Size_Filt(1)/2)/pixel_size))));
        else
            p = FastPeakFind(img_temp2);
        end
        cell_count_temp{z}(kk) = length(p)/2;
        R = p(2:2:end);
        C = p(1:2:end);
        ind = sub2ind(size(img2),R,C);
        cell_loc_temp{z}{kk} = ind;
    end
end
%L_max is the number of rafts in the image
if L_max<256
    L = uint8(L);
elseif L_max<65536
    L = uint16(L);
end