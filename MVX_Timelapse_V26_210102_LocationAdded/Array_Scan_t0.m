function [TIME,L_cell,Centroid_cell,integ,area,cell_count,idx,intensity,integ_inblue,integ_total,cell_loc,camera_error] = Array_Scan_t0(mmc,Positions,m,n,Filters,Exposures,DIR,pixel_size,q,timepoint,p,Size_Filt,Width,Gap,Mask_dye,bfIntensity)
global LED user sys

STRING = ['Array ',num2str(q),': '];
count = 0;
waitb = waitbar(0,'Initializing...');
img = cell(6,1);
timerrr = tic;
camera_error = [];
%Raster scan array
for i = 1:m
    if mod(i,2)
        for j = 1:n
            count = count + 1;
            waitbar(count/(m*n),waitb,[STRING,sprintf('image %d of %d',count,m*n)]);
            if sum(Filters)>0
                powervalue=num2str(bfIntensity);
                fwrite(LED,['SET:LEVEL:CHANNEL1,' powervalue ';SET:LEVEL:CHANNEL2,' powervalue ';'])
            end
            
            mmc.setXYPosition('XYStage',Positions{6}(i,j,1),Positions{6}(i,j,2));
            mmc.setPosition('ZStage',Positions{6}(i,j,3));
            
            if sum(Filters)>0
                while LED.BytesAvailable==0
                end
                while LED.BytesAvailable>0
                    LEDMessage = fscanf(LED);
                end
            end
            
            %If not the first image
            if count>1
                %For scan, saves images
                %While moving to next positoin, run image analysis in
                %parallel 
                if ~isempty(DIR)
                    F(count-1) = parfeval(p,@ANALYZE_IMG_t0,12,img,ii,jj,timepoint,[DIR,'\',num2str(q)],Filters,pixel_size,Size_Filt,Width,Gap,Mask_dye);
                %During re-scan (pre-release of rafts) does not save images
                else
                    F(count-1) = parfeval(p,@ANALYZE_IMG_rescan,3,img,ii,jj,pixel_size,Width,Gap);
                end
            end
            mmc.waitForSystem()
            %
            for k = [6 find(Filters)]
                if k<6
                    powervalue=num2str(0);
                    fwrite(LED,['SET:LEVEL:CHANNEL1,' powervalue ';SET:LEVEL:CHANNEL2,' powervalue ';'])
                    while LED.BytesAvailable==0
                    end
                    while LED.BytesAvailable>0
                        LEDMessage = fscanf(LED);
                    end
                end
                mmc.setPosition('ZStage',Positions{k}(i,j,3));
                SetFilter(mmc,k)
                mmc.setExposure(Exposures(k))
                mmc.waitForSystem()
                [img{k}, camera_error_count] = Snap2(mmc);
                camera_error = [camera_error, camera_error_count];
            end
            %Saves old index for image analysis (F) --> can reference last
            %index while moving the index forward
            ii = i;
            jj = j;
        end
    else
        for j = n:-1:1
            count = count + 1;
            waitbar(count/(m*n),waitb,[STRING,sprintf('image %d of %d',count,m*n)]);
            if sum(Filters)>0
                powervalue=num2str(bfIntensity);
                fwrite(LED,['SET:LEVEL:CHANNEL1,' powervalue ';SET:LEVEL:CHANNEL2,' powervalue ';'])
            end
            
            mmc.setXYPosition('XYStage',Positions{6}(i,j,1),Positions{6}(i,j,2));
            mmc.setPosition('ZStage',Positions{6}(i,j,3));
            
            if sum(Filters)>0
                while LED.BytesAvailable==0
                end
                while LED.BytesAvailable>0
                    LEDMessage = fscanf(LED);
                end
            end
            if count>1
                %                 imshow(img{6},[],'Parent',View)
                if ~isempty(DIR)
                    F(count-1) = parfeval(p,@ANALYZE_IMG_t0,12,img,ii,jj,timepoint,[DIR,'\',num2str(q)],Filters,pixel_size,Size_Filt,Width,Gap,Mask_dye);
                else
                    F(count-1) = parfeval(p,@ANALYZE_IMG_rescan,3,img,ii,jj,pixel_size,Width,Gap);
                end
            end
            mmc.waitForSystem()
            for k = [6 find(Filters)]
                if k<6
                    powervalue=num2str(0);
                    fwrite(LED,['SET:LEVEL:CHANNEL1,' powervalue ';SET:LEVEL:CHANNEL2,' powervalue ';'])
                    while LED.BytesAvailable==0
                    end
                    while LED.BytesAvailable>0
                        LEDMessage = fscanf(LED);
                    end
                end
                mmc.setPosition('ZStage',Positions{k}(i,j,3));
                SetFilter(mmc,k)
                mmc.setExposure(Exposures(k))
                mmc.waitForSystem()
                [img{k},camera_error_count] = Snap2(mmc);
                camera_error = [camera_error,camera_error_count];
            end
            ii = i;
            jj = j;
        end
    end
end

%Turn off LED
powervalue=num2str(0);
fwrite(LED,['SET:LEVEL:CHANNEL1,' powervalue ';SET:LEVEL:CHANNEL2,' powervalue ';'])
while LED.BytesAvailable==0
end
while LED.BytesAvailable>0
    LEDMessage = fscanf(LED);
end

TIME = toc(timerrr);
% imshow(img{6},[],'Parent',View)
%Save and analyze the very last image set
if ~isempty(DIR)
    F(count) = parfeval(p,@ANALYZE_IMG_t0,12,img,ii,jj,timepoint,[DIR,'\',num2str(q)],Filters,pixel_size,Size_Filt,Width,Gap,Mask_dye);
else
    F(count) = parfeval(p,@ANALYZE_IMG_rescan,3,img,ii,jj,pixel_size,Width,Gap);
end

%Creates cell arrays for each parameter for each field of view
L_cell = cell(m,n);
Centroid_cell = cell(m,n);
integ = cell(m,n,5);
area = cell(m,n,5);
cell_count = cell(m,n,5);
idx = cell(m,n,5);
intensity = cell(m,n,5);
integ_inblue = cell(m,n,5);
integ_total = cell(m,n,5);
cell_loc = cell(m,n,5);
%For timelapse scan
if ~isempty(DIR)
    for k = 1:length(F)
        waitbar(k/(m*n),waitb,['Processing Images: ',sprintf('image %d of %d',k,length(F))])
        [i,j,integ_temp,area_temp,cell_count_temp,idx_temp,intensity_temp,integ_inblue_temp,integ_total_temp,L,CENTROID,cell_loc_temp] = fetchOutputs(F(k));
        L_cell{i,j} = L;
        Centroid_cell{i,j} = CENTROID;
        %Fills in parameter cell arrays from parallel jobs with info
        %from F
        for z = find(Filters)
            integ{i,j,z} = integ_temp{z};
            area{i,j,z} = area_temp{z};
            cell_count{i,j,z} = cell_count_temp{z};
            idx{i,j,z} = idx_temp{z};
            intensity{i,j,z} = intensity_temp{z};
            integ_inblue{i,j,z} = integ_inblue_temp{z};
            integ_total{i,j,z} = integ_total_temp{z};
            cell_loc{i,j,z} = cell_loc_temp{z};
        end
    end
else
    %For re-scan (pre-release)
    %Outputs centroid of each raft for re-registering where each raft is
    %for release
    for k = 1:length(F)
        waitbar(k/(m*n),waitb,['Processing Images: ',sprintf('image %d of %d',k,length(F))])
        [i,j,CENTROID] = fetchOutputs(F(k));
%         L_cell{i,j} = L;
        Centroid_cell{i,j} = CENTROID;
    end
end
%Check memory using global counter
[user{timepoint+1}, sys{timepoint+1}] = memory;
memory
delete(F)
clear F
close(waitb)