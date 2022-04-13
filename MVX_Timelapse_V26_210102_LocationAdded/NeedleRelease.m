function NeedleRelease(DIR,Rel_Loc,Raft_locations_cell,Plane_cell,fittp_cell,Positions,Array_Dim,pixel_size,array_select,WandPos,StageWandPos,NeedlePos,c1,c2,r1,r2,WPXCords,WPYCords,MagFocus,appButton,releaseflag,collectflag,releasevideo,arrayType)
global mmc sm1 sm2 s
%% LIST OF HARD CODED PARAMETERS TO CHANGE
%buffer size
% coordinate system is positive X to the right, positive y upwards, (0,0) is in the bottom left

%%
%data=guidata(hObject);
%% Parameters
idlist = [0.75 1 1.2 1 1 1 1 1 1 1 1 1 1 1 1 ]; % displacement multiplier
% idlist = [0.75 0.75 0.75 0.75 0.5 0.5 0.25 0.25 0.5 1 1 1 1 1 1 ];
% idlist = [0.75 1 1 1 0.75 0.75 0.75 0.75 1 1 1 1 1 1 ];
% displthresh = 45/pixel_size; % pixel displacement threshold, if higher, the needle will not try to go farther for fear of off-target
displthresh = 35/pixel_size;
NeedleStep = 0.051; %Step size of needle actuator = 0.051 mm/step
WandStep = 51; %Step size of wand actuator = 0.051 um/step
%Release distance --> Distance below the array to hold the needle (hold
%position)
ReleaseDist = 6.018; %Units are mm; 1900 series motor: 51um per stept
NeedleStroke = 8.313; % Distance needle travels to displace a microraft
NeedleExtension = NeedlePos(3)-ReleaseDist*1000;  %Distance from fully retracted needle to hold position; units are microns
ZCorrection = 0.04323; %FudgeFactor = .04323; z correction factor of steps to microns
xoffset = -87; yoffset = 58; %to fix inprecise calibrations for needle offset  --> To adjust for non-paracentricicty; need to shift to left (-) 100 um
FocusOffset = -0.5; %normally -1
movedlist = NeedleStroke+[0 0 0.1 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2]; % travel distance
aimoffx = 0;aimoffy = 0; %for other release patterns
maxAttempts = 8;
%binning='1x1'; %do not adjust
vidfps = 10;
%Raft geometry  --> load in
gap = Array_Dim(2)/pixel_size;%size dependant
raft = Array_Dim(1)/pixel_size;%size dependant
width = Array_Dim(1)/pixel_size;
debug_window = (Array_Dim(1)+1.5*2.2*Array_Dim(2))/pixel_size; %****analysis window size in microns
half_debug_window = debug_window/2; %****half analysis window size in microns
% % releasevideo = 0; %If true, saves a video of the isolation process
%Save option: if true, saves debugging images and data
saveoption = 1;
% % releaseflag = 0; %If zero then won't actuate needle
% % collectflag = 0; %If zero then won't collect

%% Initializations
dx2s = [];dy2s = []; isaborted = 0;altdx = [];altdy = [];
dxs = []; dys = [];
th = [];
autofocustimes = [];
largemovetimes = [];
needletimes = [];
FindRaftTimes = [];
CollectionTimes = [];
SingleRaftTimes = [];
Totaltime = [];
count = [];
if ~mmc.isSequenceRunning()
    %mmc.setProperty('Binning',binning)
    mmc.initializeCircularBuffer();
    mmc.setCircularBufferMemoryFootprint(250);% // size in MB
    mmc.startContinuousSequenceAcquisition(1);
    pause(0.2)
end
%% initialize timer (potentially possible to not store in data)
%ct = camera timer
ct = timer;
ct.BusyMode = 'drop';
ct.ExecutionMode = 'fixedSpacing';
ct.Period = 1/vidfps;
ct.TasksToExecute = Inf;
%p=plot(data.gui.plotax,[1],[1],'b*');
camera_count = 1;
image_data = [];
mkdir([DIR,'\ReleaseData']);
folder = [DIR,'\ReleaseData'];% save folder directory ("pdir")
mkdir([DIR,'\ReleaseData\ReleaseVideo']);
ct.UserData = {camera_count, image_data releasevideo, DIR}; %current timer count, last circular buffer frame, video saving flag.
ct.TimerFcn = @step_video;
start(ct) 

%% more initializations
%Make a directory to save release data into, such as images etc.
raftpix = raft+gap;%size dependant
analysis_window = raftpix/2;%size dependant
half_analysis_window = analysis_window/2;
%Windows for analyzing differences between images
diff_window = width+sqrt(2)*width;%size dependant
half_diff_window = diff_window/2;
%Set initial row
r = r1;
c = c1;
%Create waitbar
waitb = waitbar(0,'','Position',[75, 700, 280, 60]);
%% build list of target rafts
for q = array_select %Array number list --> loaded from Rel_Loc.mat
    clist = Rel_Loc(find(Rel_Loc(:,2)==q),1); % load index numbers of target microrafts for this array
    well_address = cell(size(clist,1),3);
    %% initialize needle to base, then raise to position below array
    if releaseflag && NeedleExtension >= 14600
        [~] = NeedleZero(NeedleExtension);
        pause(0.5);
    else
        releaseflag = 0;
        %f = msgbox('Needle Extension Value Not Valid: Releases Cancelled');
    end
    Totaltimer=tic; % total process timer
    %% Main loop
    clistSize = size(clist,1);
    for i=1:size(clist,1)
        %create cell arrays to hold all images
        preRelease = cell(1,maxAttempts);
        postRelease = cell(1,maxAttempts);
        detectionMask = cell(1,maxAttempts);
        diffRemain = cell(1,maxAttempts);
        diffGone = cell(1,maxAttempts);
        all_images = [];
        grid_size = [3,maxAttempts];
        %Display waitbar
        waitbar((i/clistSize),waitb,['Isolating Raft ',num2str(i),' of ',num2str(clistSize)]);
        SingleRaftTimer = tic;
        % Reset needle Z every 50 targets
        if mod(i,50)==0
            %reset needle
            if releaseflag  && NeedleExtension>= 14600
                [~] = NeedleZero(NeedleExtension);
                pause(1);
            end
        end
        % AUTOFOCUS GOES HERE -- if updating periodically
        
        %% curve-correction for needle translation distances
        if i==1
            %Find the travel distance to add to compensate for play in
            %membrane
            [needle_Travel_Distances] = FlippedArrayCurve(Positions{q}(:,:,1),Positions{q}(:,:,2),Positions{q}(:,:,3),ZCorrection);
            xP = Positions{q}(:,:,1);
            yP = Positions{q}(:,:,2);
            %creates a thin plate spline to be used for determining needle
            %release distances
            [fittp_needle, ~] = fit( [xP(:), yP(:)], needle_Travel_Distances(:), 'thinplateinterp' );
        end
        %% Enter Main loop
        RaftIndex = clist(i); %the ith index number for this array
        %ll=Raft_locations_cell{1}; % x and y position
        timestr = ['Index ' num2str(Rel_Loc(i,1)) ', Array ',num2str(Rel_Loc(i,2))]; %save string
        if ~isequal(Raft_locations_cell{q}(RaftIndex,:),[0,0])
            %Indicates that the raft coordinates were identified
            RaftFound = 1;
            %% Find microraft and go to it
            t=tic; %individual release timer
            micron_loc = Raft_locations_cell{q}(RaftIndex,:);% the x and y coordinate of the ith target for the qth array
            Needle_Loc = NeedlePos(1:2)*pixel_size; %calibxy= "NeedlePos(1:2)" coordinate system is for: show image with imshow get(hObject,'CurrentPoint');
            %Time it takes to move the stage and focus
            largemovetimer = tic;
            newposx = round(micron_loc(1) - (Needle_Loc(1) - 1024*pixel_size) + xoffset);
            newposy = round(micron_loc(2) - (Needle_Loc(2) - 1024*pixel_size) + yoffset);
            mmc.setXYPosition('XYStage',newposx,newposy);
            zn = fittp_cell{q}(newposx,newposy);
            mmc.setPosition('ZStage',zn+FocusOffset); % offset helps with rafts partially released
            mmc.waitForSystem
            largemovetimes(i) = toc(largemovetimer);
            
            %% Autofocus
            %Run autofocus after 10th release
%             if i == 10
%                 CurrentZ = mmc.getPosition('ZStage');
%                 mmc.waitForSystem;
%                 %Add in Mag code
%                 %Run rough autofocus code on current position
%                 [ZF] = AFocusIsolation(CurrentZ, Array_Dim(1), Array_Dim(2), pixel_size, MagFocus);
%                 %Only set new focus if within +/- 15 "microns" on prior cont.
%                 if ZF<(CurrentZ+10) && ZF >(CurrentZ-10)
%                     mmc.setPosition('ZStage',ZF); % Set z stage to ZF
%                     mmc.waitForSystem
%                 end
%             end  
            %% Release loop
            count(i) = 1; %needle attempts
            met(i,1) = 1; %fraction of microraft well occupied by identified microraft material
            met2(i,1) = 0; % fraction of removed microraft material
            met3(i,1) = 1; % fraction of microraft material remaining
            diffdoneflag = 0; %experimental metric probably not useful
            decision=true;%decision= %fasle means done, true means release more
            while decision
                %% Allow for abort
                %Check whether isolation has been aborted --> app.Button.Value will
                %be equal to 1 after pressing begin and 0 after pressing abort
                isaborted = ~appButton.Value;
                if isaborted
                    break
                end
                %% Calculate a retarget (only for second poke or later)
                if count(i)== 1
                    dx = 0;
                    dy = 0;
                else
                    %Finds rafts that are partially released --> if they
                    %have moved from a previous puncture
                    [dx,dy] = modifyRetarget(littlei,th,i,count,idlist,displthresh); % compute the displacements dx and dy that are required to adjust the microneedle from previous to the actual microraft location
                    dxs(i,count(i)) = dx; %Check to see if recording data
                    dys(i,count(i)) = dy;
                end
                %Move to adjusted position after locating new raft centroid
                mmc.setXYPosition('XYStage',newposx+pixel_size*dx,newposy+pixel_size*dy);
                mmc.waitForDevice('XYStage');
                
                %% Initial drift compensation, uses geometry-specific segmentation to find a very accurate initial microraft location
                if count(i) == 1
                    %% gather 3 new images
                    tempd = get(ct,'Userdata');%acquire timer userdata, containing the last image frame (as well as some other stuff)
                    imgc = tempd(2); % access, specifically, the last image frame
                    img = imgc{1}; % access, specifically, the cData of the last image frame
                    notchanged = 3; %number of fresh images needed to acquire --> might be able to change to 2
                    fixstopcount = 0;
                    while notchanged~=0 % until 3 "fresh" images are acquired, continue pinging buffer for last frames
                        tempd = get(ct,'Userdata');
                        imgc = tempd(2); %Image data
                        imgn = imgc{1}; %image data --> read from cell
                        if ~isequal(img, imgn)
                            notchanged = notchanged-1;
                        end
                        img = imgn;
                    end
                    img = rot90(img,1) ;%flipud(double(img));
                    
                    bigio = img; % Big image original; hold onto initial timepoint image
                    %% get threshold and save pre-release image, 300 x 300 pix centered on needle/raft (should be coaxial)
                    if arrayType == 1  %Standard arrays
                        tempth = multithresh(img,2);
                    elseif arrayType == 2  %Clear Microraft arrays 100x100x200
                        tempth = multithresh(img,2); %Normal thresh: 2 -->Increasing levels decreases sensitivity
                    else   %Clear microraft arrays 100x100x60
                        tempth = multithresh(img,4);
                    end
                    %tempth = multithresh(img,3);
                    %%th(i)=tempth(2);
                    th(i) = tempth(1);
                    %Set the curent threhsold to the medium (of non-zero
                    %values)
                    th(i) = median(th(th>0));
                    if saveoption
                        si = imcrop(img,[ NeedlePos(1)-half_debug_window-dx NeedlePos(2)-half_debug_window-dy debug_window debug_window]);
                        %imwrite(uint16(si), [folder '\Pre-Release ' timestr '.png'])
                        preRelease{1,1} = si;
                        all_images = imtile({preRelease{1,1:maxAttempts}, postRelease{1,1:maxAttempts}, detectionMask{1,1:maxAttempts}} , 'GridSize', grid_size, 'BorderSize', [3 3]);
                        imwrite(all_images, [folder '\' timestr '.tif'])
                    end
                    FindRaftTimer = tic;
                    %% Need to determine dx and dy required to shift the needle to the actual microraft location. tt= [dx,dy]
                    tt = findRaftInCenter(img,th(i),NeedlePos(1:2),Array_Dim(1),Array_Dim(2),pixel_size);
                    FindRaftTimes(i,count(i)) = toc(FindRaftTimer);
                    %% Shift centroid via geometry
                    %Check to make sure tt is not empty or does not have
                    %too many points
                    if isempty(tt) || size(tt,1)>1
                        %don't adjust
                        'couldnt adjust'
                    else
                        %If shifts are good then create new shifts
                        dx2 = tt(1);dy2 = tt(2); %the shift
                        dx2s(i,count(i)) = dx2; %keep track of shift
                        dy2s(i,count(i)) = dy2;
                        %% check for raft, maybe useful to ignore releases altogether
                        newposx = newposx+(pixel_size*dx2);
                        newposy = newposy+(pixel_size*dy2);
                        mmc.setXYPosition('XYStage',newposx,newposy);
                        mmc.waitForSystem
                    end
                end
                %% Release Raft
                %Be very careful --> dt is in millimeters
                %add in limits of no more than 4 millimeters (or something)
                dt = (movedlist(count(i))+(fittp_needle(mmc.getXPosition('XYStage'),mmc.getYPosition('XYStage'))/1000)); % compute needle translation distance for this particular target %Distance needs to be in mm
                moved = round(dt/NeedleStep); %0.051mm/step
                needletimer = tic;
                %Do not release if release distance is greater than 9
                if releaseflag && dt<9
                    move(sm2,-(moved)); %Extend (-)
                    move(sm2,moved);  %Retract (+)
                    release(sm2);
                end
                needletimes(i,count(i)) = toc(needletimer);
                %% Image Saving and detecting release
                pause(0.3) % allow raft to float up
                mmc.waitForSystem
                %%  Acquire final image
                tempd = get(ct,'Userdata');
                imgc = tempd(2);
                img = imgc{1};
                fixstopcount = 0;
                notchanged = 1;
                while notchanged & fixstopcount<100
                    tempd = get(ct,'Userdata');
                    imgc = tempd(2);
                    imgn = imgc{1};
                    if ~isequal(img, imgn)
                        notchanged=0;
                    end
                    img = imgn;
                    fixstopcount = fixstopcount+1;
                end
                img = rot90(img,1) ;%flipud(double(img));
                %%      Save Post-Release image
                diffi = medfilt2(double((img)))-medfilt2(double(imtranslate(bigio,[-dx-dx2,-dy-dy2]))); % dy2 may need negation
                diffi = imcrop(diffi,[ NeedlePos(1)-half_diff_window-dx NeedlePos(2)-half_diff_window-dy diff_window diff_window]);
                diffmaskgone = diffi>1000;%need the threshold you would use for a backgorund-subtracted image
                diffmaskremain = diffi<-1000;%thdiff(i);
                if saveoption
                    si = imcrop(img,[ NeedlePos(1)-half_debug_window-dx NeedlePos(2)-half_debug_window-dy debug_window debug_window]);
                    %imwrite(uint16(si), [folder '\Post-release ' timestr ' Try ' num2str(count(i)-1) '.png']) %Make sure to check image to make sure it is an image from the intended timepoint of acquisition
                    %Save the difference file --> used for debugging
                    %save([folder '\' timestr ' Try ' num2str(count(i)-1) ' diffi.mat'],'diffi')
                    %Difference image saved to show material removed
                    %imwrite(uint16(diffmaskgone), [folder '\DiffGone ' timestr ' Try ' num2str(count(i)-1) '.tif'])
                    %Difference image saved to show material remaining
                    %imwrite(uint16(diffmaskremain), [folder '\DiffRemain ' timestr ' Try ' num2str(count(i)-1) '.tif'])
                    postRelease{1,count(i)} = si;
                    diffRemain{1,count(i)} = diffmaskremain;
                    diffGone{1,count(i)} = diffmaskgone;
                
                end
                count(i) = count(i)+1;
                %%      Crop Image and Save
                %Create cropped post-release image
                littlei = imcrop((img),[ NeedlePos(1)-analysis_window-dx NeedlePos(2)-analysis_window-dy raftpix raftpix]);
                %Apply area metrics to determine remaining/removed/gone
                met(i,count(i)) = sum(littlei(:)<th(i))/(raftpix*raftpix);
                met2(i,count(i)) = sum(diffmaskgone(:))/(debug_window*debug_window);
                met3(i,count(i)) = sum(diffmaskremain(:))/(debug_window*debug_window);
                %Another possible option: dermine whether centroid of raft
                %is outside of well ---> see Matt's original NeedleRelease
                if saveoption
                    timevec = clock;
                    %Save the difference segmented image mask
                    %imwrite(littlei<th(i), [folder '\Detection mask ' timestr ' Try ' num2str(count(i)-1) '.png'])
                    detectionMask{1,(count(i)-1)} = littlei<th(i);
                    all_images = imtile({preRelease{1,1:maxAttempts}, postRelease{1,1:maxAttempts}, detectionMask{1,1:maxAttempts}} , 'GridSize', grid_size, 'BorderSize', [3 3]);
                    imwrite(all_images, [folder '\' timestr '.tif'])
                end
                %% make a decision.
                %             if met3(i,count(i))<0.01 & met2(i,count(i))>0.02 % if the remaining mat is really low
                %                 diffdoneflag=1;
                %             end
                %Display attempt number in command window
                AttemptNumber = count(i)-1
                %Make a decision
                if count(i)>maxAttempts
                    decision = false;
                end
                %Stop release attempts if less than 4% of segmented raft is
                %left AND if at least three attempts have been made
                if met(i,count(i))<0.02 && count(i) > 5
                    decision = false;
                end
                %             if (diffdoneflag & count(i)>3)
                %                 decision=false;
                %             end
            end
            tims(i) = toc(t);
            save([folder '\ReleaseStatistics at ' timestr '.mat'],'count','tims','met','met2','met3','th','clist','dx2s','dy2s','largemovetimes','needletimes','camera_count','FindRaftTimes','dxs','dys')
        else
            'Raft coordinates not found: moving to next raft'
            RaftFound = 0;
        end
        %Check whether isolation has been aborted --> app.Button.Value will
        %be equal to 1 after pressing begin and 0 after pressing abort
        isaborted = ~appButton.Value;
        if isaborted
            break
        end
        %%
        %If true, collect microrafts
        if collectflag && RaftFound == 1
            if r>r2
                r = r1;
                c = c+1;
            end
            if c>c2
                z = warndlg('The end of the well plate has been reached! Please insert another plate.');
                uiwait(z)
                r = r1;
                c = c1;
            end
            %Record the row, column, and array adress of isolated raft
            well_address(i,1:3) = {r,c,timestr};
            %%% ORIENT PICKANDPLACE WAND ~~~~~~~~~~~~~~~~~~~~~~
            %Start timer
            CollectionTimer = tic;
            %Move wand to center of FOV and extend wand into media
            %clear raft location coordinates
            currentX_temp = [];
            currentY_temp = [];
            currentX = [];
            currentY = [];
            
            currentX_temp = Raft_locations_cell{q}(RaftIndex,1);
            currentY_temp = Raft_locations_cell{q}(RaftIndex,2);
            %For large wands, keep wand within array boundary with
            %correctoin
            [currentX,currentY] = WandCorrection(currentX_temp, currentY_temp, Raft_locations_cell);
            %Calculate position for wand; round to a multiple of 3 (3
            %microns/step
            MX_temp = WandPos(5,1)-(currentX-StageWandPos(1,1));
            RatioX = floor(MX_temp/3);
            MX = RatioX*3;
            %Calculate position for wand; round to a multiple of 3 (3
            %microns/step
            MY_temp = WandPos(5,2)+(currentY-StageWandPos(1,2));
            RatioY = floor(MY_temp/3);
            MY = RatioY*3;
            if MX >=0 && MX <=300000 && MY >=0 && MY<=210000
                MoogMove(MX, MY,s)
                pause(3)
                %Extend wand into media, pause to capture raft, retract
                %(-)Extend wand (+) Retract wand
                WandMove = WandPos(5,3)/WandStep;
                move(sm1,-WandMove);
                release(sm1);
                pause(15);     %25     %%%%%%%%%%%% Wait Time %%%%%%%%%%%%%%%%%
                sm1.RPM = 7; %Slow wand retraction to keep raft: Speed (fast-slow:1x-6x): 8,2,5,1,7
                move(sm1,WandMove);
                release(sm1);
                sm1.RPM = 2; %Increase wand extension for deposition
            else
                'Move beyond Moog (xy) actuator limits'
            end

            %Move actuators to well plate position --> row and column #
            %Adjust X and Y to be multiples of 3 (3 microns/step for Moog)
            MX_temp = floor(WPXCords(r));
            RatioX = floor(MX_temp/3);
            MX = RatioX*3;
            MY_temp = floor(WPYCords(c));
            RatioY = floor(MY_temp/3);
            MY = RatioY*3;
            
            if MX >=0 && MX <=300000 && MY >=0 && MY<=210000
                MoogMove(MX, MY,s)
                pause(2);
                %Extend wand into media, wait, then retract: (-)Extend wand (+) Retract wand
                WandMove = WandPos(1,3)/WandStep;
                move(sm1,-WandMove);
                release(sm1);
                pause(12);     %12    %%%%%%%%%%%% Wait Time %%%%%%%%%%%%%%%%%
                move(sm1,WandMove);
                release(sm1);
                %Increase row count
                r = r+1;
            else
                'Move beyond Moog (xy) actuator limits'
            end

            %Move wand to home position to avoid lighting issues during
            %release
            if MX >=0 && MX <=300000 && MY >=0 && MY<=210000
                MoogMove(0,0,s)
                pause(1.5)
            else
                'Move beyond Moog (xy) actuator limits'
            end
            CollectionTimes(i) = toc(CollectionTimer);
        end
        SingleRaftTimes(i) = toc(SingleRaftTimer);
        %Check whether isolation has been aborted --> app.Button.Value will
        %be equal to 1 after pressing begin and 0 after pressing abort
        isaborted = ~appButton.Value;
    end
    %Save release statistics 
    if ~isaborted
        save([folder '\ReleaseStatistics at end for array ',num2str(q),'.mat'],'count','tims','met','met2','met3','th','clist','dx2s','dy2s','largemovetimes','needletimes','camera_count','FindRaftTimes','CollectionTimes','SingleRaftTimes','dxs','dys')
    end
    
    save([folder '\Well Addresses of Isolated Rafts for Array ', num2str(q),'.mat'], 'well_address');
    %Clear variables
    count= []; tims=[];met=[];met2=[];met3=[];th=[]; clist=[]; dxs=[];dys=[];dx2s=[];dy2s=[];largemovetimes=[];needletimes=[];
    camera_count = 1;
    image_data = [];
    %show i
    CurrentRaftNumber = i
    %If releasing from more than one array, give a warning dialogue for
    %next array and pause until click 'okay'
    if size(array_select,2)>1 && q <array_select(end)
        warning1 = warndlg(['Finished with releases from Array ',num2str(q), ', Click OK when ready for next array']);
        uiwait(warning1);
    end
end
close(waitb)
Totaltime = toc(Totaltimer);
%% Clean up/Save results/abort
%Zero needle and then set to an extension value of ~ 20 mm --> convert to
%microns before 
% if releaseflag && NeedleExtension>= 14600
%     WandFinalMove = (NeedleExtension-14600)/WandStep;
%     move(sm1,-WandFinalMove);
%     release(sm1);
% end
%% stop feed
data.releasevideo=0;
tempt = ct.UserData;
stop(ct)
delete(ct)
mmc.stopSequenceAcquisition();
%mmc.setProperty('HamamatsuHam_DCAM','Binning','1x1')
camera_count_temp = tempt(1); % if you wanted to get the frame, it would be tempt(2)
camera_count = camera_count_temp{1}; % you always need to get out the first cell element to actually get the original variable
%save([folder],'camera_count');

'Isolation complete'
end