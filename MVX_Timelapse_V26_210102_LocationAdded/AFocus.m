function [fittp_BF,fittp_CTL,Plane,POS_array] = AFocus(Plane,Mag,correction_factor,pixel_size,Width,Gap,POS,Filters,Exposures,objective,bfIntensity)
global Div BB mmc LED zFScore

POS_array = cell(6,1);
flag = ones(6,1);
%make sure there are at least 3 rafts in each field of view --> cropping
%image so you only see a 3x3 --> minimum for good focus score (speeds up
%analysis)
Div = floor(2048/(Width+Gap)*pixel_size/3);

switch objective
    %1x objective
    case 1
        NA = 0.25;
        MagArray = [0.63,0.8,1,1.25,1.6,2,2.5,3.2,4,5,6.3];
    %2x objective
    case 2
        NA = 0.5;
        MagArray = 2*[0.63,0.8,1,1.25,1.6,2,2.5,3.2,4,5,6.3];
end

Lambda = .45;
e = 6.5;
n = 1.33;
%Fudge factor for z stage movement; correction factor for depth of well
%FudgeFactor = .023;
FudgeFactor = .04323;
correction = correction_factor*FudgeFactor*2;

M = MagArray(Mag);
%equation for depth of focus for a microscope
df = (((Lambda*n)/(NA^2))+((n*e)/(M*NA)))*FudgeFactor;

if length(Filters)>1
    test_filters = find(Filters(1:5));
    if Filters(6)>0
        test_filters = [6 test_filters];
    end
else
    test_filters = Filters;
end

string_cell = {'Autofocus on top left corner...';'Autofocus on top right corner...';'Autofocus on bottom right corner...';'Autofocus on bottom left corner...';'Autofocus on center...'};
new_planes = cell(5,1);
waitb = waitbar(0);
%For brightfield
for num = 1:length(string_cell)
    waitb = waitbar(num/6,waitb,string_cell{num});
    mmc.setXYPosition('XYStage',Plane(num,1),Plane(num,2));
    mmc.setExposure(Exposures(6));
    SetFilter(mmc,0);
    %turn on LED
    powervalue=num2str(bfIntensity);
    fwrite(LED,['SET:LEVEL:CHANNEL1,' powervalue ';SET:LEVEL:CHANNEL2,' powervalue ';'])
    pause(0.5)
    while LED.BytesAvailable>0
        LEDMessage = fscanf(LED);
    end
    %Uses optimization function to find the optimum z plane for focusing
    opt = optimset('MaxFunEval',20,'Display','off','TolX',4*df);
    %Experimentally found minimums and maximums of z plane
    %Note that moving 4 df's doesn't change much, so 40*df is about 10*df
    %experimentally
%     Z_lower = Plane(num,3)-40*df;
%     Z_upper = Plane(num,3)+40*df;
    Z_lower = Plane(num,3)-20*df;
    Z_upper = Plane(num,3)+20*df;
    BB = [];
    mmc.waitForSystem()
    Flag_BF = 0;
    %While an actual minimum has not been found, Flag_BF is 0 --> 
    %ensure that you find the true minimum
    while Flag_BF < 1
        [ZF,~,Flag_BF] = myfminbnd(@FScore,Z_lower,Z_upper,opt);
        Z_lower = Z_lower-30*df;
        Z_upper = Z_upper+30*df;
        opt.MaxFunEvals = opt.MaxFunEvals*2;
    end
    zfOffset = 2;  %Corrects for curve offset from best focus plane
    ZF = ZF + zfOffset;
    %Sets the z position to the optimal z plane, ZF
    mmc.setPosition('ZStage',ZF)
    Plane(num,3) = ZF;
    mmc.waitForSystem()
    img = Snap2(mmc);
    

    %uraft segments the microraft --> creates a black backround with white
    %squares where the rafts are
    bw = uraft(img,Width,Gap,pixel_size);
    seBW = strel('disk',round((Width/4)/pixel_size));
    bw = imerode(bw,seBW);
    %Crops out individual microrafts  --> for fluorescence analysis
    BB = regionprops(bw,'BoundingBox');
    mmc.waitForSystem()
    
    New_corner = repmat(Plane(num,:),6,1);
    %Uses optimization function to find the optimum z plane for focusing
    opt = optimset('MaxFunEval',40,'Display','off','TolX',df);
    Z_lower = Plane(num,3)-20*df;
    Z_upper = Plane(num,3)+40*df;
    %Find optimum z plane for focusing in each channel selected
    for z = test_filters
        zFScore = z;
        mmc.setExposure(Exposures(z))
        if z==6
            SetFilter(mmc,0);
        else
            powervalue=num2str(0);
            fwrite(LED,['SET:LEVEL:CHANNEL1,' powervalue ';SET:LEVEL:CHANNEL2,' powervalue ';'])
            pause(0.5)
            while LED.BytesAvailable>0
                LEDMessage = fscanf(LED);
            end
            SetFilter(mmc,z);
        end
        mmc.waitForSystem()
        %Optimization of z plane within the segemented rafts
        [ZF,~,flag(z)] = myfminbnd(@FScore,Z_lower,Z_upper,opt);
        if flag<1
            opt = optimset('MaxFunEval',40,'Display','off','TolX',df);
            Z_lower = Plane(num,3)-30*df;
            Z_upper = Plane(num,3)+40*df;
            [ZF,~,flag(z)] = myfminbnd(@FScore,Z_lower,Z_upper,opt);
            %Reset upper and lower bounds for next filter
            opt = optimset('MaxFunEval',40,'Display','off','TolX',df);
            Z_lower = Plane(num,3)-20*df;
            Z_upper = Plane(num,3)+40*df;
        end
        mmc.waitForSystem()
        New_corner(z,3) = ZF;
    end
    new_planes{num} = New_corner;
end

powervalue=num2str(0);
fwrite(LED,['SET:LEVEL:CHANNEL1,' powervalue ';SET:LEVEL:CHANNEL2,' powervalue ';'])
pause(0.5)
while LED.BytesAvailable>0
    LEDMessage = fscanf(LED);
end

waitbar(1,waitb,'Updating Positions...');

%If there are more than one filter selected
%Use the grid of optimum z planes to find a best fit for the z position of
%focus
if length(Filters)>1
    for z = find(Filters)
        if flag(z) == 1
            Temp_Plane = [new_planes{1}(z,:);new_planes{2}(z,:);new_planes{3}(z,:);new_planes{4}(z,:);new_planes{5}(z,:)];
            fittp = fit([Temp_Plane(:,1) Temp_Plane(:,2)], Temp_Plane(:,3), 'thinplateinterp');
            
            [m,n,~] = size(POS);
            
            minZ = min(Temp_Plane(:,3))-1000;
            
            NewPOS = POS;
            %Update the z position to match the best fit of optimum z
            %planes
            for ii = 1:m
                for j = 1:n
                    NewPOS(ii,j,3) = fittp(NewPOS(ii,j,1),NewPOS(ii,j,2));
                    %If the new position is less than the minimum z plane
                    %then set it to the minimum z plane
                    if NewPOS(ii,j,3)<minZ
                        NewPOS(ii,j,3) = minZ;
                    end
                end
            end
            POS_array{z} = NewPOS;
        else
            POS_array{z} = [];
        end
    end
  %Autopopulates unselected channels with z positions of selected autofocus channels
  %Because there has to be a position array for every channel
    ind = find(~cellfun('isempty',POS_array),1,'first');
    if isempty(ind)
        NewPOS = POS;
    else
        NewPOS = POS_array{ind};
    end
    if isempty(POS_array{1})
        POS_array{1} = NewPOS;
    end
    if isempty(POS_array{2})
        POS_array{2} = NewPOS;
    end
    if isempty(POS_array{3})
        POS_array{3} = NewPOS;
    end
    if isempty(POS_array{4})
        POS_array{4} = NewPOS;
    end
    if isempty(POS_array{5})
        POS_array{5} = NewPOS;
    end
    if isempty(POS_array{6})
        POS_array{6} = NewPOS;
    end
    fittp_BF = [];
    fittp_CTL = [];
%If there is one filter or less selected --> find best fit, update z position
%values for brightfield
else
    fittp_BF = fit([Plane(:,1) Plane(:,2)], Plane(:,3), 'thinplateinterp');
    [m,n,~] = size(POS);
    
    minZ = min(Plane(:,3))-1000;
    
    NewPOS = POS;
    for ii = 1:m
        for j = 1:n
            NewPOS(ii,j,3) = fittp_BF(NewPOS(ii,j,1),NewPOS(ii,j,2));
            if NewPOS(ii,j,3)<minZ
                NewPOS(ii,j,3) = minZ;
            end
        end
    end
    %If there Filters is empty, then set CTL filter z position to
    %brithtfield z position
    %Used to create a fit for microraft release --> ensures equal distance
    %for releasing microrafts
    if isempty(Filters)
        fittp_CTL = fittp_BF;
    else
        if flag(Filters) == 1
            Temp_Plane = [new_planes{1}(Filters,:);new_planes{2}(Filters,:);new_planes{3}(Filters,:);new_planes{4}(Filters,:);new_planes{5}(Filters,:)];
            fittp_CTL = fit([Temp_Plane(:,1) Temp_Plane(:,2)], Temp_Plane(:,3), 'thinplateinterp');
        else
            fittp_CTL = fittp_BF;
        end
    end
    POS_array = NewPOS;
end

% Plane = [topleft;topright;bottomleft;bottomright;center];
close(waitb)



