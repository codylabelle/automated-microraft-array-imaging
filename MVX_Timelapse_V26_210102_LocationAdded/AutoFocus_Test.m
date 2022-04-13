%%% Use this code to test autofocus settings in brightfield or fluorescence

%Find focal plane for BF
global BB mmc LED zFScore Div  

Div = floor(2048/(Width+Gap)*pixel_size/3);
NA = 0.5;

Lambda = .45;
e = 6.5;
n = 1.33;
%Fudge factor for z stage movement; correction factor for depth of well
%FudgeFactor = .023;
FudgeFactor = .04323;
Well_Depth = 175;
correction = Well_Depth*0.64*FudgeFactor*2;
%If testing BF input 1, if fluorescence input 0;
BF = 0;
powerLED = 120; % set power level of brightfield

M = 8;
%equation for depth of focus for a microscope
df = (((Lambda*n)/(NA^2))+((n*e)/(M*NA)))*FudgeFactor;

powervalue=num2str(powerLED);
fwrite(LED,['SET:LEVEL:CHANNEL1,' powervalue ';SET:LEVEL:CHANNEL2,' powervalue ';'])
pause(0.5)
while LED.BytesAvailable>0
    LEDMessage = fscanf(LED);
end
X = mmc.getXPosition;
Y = mmc.getYPosition;
Z = mmc.getPosition('ZStage');
%Uses optimization function to find the optimum z plane for focusing
opt = optimset('MaxFunEval',20,'Display','off','TolX',2*df);
%Experimentally found minimums and maximums of z plane
%Note that moving 4 df's doesn't change much, so 40*df is about 10*df
%experimentally

Z_lower = Z-20*df;  %40
Z_upper = Z+20*df;
BB = [];
mmc.waitForSystem()
Flag_BF = 0;
%While an actual minimum has not been found, Flag_BF is 0 -->
%ensure that you find the true minimum
tic
while Flag_BF < 1
    [ZF,~,Flag_BF] = myfminbnd(@FScore,Z_lower,Z_upper,opt);
    %         Z_lower = Z_lower-20*df;
    %         Z_upper = Z_upper+40*df;
    Z_lower = Z_lower-30*df;  %30
    Z_upper = Z_upper+30*df;   %30 
    opt.MaxFunEvals = opt.MaxFunEvals*2;
end
display(['Normal ZF: ',num2str(ZF)]); display(['Adjusted ZF: ',num2str(ZF+2)]);
toc
%Turn off LED
powervalue=num2str(0);
fwrite(LED,['SET:LEVEL:CHANNEL1,' powervalue ';SET:LEVEL:CHANNEL2,' powervalue ';'])
pause(0.5)
while LED.BytesAvailable>0
    LEDMessage = fscanf(LED);
end

%% Use this step to create bounding box; change after if in fluorescence
powervalue=num2str(powerLED);
fwrite(LED,['SET:LEVEL:CHANNEL1,' powervalue ';SET:LEVEL:CHANNEL2,' powervalue ';'])
pause(0.5)
while LED.BytesAvailable>0
    LEDMessage = fscanf(LED);
end
zfOffset = 2;
ZF = ZF + zfOffset;
mmc.setPosition('ZStage',ZF)
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

if BF == 0
    powervalue=num2str(0);
    fwrite(LED,['SET:LEVEL:CHANNEL1,' powervalue ';SET:LEVEL:CHANNEL2,' powervalue ';'])
    pause(0.5)
    while LED.BytesAvailable>0
        LEDMessage = fscanf(LED);
    end
end
imshowpair(img,bw)
%%
if BF == 1
    powervalue=num2str(300);
    fwrite(LED,['SET:LEVEL:CHANNEL1,' powervalue ';SET:LEVEL:CHANNEL2,' powervalue ';'])
    pause(0.5)
    while LED.BytesAvailable>0
        LEDMessage = fscanf(LED);
    end
end

%Change 'z' to determine which measure to use in FScore
if BF == 1
    zFScore=6; %6 for brightfield
else
    zFScore=1; %1 is abritrary
end
%Find focus score for brigthfield inside of rafts
tic
%Uses optimization function to find the optimum z plane for focusing
opt = optimset('MaxFunEval',40,'Display','off','TolX',1*df);
Z_lower = ZF-20*df;
Z_upper = ZF+40*df;
%Find optimum z plane for focusing in each channel selected

%Optimization of z plane within the segemented rafts
[ZF2,~,flag] = myfminbnd(@FScore,Z_lower,Z_upper,opt);
if flag<1
    opt = optimset('MaxFunEval',40,'Display','off','TolX',df);
    Z_lower = ZF-30*df;
    Z_upper = ZF+40*df;
    [ZF2,~,flag] = myfminbnd(@FScore,Z_lower,Z_upper,opt);
    %Reset 
    opt = optimset('MaxFunEval',40,'Display','off','TolX',df);
    Z_lower = ZF-20*df;
    Z_upper = ZF+40*df;
end
toc

%Turn off LED
powervalue=num2str(0);
fwrite(LED,['SET:LEVEL:CHANNEL1,' powervalue ';SET:LEVEL:CHANNEL2,' powervalue ';'])
pause(0.5)
while LED.BytesAvailable>0
    LEDMessage = fscanf(LED);
end

ZF2