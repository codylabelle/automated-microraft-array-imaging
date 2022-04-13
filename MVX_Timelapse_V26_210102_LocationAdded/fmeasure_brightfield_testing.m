global mmc LED BB
% Use this code to test autofocus in fluorescence channels --> must be used
% in conjuction with the invertedMVXGUI

%%%%%%%%%%% Change these parameters
% pixel_size = 1.0257;
pixel_size = 0.8203;
M = 8    %6.4;
powerLED = 120;

%%%%%%%%%%%
Array_Dim=[100 50 140];
Width=100;
Gap=50;
Div = floor(2048/(Width+Gap)*pixel_size/3);

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


NA =0.5;
%equation for depth of focus for a microscope
df = (((Lambda*n)/(NA^2))+((n*e)/(M*NA)))*FudgeFactor;

%First set positions and create bounding boxes
value = powerLED;
powervalue=num2str(value);
fwrite(LED,['SET:LEVEL:CHANNEL1,' powervalue ';SET:LEVEL:CHANNEL2,' powervalue ';'])
pause(0.5)
while LED.BytesAvailable>0
    LEDMessage = fscanf(LED);
end

X = mmc.getXPosition;
Y = mmc.getYPosition;
Z = mmc.getPosition('ZStage');

mmc.waitForSystem()
img = Snap2(mmc);
bw = uraft(img,Width,Gap,pixel_size);
seBW = strel('disk',round((Width/4)/pixel_size));
bw = imerode(bw,seBW);
% imshowpair(img,bw2)
BB = regionprops(bw,'BoundingBox');

value = 0;
powervalue=num2str(value);
fwrite(LED,['SET:LEVEL:CHANNEL1,' powervalue ';SET:LEVEL:CHANNEL2,' powervalue ';'])
pause(0.5)
while LED.BytesAvailable>0
    LEDMessage = fscanf(LED);
end
% imtool(bw)
imshowpair(img,bw)

%%
tic

value = powerLED;
powervalue=num2str(value);
fwrite(LED,['SET:LEVEL:CHANNEL1,' powervalue ';SET:LEVEL:CHANNEL2,' powervalue ';'])
pause(0.5)
while LED.BytesAvailable>0
    LEDMessage = fscanf(LED);
end

%Second: switch to fluorescence channel and take z stack images
distance_limit_lower = 40*df;
distance_limit_upper = 40*df;
ZPOS = Z-(distance_limit_lower):10*FudgeFactor:Z+(distance_limit_upper);
img = cell(size(ZPOS));
img_d = cell(size(ZPOS));

for i = 1:length(ZPOS)
    mmc.setPosition('ZStage',ZPOS(i))
    mmc.waitForSystem()
    img{i} = Snap2(mmc);
    img_d{i} = double(img{i});
end
toc

value = 0;
powervalue=num2str(value);
fwrite(LED,['SET:LEVEL:CHANNEL1,' powervalue ';SET:LEVEL:CHANNEL2,' powervalue ';'])
pause(0.5)
while LED.BytesAvailable>0
    LEDMessage = fscanf(LED);
end


%%
% run fmeasure on brightfield
BB = [];
% measure = {'ACMO','BREN','CURV','GDER','GLVA','GLLV','GLVN','GRAE','GRAT','GRAS','HELM','HISE','HISR','LAPE','LAPM','LAPV','LAPD','SFIL','SFRQ','TENG','TENV','VOLA'};
measure = {'LAPM','LAPE','GLVA','LAPV','LAPD','GRAT','GRAE','VOLA','GLVN','HISR'};
Time = zeros(length(measure),1);
fm = zeros(length(img),length(measure));
for i = 1:length(measure)
    i
    tic
    for j = 1:length(img)
        fm_temp = zeros(Div,Div);
        count = 0;
        for ii = 0:(Div-1)
            for jj = 0:(Div-1)
                count = count+1;
                %Takes standard deviation of image (in vector form)
                %could replace fmeasure with std2
                fm_temp(ii+1,jj+1) = fmeasure(img{j},measure{i},[100+(ii*(1848/Div)),100+(jj*(1848/Div)),(1848/Div),(1848/Div)]);
            end
        end
        fm(j,i) = -median(fm_temp,'all');
%         fm(j,i) = -mean2(fm_temp);
    end
    Time(i) = toc;
end

lowerZ = ceil(length(ZPOS)/2-30);
upperZ = ceil(length(ZPOS)/2+30);

%Find the location of the minimum f measurement
[~,minfm] = min(fm);
%Find which measures are greater than the minimum and maximum position
IND = find(minfm<=upperZ & minfm>=lowerZ);
[~,minTime] = min(Time(IND));
measure{IND(minTime)}

[TimeSort,minTime] = sort(Time(IND));
measure(IND(minTime))
TimeSort'
minfm(IND(minTime))
%Notes: only displays measures for which the best score was between the
%limits
%Compare the values given (ex. 22 or 26) to the positions in ZPOS to see
%the z location 

%%
%Ideal = 22;
Ideal = round(length(ZPOS)/2);
count2 = 1;
for i = IND
    plot(fm(:,i))
    hold on
    plot(Ideal,fm(Ideal,i),'-ro')
    hold on
    plot(minfm(count2),fm(minfm(count2),i),'b*')
    hold off
    title([measure{i},', Time: ',num2str(Time(i))])
    waitforbuttonpress
    count2 = count2 +1;
end

%%
% run fmeasure on brightfield with two measures
BB = [];
measure = {'GLVN','HISR'};
Time = zeros(length(measure),1);
fm = zeros(length(img),length(measure));
fm = []; fm_temp = []; fm_temp2=[];
for j = 1:length(img)
    fm_temp = zeros(Div,Div);
    count = 0;
    for ii = 0:(Div-1)
        for jj = 0:(Div-1)
            count = count+1;
            fm_temp(ii+1,jj+1) = fmeasure(img{j},measure{1},[100+(ii*(1848/Div)),100+(jj*(1848/Div)),(1848/Div),(1848/Div)]);
            fm_temp2(ii+1,jj+1) = fmeasure(img{j},measure{2},[100+(ii*(1848/Div)),100+(jj*(1848/Div)),(1848/Div),(1848/Div)]);
        end
    end
%     meanFM1 = mean2(fm_temp); meanFM2 = mean2(fm_temp2);  
%     ratioFM = meanFM1/meanFM2; %Don't use moving ratio: creates badcurve
    fm(j) = -median(fm_temp,'all')+median(fm_temp2,'all')*0.1;
end

[~,minfm] = min(fm)
plot([1:length(img)],fm(1,:));
hold on
plot(minfm,fm(minfm),'r*')
display(num2str(ZPOS(minfm)))
display(num2str(ZPOS(minfm)+2))

