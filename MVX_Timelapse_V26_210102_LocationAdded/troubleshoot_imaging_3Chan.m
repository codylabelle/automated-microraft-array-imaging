global mmc

Filters = [1 0 1 1];
Exposures = [100 200 150];

for k = [find(Filters)]
    SetFilter(mmc,k)
    mmc.setExposure(Exposures(k))
    mmc.waitForSystem()
    [img{k}, camera_error_count] = Snap2(mmc);
end

%reset
SetFilter(mmc,0)

%%
imtool(img{1},[])
imtool(img{2},[])
imtool(img{3},[])

%% testing a pause

global mmc

Filters = [1 0 1 1];
Exposures = [100 200 150];

%Flag changes to 1 after first image taken; only pause after first
flag = 0;
for k = [find(Filters)]
    SetFilter(mmc,k)
    mmc.setExposure(Exposures(k))
    mmc.waitForSystem()
    
    if length(k) > 2 && flag > 0
        pause(0.2)
    end
    
    [img{k}, camera_error_count] = Snap2(mmc);
    flag = flag + 1;
end

%reset
SetFilter(mmc,0)

%%%%% Remember that flag count and length needed to be adjusted to acount
%%%%% for brightfield image: use length(k) > 3 and flag > 1


