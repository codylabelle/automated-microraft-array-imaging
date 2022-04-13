function fm = FScore(Z_Pos)
global mmc Div BB zFScore
mmc.setPosition('ZStage',Z_Pos);
mmc.waitForSystem()
%Uses Snap2 function to take an image with camera
img = Snap2(mmc);

%BB = bounding box; if brightfield makes boundaries of 3x3 box
if isempty(BB)
    fm_temp = zeros(1,Div^2);
    fm_temp2 = zeros(1,Div^2);
    count = 0;
    for ii = 0:(Div-1)
        for jj = 0:(Div-1)
            count = count+1;
            %Takes standard deviation of image (in vector form)
            %could replace fmeasure with std2  *GRAT*, *GLVA*
%             fm_temp(count) = fmeasure(img,'GLVA',[100+(ii*(1848/Div)),100+(jj*(1848/Div)),(1848/Div),(1848/Div)]);
            fm_temp(count) = fmeasure(img,'GLVN',[100+(ii*(1848/Div)),100+(jj*(1848/Div)),(1848/Div),(1848/Div)]);
            fm_temp2(count) = fmeasure(img,'HISR',[100+(ii*(1848/Div)),100+(jj*(1848/Div)),(1848/Div),(1848/Div)]);
        end
    end
else
    %For brightfield, use GLVA --> finds cells the best
    if zFScore == 6
        %For fluorescence focus scores, uses individual microraft positions -->
        %only looks at cells inside of microrafts
        %looks at each raft individually
        fm_temp = zeros(1,length(BB));
        for ii = 1:length(BB)
            fm_temp(ii) = fmeasure(img,'GLVA',BB(ii).BoundingBox);
        end
    %For fluorescence images, use other measures; GRAE, HISR,VOLA, SFRQ, *GRAT*     
    else
        fm_temp = zeros(1,length(BB));
        for ii = 1:length(BB)
            fm_temp(ii) = fmeasure(img,'HISR',BB(ii).BoundingBox);
        end
    end
end
%Finds the negative median focus score from the different 3x3 microraft areas
if zFScore == 6  %Brightfield: only if selected for autofocus
    fm = -median(fm_temp);
elseif isempty(BB) %Brightfield: autofocus prior to fluorescence autofocus
    fm = -median(fm_temp,'all')+median(fm_temp2,'all')*0.1; %0.1 adjusts magnitude of HISR to GLVN
else %Fluorescence autofocus
    fm = -mean(fm_temp);
end
