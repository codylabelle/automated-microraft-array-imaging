function [dx,dy]=modifyRetarget(littlei,th,i,count,idlist,displthresh)
%littlei= cropped intensity image
%th=intensity threshold value
%dxs, dys= list of all dx and dys taken (for debugging purposes)
%i = current raft number
%count=ith tagret's cumulative number of release attempts
%idlist = displacement adjustment parameter from beginning of NeedleRelease_V2
% displthrsh = maximum allowable displacement, in pixels, from beginning of NeedleRelease_V2

lm = littlei<th(i); %threshold image intensity
s = regionprops(lm,'Area','Centroid'); %check to make sure there are objects
if ~isempty(s) % if there are objects, find the area-weighted center of mass of all segmented objects
    [y, x] = ndgrid(1:size(lm, 1), 1:size(lm, 2));
    c = median([x(logical(lm)), y(logical(lm))]);
    try
        dx = c(1)-(size(littlei,1)/2);%positive, to right
        dy = c(2)-(size(littlei,1)/2);%positive, down
    catch
        dx = 0;
        dy = 0;
    end
else %not possible but w/e
    dx = 0;dy = 0;
end
%% restrict re-aiming to prevent overshoots
if abs(dx)>3
    dx=dx*idlist(count(i));
end
if abs(dy)>3
    dy=dy*idlist(count(i));
end
%% remove restrictions for later releases BUT add a max limit
if (abs(dx)>10 &  count(i)>6) | abs(dx)>displthresh  %Original valus: dx>10, count >5
    dx=sign(dx)*displthresh;
end
if (abs(dy)>10 & count(i)>6) | abs(dy)>displthresh   %Original valus: dx>10, count >5
    dy=sign(dy)*displthresh;
end
