function step_video(mtimer,~)
global mmc
d=get(mtimer,'UserData');
img=mmc.getLastImage();
img=typecast(img,'uint16');
img=reshape(img,[2048,2048]);
img=rot90(img,2);
d{2}=img;
if d{3}
    %'Saving'
    t=d{1};
    if mod(t,5)==0 %4 to 32 are good values
        imwrite(rot90(img,1),[d{4} '\ReleaseData\ReleaseVideo\' num2str(t) '.tif'])
    end
end
d{1}=d{1}+1;
set(mtimer,'UserData',d)
end
