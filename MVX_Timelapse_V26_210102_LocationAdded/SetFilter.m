function [] = SetFilter(mmc,Filter)

%0: Bright Field
%1: UV EX / Blue EM
%2: Violet EX / Blue EM
%3: Blue EX / Green EM
%4: Green EX / Red EM
%5: Red EX / Far Red EM

%C = mmc.getCameraDevice();
%S = mmc.getShutterDevice();

switch Filter
    case 0
        mmc.setShutterOpen('Shutter',0)
        mmc.setAutoShutter(0)
        mmc.setProperty('Wheel_EM','State',0)
        mmc.setProperty('Wheel_EX','State',0)
    case 1 %UV/Blue
        mmc.setAutoShutter(1);
        mmc.setProperty('Wheel_EX','State',0)
        mmc.setProperty('Wheel_EM','State',1)
    case 2 %Violet/Blue
        mmc.setAutoShutter(1);
        mmc.setProperty('Wheel_EX','State',1)
        mmc.setProperty('Wheel_EM','State',1)
    case 3  %Blue/Green
        mmc.setAutoShutter(1);
        mmc.setProperty('Wheel_EX','State',2)
        mmc.setProperty('Wheel_EM','State',2)
    case 4  %Green/Red
        mmc.setAutoShutter(1);
        mmc.setProperty('Wheel_EX','State',3)
        mmc.setProperty('Wheel_EM','State',3)
    case 5 %Red/Far red
        mmc.setAutoShutter(1);
        mmc.setProperty('Wheel_EX','State',4)
        mmc.setProperty('Wheel_EM','State',4)
    otherwise
        mmc.setShutterOpen('Shutter',0)
        mmc.setAutoShutter(0)
        mmc.setProperty('Wheel_EM','State',0)
end

mmc.waitForSystem()
