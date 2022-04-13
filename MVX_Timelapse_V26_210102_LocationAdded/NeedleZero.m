function [FinalZ] = NeedleZero(NeedleExtension)
global mmc sm2
%This function zeroes the needle
%Input: Needle Extension value in microns
%Output: Final extension value in microns

%Focus on the array, then lower the objective ~9.5 mm (for 2x objective) --> moves objective out of the way for zeroing needle device
CurrentZ = mmc.getPosition('ZStage');
%Adding to current z location lowers the objective: correction factor is 0.0432 --> 400/.04323 = 9456 um
MoveZ = CurrentZ + 400;
mmc.setPosition('ZStage',MoveZ);
mmc.waitForSystem
NewZ = mmc.getPosition('ZStage');
pause(1.5);

%Check to see if objective is out of the way (+/- 10 um on Prior display), then zero the needle release device --> retract 3000um (30mm)
NeedleStep = 51;
NeedleExtension = round(NeedleExtension);
if NewZ < (MoveZ+10) && NewZ > (MoveZ-10)
    %Retract (+)
    move(sm2,29988/NeedleStep);
    release(sm2);
    %Next, extend the needle 20mm and set current Z height to 20mm
    %Extend (-)
    move(sm2,-NeedleExtension/NeedleStep);
    release(sm2);
    FinalZ = NeedleExtension;
    %Last, refocus on array
    mmc.setPosition('ZStage',CurrentZ);
    mmc.waitForSystem
    pause(1.5);
else
    f = msgbox('Could not zero needle device due to Z stage')
end