function [ZF] = AFocusIsolation(Z, Width, Gap, pixel_size, MagFocus) 
global mmc Div

Div = floor(2048/(Width+Gap)*pixel_size/3);
NA = 0.5;
Lambda = .45;
e = 6.5;
n = 1.33;
%Fudge factor for z stage movement; correction factor for depth of well
%FudgeFactor = .023;
FudgeFactor = .04323;
Well_Depth = 175;
%%correction = Well_Depth*0.64*FudgeFactor*2;

M = MagFocus;
%equation for depth of focus for a microscope
df = (((Lambda*n)/(NA^2))+((n*e)/(M*NA)))*FudgeFactor;

%Uses optimization function to find the optimum z plane for focusing
opt = optimset('MaxFunEval',20,'Display','off','TolX',2*df);
%Experimentally found minimums and maximums of z plane
%Note that moving 4 df's doesn't change much, so 40*df is about 10*df
%experimentally

Z_lower = Z-40*df;
Z_upper = Z+40*df;

mmc.waitForSystem()
Flag_BF = 0;
%While an actual minimum has not been found, Flag_BF is 0 -->
%ensure that you find the true minimum
while Flag_BF < 1
    [ZF,~,Flag_BF] = myfminbnd(@FScoreIsolation,Z_lower,Z_upper,opt);
    %         Z_lower = Z_lower-20*df;
    %         Z_upper = Z_upper+40*df;
    Z_lower = Z_lower-30*df;
    Z_upper = Z_upper+30*df;
    opt.MaxFunEvals = opt.MaxFunEvals*2;
end