function MoogMove(x,y,s)
%Conversion factor from step size to mm; 3um/step
%Requires integer values of steps
X = x/3;
Y = y/3;
%Set position of indivdual motors
fprintf(s,['PT:1=',num2str(Y),' G:1 ']);
fprintf(s,['PT:2=',num2str(X),' G:2 ']);

% 
% %Set position of indivdual motors
% fprintf(s,'PT:1=%d G:1 ',y);
% fprintf(s,'PT:2=%d G:2 ',x);
