addpath('../mStuff')
for Ni=57:64
fn1=strcat('wagBoundary30m9dDay',num2str(Ni),'F.mat');
fn2=strcat('wagBoundary30m9dDay',num2str(Ni),'B.mat');
load(fn1);
load(fn2);
disp('load done')
disp(num2str(Ni))
figure; plot(lontrF(:,1:3:end),lattrF(:,1:3:end),'r')
hold on; plot(lontrB(:,3:3:end),lattrB(:,3:3:end),'b')
plot(lonCoast,latCoast,'k')
axis([-6 0 35 37])
title(strcat('9 days forward and backward starting day ',num2str(Ni)))
fn3=strcat('wagBoundaryPlot30mStartday',num2str(Ni),'.pdf');
disp('saving')
save2pdf(fn3)
close all
clear
end
