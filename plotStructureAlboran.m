%plots an example of the structures in the alboran sea via manifolds of
%hyperbolic points outlining the WAG and EAG

load('C:\Users\JayB\Documents\MATLAB\MITgcm\hyperbolic10daySurface.mat')
lon1=lontr;
lat1=lattr;
lon1b=lontr2;
lat1b=lattr2;

load('C:\Users\JayB\Documents\MATLAB\MITgcm\hyperbolicEAGeast10daySurface.mat')
lon2=lontr;
lat2=lattr;
lon2b=lontr2;
lat2b=lattr2;

load('C:\Users\JayB\Documents\MATLAB\MITgcm\hyperbolicEAG10daySurface.mat')
lon3=lontr;
lat3=lattr;
lon3b=lontr2;
lat3b=lattr2;

load('C:\Users\JayB\Documents\MATLAB\MITgcm\hyperbolicWAG10daySurface.mat')
lon4=lontr;
lat4=lattr;
lon4b=lontr2;
lat4b=lattr2;

ploting4=zeros(20);
ploting4(1:5:end)=1;
ploting4(1:8,:)=0;
ploting4=logical(ploting4);
ploting4b=ploting4;
ploting4b([1:6 15:20],:)=0;


load('C:\Users\JayB\Documents\MATLAB\MITgcm\hyperbolicStrait10daySurface.mat')
lon5=lontr;
lat5=lattr;
lon5b=lontr2;
lat5b=lattr2;

ploting5=zeros(20);
ploting5(512,1:10)=1;
ploting5=logical(ploting5);


figure
plot(lonCoast,latCoast,'k')
hold on
plot(lon1(:,[7 90:150]),lat1(:,[7 90:150]),'r',lon2(:,[110:20:300 240 241]),lat2(:,[110:20:300 240 241]),'r',lon3(:,34:20:360),lat3(:,34:20:360),'r',lon4(:,ploting4(:)),lat4(:,ploting4(:)),'r',lon5(:,ploting5(:)),lat5(:,ploting5(:)),'r')
plot(lon1b(:,[7 90:150]),lat1b(:,[7 90:150]),'b',lon2b(:,[110:20:300 240 241]),lat2b(:,[110:20:300 240 241]),'b',lon3b(:,16:20:366),lat3b(:,16:20:366),'b',lon4b(:,ploting4b(:)),lat4b(:,ploting4b(:)),'b',lon5b(:,ploting5(:)),lat5b(:,ploting5(:)),'b')

figure
plot(lonCoast,latCoast,'k')
hold on
scatter(lon1(:),lat1(:),16,'ro')
scatter(lon2(:),lat2(:),16,'ro')
scatter(lon3(:),lat3(:),16,'ro')
scatter(lon4(:),lat4(:),16,'ro')
scatter(lon5(:),lat5(:),16,'ro')
plot(lon1b,lat1b,'b',lon2b,lat2b,'b',lon3b,lat3b,'b',lon4b,lat4b,'b',lon5b,lat5b,'b')