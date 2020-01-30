%cross-isopycnal flux estimate from 3D trajectories
load('iso275depthNFrev.mat')
load('geometrySpinupSteady')
load('traj3Diso27sInt1day32.mat')

%%
xmin=min(XC(:)); ymin=min(YC(:));
xcM=111000*cosd(YC).*(XC-xmin*ones(size(XC)));
ycM=111000*(YC-ymin.*ones(size(YC)));
xuM=111000*cosd(YU).*(XU-xmin*ones(size(XU)));
yuM=111000*(YU-ymin.*ones(size(YU)));
xvM=111000*cosd(YV).*(XV-xmin*ones(size(XV)));
yvM=111000*(YV-ymin.*ones(size(YV)));
%%
x0=xf1(1,:);
y0=yf1(1,:);
z0=zf1(1,:);
figure; plot3(x0,y0,z0,'o'); hold all; plot3(xf1,yf1,zf1)
%% vertical error
zsigma=griddata(xcM,ycM,-isoDepth(:,:,33),xf,yf);
zErr=zf-zsigma;
figure; pcolor(xcM,ycM,vertVelCross(:,:,31)); shading 'flat'; hold on; scatter(xf,yf,36,zErr./86400,'linewidth',2)