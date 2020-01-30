%mean 148 day plots

% load('geometrySpinupSteady.mat','XC','YC','dInterface')
% dBin=0.5*(dInterface(2:end))+0.5*dInterface(1:end-1);
% load('saTCtCpSigma162NR.mat','Sigma')
% SigmaM=nanmean(Sigma(:,:,:,1:148),4);
% load('uvwDailyNativeNF.mat','U')
% UmeanN=nanmean(U(:,:,:,1:148),4);
% load('uvwDailyNativeNF.mat','V')
% VmeanN=nanmean(V(:,:,:,1:148),4);
% load('uvwDailyNativeNF.mat','W')
% WmeanN=nanmean(W(:,:,:,1:148),4);
% load('tsPress162NF.mat','PotTave')
% T=nanmean(PotTave(:,:,:,1:148),4); clear PotTave
% load('tsPress162NF.mat','PractSave')
% S=nanmean(PractSave(:,:,:,1:148),4); clear PractSave
% save('uvwtsSigma148meanNF.mat')

%%
load('uvwtsSigma148meanNF.mat')

figure; pcolor(YC(225,:),-dBin,squeeze(Sigma(225,:,:)).'); shading 'flat'; colorbar
axis([35.8 36.1 -1000 0])
caxis([24 30])
%saveas(gcf,'meanSigmaX225nf','fig')

figure; pcolor(YC(350,:),-dBin,squeeze(Sigma(350,:,:)).'); shading 'flat'; colorbar
axis([35.1 36.7 -1600 0])
caxis([24 30])

%%
figure; pcolor(YC(225,:),-dBin,squeeze(T(225,:,:)).'); shading 'flat'; colorbar
hold on; [c1,h1]=contour(YC(225,:),-dBin,squeeze(Sigma(225,:,:)).',[26.3 26.5 26.75 27 27.5 28 28.5 28.9 29 29.1],'k');
clabel(c1,h1)
axis([35.8 36.1 -1000 0])
caxis([6 20]); title('Mean Potentail Temperature, Strait')
xlabel('latitude')

figure; pcolor(YC(350,:),-dBin,squeeze(T(350,:,:)).'); shading 'flat'; colorbar
hold on; [c1,h1]=contour(YC(350,:),-dBin,squeeze(Sigma(350,:,:)).',[26.3 26.5 26.75 27 27.5 28 28.5 28.9 29 29.1],'k');
clabel(c1,h1); axis([35.1 36.7 -1600 0])
caxis([6 20]); title('Mean Potential Temperature, WAG')
xlabel('latitude')

figure; pcolor(XC(:,120),-dBin,squeeze(T(:,120,:)).'); shading 'flat'; colorbar
hold on; [c1,h1]=contour(XC(:,120),-dBin,squeeze(Sigma(:,120,:)).',[26.3 26.5 26.75 27 27.5 28 28.5 28.9 29 29.1],'k');
clabel(c1,h1); axis([-8 0 -1600 0])
caxis([6 20]); title('Mean Potential Temperature and Potential Density')
xlabel('longitude')

figure; pcolor(YC(225,:),-dBin,squeeze(S(225,:,:)).'); shading 'flat'; colorbar
hold on; [c1,h1]=contour(YC(225,:),-dBin,squeeze(Sigma(225,:,:)).',[26.3 26.5 26.75 27 27.5 28 28.5 28.9 29 29.1],'k');
clabel(c1,h1)
axis([35.8 36.1 -1000 0])
caxis([32 40]); title('Mean Salinity, Strait')
xlabel('latitude')

figure; pcolor(YC(350,:),-dBin,squeeze(S(350,:,:)).'); shading 'flat'; colorbar
hold on; [c1,h1]=contour(YC(350,:),-dBin,squeeze(Sigma(350,:,:)).',[26.3 26.5 26.75 27 27.5 28 28.5 28.9 29 29.1],'k');
clabel(c1,h1); axis([35.1 36.7 -1600 0])
caxis([32 40]); title('Mean Salinity, WAG')
xlabel('latitude')

figure; pcolor(XC(:,120),-dBin,squeeze(S(:,120,:)).'); shading 'flat'; colorbar
hold on; [c1,h1]=contour(XC(:,120),-dBin,squeeze(Sigma(:,120,:)).',[26.3 26.5 26.75 27 27.5 28 28.5 28.9 29 29.1],'k');
clabel(c1,h1); axis([-8 0 -1600 0])
caxis([32 40]); title('Mean Salinity and Potential Density')
xlabel('longitude')
%%
figure; pcolor(YC(225,:),-dBin,86400*squeeze(WmeanN(225,:,:)).'); shading 'flat'; colorbar
hold on; [c1,h1]=contour(YC(225,:),-dBin,squeeze(Sigma(225,:,:)).',[26.3 26.5 26.75 27 27.5 28 28.5 28.9 29 29.1],'k');
clabel(c1,h1)
axis([35.8 36.1 -1000 0])
colormap(b2r(-150,150)); 
title('Mean Vertical Velocity in m/day, Strait')
xlabel('latitude')

figure; pcolor(YC(350,:),-dBin,86400*squeeze(WmeanN(350,:,:)).'); shading 'flat'; colorbar
hold all; [c1,h1]=contour(YC(350,:),-dBin,squeeze(Sigma(350,:,:)).',[26.3 26.5 26.75 27 27.5 28 28.5 28.9 29 29.1],'k');
clabel(c1,h1); 
axis([35.1 36.7 -1600 0])
colormap(b2r(-20,20)); 
title('Mean Vertical Velocity in m/day, WAG')
xlabel('latitude')

figure; pcolor(XC(:,120),-dBin,86400*squeeze(WmeanN(:,120,:)).'); shading 'flat'; colorbar
hold all; [c1,h1]=contour(XC(:,120),-dBin,squeeze(Sigma(:,120,:)).',[26.3 26.5 26.75 27 27.5 28 28.5 28.9 29 29.1],'k');
clabel(c1,h1); axis([-8 0 -1600 0])
colormap(b2r(-50,50)); 
title('Mean Vertical Velocity (m/day) and Potential Density')
xlabel('longitude')
%%
load('geometrySpinupSteady','Angle*')
Urot=UmeanN.*repmat(AngleCS,[1 1 46]) - VmeanN.*repmat(AngleSN,[1 1 46]);  
Vrot=UmeanN.*repmat(AngleSN,[1 1 46]) + VmeanN.*repmat(AngleCS,[1 1 46]); 

figure; pcolor(YC(225,:),-dBin,squeeze(UmeanN(225,:,:)).'); shading 'flat'; colorbar
hold on; [c1,h1]=contour(YC(225,:),-dBin,squeeze(Sigma(225,:,:)).',[26.3 26.5 26.75 27 27.5 28 28.5 28.9 29 29.1],'k');
clabel(c1,h1)
axis([35.8 36.1 -1000 0])
caxis([-1 1]); colormap(b2r(-1,1))
title('Mean Zonal Velocity and Potential Density, Strait')

figure; pcolor(YC(350,:),-dBin,squeeze(UmeanN(350,:,:)).'); shading 'flat'; colorbar
hold on; [c1,h1]=contour(YC(350,:),-dBin,squeeze(Sigma(350,:,:)).',[26.3 26.5 26.75 27 27.5 28 28.5 28.9 29 29.1],'k');
clabel(c1,h1)
axis([35.1 36.7 -1600 0]); colormap(b2r(-1,1))
caxis([-1 1]); colormap(b2r(-1,1))
title('Mean Zonal Velocity and Potential Density, WAG')

figure; pcolor(YC(225,:),-dBin,squeeze(UmeanN(225,:,:)).'); shading 'flat'; colorbar
hold on; [c1,h1]=contour(YC(225,:),-dBin,squeeze(Sigma(225,:,:)).',[26.3 26.5 26.75 27 27.5 28 28.5 28.9 29 29.1],'k');
clabel(c1,h1)
axis([35.8 36.1 -1000 0])
caxis([-1 1]); colormap(b2r(-1,1))
title('Mean Eastward Velocity and Potential Density, Strait')

figure; pcolor(YC(350,:),-dBin,squeeze(UmeanN(350,:,:)).'); shading 'flat'; colorbar
hold on; [c1,h1]=contour(YC(350,:),-dBin,squeeze(Sigma(350,:,:)).',[26.3 26.5 26.75 27 27.5 28 28.5 28.9 29 29.1],'k');
clabel(c1,h1)
axis([35.1 36.7 -1600 0]); colormap(b2r(-1,1))
caxis([-1 1]); colormap(b2r(-1,1))
title('Mean Eastward Velocity and Potential Density, WAG')

figure; pcolor(XC(:,120),-dBin,squeeze(VmeanN(:,120,:)).'); shading 'flat'; colorbar
hold on; [c1,h1]=contour(XC(:,120),-dBin,squeeze(Sigma(:,120,:)).',[26.3 26.5 26.75 27 27.5 28 28.5 28.9 29 29.1],'k');
clabel(c1,h1)
 colormap(b2r(-1,1))
caxis([-1 1]); colormap(b2r(-1,1))
title('Mean Meridional Velocity and Potential Density')

figure; pcolor(XC(:,120),-dBin,squeeze(Vrot(:,120,:)).'); shading 'flat'; colorbar
hold on; [c1,h1]=contour(XC(:,120),-dBin,squeeze(Sigma(:,120,:)).',[26.3 26.5 26.75 27 27.5 28 28.5 28.9 29 29.1],'k');
clabel(c1,h1)
 colormap(b2r(-1,1))
caxis([-1 1]); colormap(b2r(-1,1))
title('Mean Northward Velocity and Potential Density')

figure; pcolor(XC(:,120),-dBin,squeeze(Urot(:,120,:)).'); shading 'flat'; colorbar
hold on; [c1,h1]=contour(XC(:,120),-dBin,squeeze(Sigma(:,120,:)).',[26.3 26.5 26.75 27 27.5 28 28.5 28.9 29 29.1],'k');
clabel(c1,h1)
 colormap(b2r(-1,1))
caxis([-1 1]); colormap(b2r(-1,1))
title('Mean Eastward Velocity and Potential Density')
