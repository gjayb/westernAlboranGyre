%mean 148 day plots

load('geometrySpinupSteady.mat','XC','YC','dInterface')
dBin=0.5*(dInterface(2:end))+0.5*dInterface(1:end-1);
load('saTCtCpSigma162NR.mat','Sigma')
SigmaM=nanmean(Sigma(:,:,:,1:148),4);
load('uvwDailyNativeNF.mat','U')
UmeanN=nanmean(U(:,:,:,1:148),4);
load('uvwDailyNativeNF.mat','V')
VmeanN=nanmean(V(:,:,:,1:148),4);
load('uvwDailyNativeNF.mat','W')
WmeanN=nanmean(W(:,:,:,1:148),4);
load('tsPress162NF.mat','PotTave')
T=nanmean(PotTave(:,:,:,1:148),4); clear PotTave
load('tsPress162NF.mat','PractSave')
S=nanmean(PractSave(:,:,:,1:148),4); clear PractSave
save('uvwtsSigma148meanNF.mat')

%%
% load('uvwSigma148meanNF.mat')
% 
% figure; pcolor(YC(225,:),-dBin,squeeze(Sigma(225,:,:)).'); shading 'flat'; colorbar
% axis([35.8 36.1 -1000 0])
% caxis([24 30])
% saveas(gcf,'meanSigmaX225nf','fig')
% 
% figure; pcolor(YC(350,:),-dBin,squeeze(Sigma(350,:,:)).'); shading 'flat'; colorbar
% axis([35.1 36.7 -1000 0])
% caxis([24 30])
% saveas(gcf,'meanSigmaX350nf','fig')
% 
% figure; pcolor(YC(225,:),-dBin,squeeze(T(225,:,:)).'); shading 'flat'; colorbar
% hold on; [c1,h1]=contour(YC(225,:),-dBin,squeeze(Sigma(225,:,:)).',[26.3 26.5 27 27.5 28 28.5 29 29.1],'k');
% clabel(c1,h1)
% axis([35.8 36.1 -1000 0])
% caxis([2 20])
% saveas(gcf,'meanSigmaX225nf','fig')
% 
% figure; pcolor(YC(350,:),-dBin,squeeze(T(350,:,:)).'); shading 'flat'; colorbar
% axis([35.1 36.7 -1000 0])
% caxis([2 20])
% saveas(gcf,'meanSigmaX350nf','fig')
