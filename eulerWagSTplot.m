%figures for advective fluxes Euler WAG
load('distancesAreas','RAC')
load('geometrySpinupSteady','*C*')
load('salinityBudgetEulerWag148daysNF275S.mat', 'inWag')
load('uvwtsSigma148meanNF.mat', 'S')
load('fluxesAdvS148day.mat', 'AdvWs')
load('edgesWAGeuler2017NF.mat')
for i=1:148; AdvWs(:,:,:,i)=AdvWs(:,:,:,i).*openZd.*repmat(RAC,[1 1 46]); end
mAdvWs=mean(AdvWs,4);  clear AdvWs
inWagH3=repmat(inWag(:,:,1),[1 1 46]);
mAdvWs(inWagH3==0)=nan;


openEW=zeros([700 200 46]); openEW(:,:,1:10)=openE(:,:,1:10)+openW(:,:,1:10); openEW(openEW>1)=1;
load('fluxesAdvS148day.mat', 'AdvUs')
load('temperatureBudgetEulerWAG148dayNF275S.mat', 'cellVol')
load('temperatureBudgetEulerWAG148dayNF275S.mat', 'DYC')
load('temperatureBudgetEulerWAG148dayNF275S.mat', 'DXC')
for i=1:148; AdvUs(:,:,:,i)=AdvUs(:,:,:,i).*openEW.*cellVol./repmat(DXC,[1 1 46]); end
mAdvUs=mean(AdvUs,4); clear AdvUs;
 

%%
openNS=zeros([700 200 46]); openNS(:,:,1:10)=openN(:,:,1:10)+openS(:,:,1:10); openNS(openNS>1)=1;
load('fluxesAdvS148day.mat', 'AdvVs')
for i=1:148; AdvVs(:,:,:,i)=AdvVs(:,:,:,i).*openNS.*cellVol./repmat(DYC,[1 1 46]); end
mAdvVs=mean(AdvVs,4); clear AdvVs;


%%
holdvar=sum(mAdvUs,3); 
index1=find(holdvar~=0);
holdvar2=sum(mAdvVs,3); index2=find(holdvar2~=0); 
scaleS=5e10;
figure; pcolor(XC,YC,sum(mAdvWs,3)); shading 'flat'
caxis([-2e10 2e10])
hold on; quiver([-5; XC(index1(1:3:end))],[36.4; YC(index1(1:3:end))],[1e10; holdvar(index1(1:3:end))]/scaleS,zeros(size([1; index1(1:3:end)])),'linewidth',2,'autoscale','off')
quiver([-5; XC(index2(1:3:end))],[36.4; YC(index2(1:3:end))],zeros(size([1; index2(1:3:end)])),[1e10; holdvar2(index2(1:3:end))]/scaleS,'linewidth',2,'autoscale','off')
title('Euler WAG mean vertically integrated salinity fluxes, psu m^3/s','fontsize',14)
axis([-5.5 -2.7 35.1 36.6])
set(gca,'fontsize',12)
%%
load('uvwtsSigma148meanNF.mat', 'T')
load('fluxesAdvT148day.mat', 'AdvWt')
for i=1:148; AdvWt(:,:,:,i)=AdvWt(:,:,:,i).*openZd.*repmat(RAC,[1 1 46]); end
mAdvWt=mean(AdvWt,4); clear AdvWt
mAdvWt(inWagH3==0)=nan; 

load('fluxesAdvT148day.mat', 'AdvUt')

for i=1:148; AdvUt(:,:,:,i)=AdvUt(:,:,:,i).*openEW.*cellVol./repmat(DXC,[1 1 46]); end
mAdvUt=mean(AdvUt,4); clear AdvUt;
 
load('fluxesAdvT148day.mat', 'AdvVt')
for i=1:148; AdvVt(:,:,:,i)=AdvVt(:,:,:,i).*openNS.*cellVol./repmat(DYC,[1 1 46]); end
mAdvVt=mean(AdvVt,4); clear AdvVt;


%%
holdvar=sum(mAdvUt,3); 
index1=find(holdvar~=0);
holdvar2=sum(mAdvVt,3); index2=find(holdvar2~=0); 
scaleT=5e10;
figure; pcolor(XC,YC,sum(mAdvWt,3)); shading 'flat'
caxis([-1e10 1e10])
hold on; quiver([-5; XC(index1(1:3:end))],[36.4; YC(index1(1:3:end))],[1e10; holdvar(index1(1:3:end))]/scaleT,zeros(size([1; index1(1:3:end)])),'linewidth',2,'autoscale','off')
quiver([-5; XC(index2(1:3:end))],[36.4; YC(index2(1:3:end))],zeros(size([1; index2(1:3:end)])),[1e10; holdvar2(index2(1:3:end))]/scaleT,'linewidth',2,'autoscale','off')
title('Euler WAG mean vertically integrated potential temperature fluxes, \circ C m^3/s','fontsize',14)
axis([-5.5 -2.7 35.1 36.6])
set(gca,'fontsize',12)
%%
holdvar3=sum(mAdvWs,3);
holdvar4=sum(mAdvWt,3);
figure; pcolor(XC,YC,holdvar4./holdvar3); shading 'flat'; colorbar
%%
load('volumesfluxesWAGeuler2017NF275S.mat', 'vfluxhEs','vfluxhWs','vfluxhNs','vfluxhSs','inWag')
mvfluxhXint=mean(sum(vfluxhEs+vfluxhWs,3),4); clear vfluxhEs vfluxhWs
mvfluxhYint=mean(sum(vfluxhNs+vfluxhSs,3),4); clear vfluxhNs vfluxhSs
load('volumesfluxesWAGeuler2017NF275S.mat', 'vfluxvD')
mvfluxD=mean(vfluxvD,4); clear vfluxvD;
%vfluxvint=squeeze(sum(vfluxvD,3)); clear vfluxvD
%%
index1=find(mvfluxhXint~=0 &mvfluxhYint~=0);
%index2=find();
inWagH1=repmat(inWag(:,:,1),[1 1 16]);
mvfluxD(inWagH1==0)=nan;
mvfluxvint=sum(mvfluxD,3);
%%
scaleV=3e5;
mvfluxvint(mvfluxvint==0)=nan;
figure; pcolor(XC,YC,mvfluxvint); shading 'flat'
axis([-6 -1 35 37.5])
caxis([-500 500])
hold on
quiver([-5;XC(index1(1:4:end))],[36.2;YC(index1(1:4:end))],[1e5;mvfluxhXint(index1(1:4:end))]/scaleV,[1e5;mvfluxhYint(index1(1:4:end))]/scaleV,'linewidth',2,'autoscale','off')
% hold all
% quiver([-5;XC(index1(1:3:end))],[36.4;YC(index1(1:3:end))],zeros(size([1;mvfluxhYint(index1(1:3:end))])),[1e5;mvfluxhYint(index1(1:3:end))]/scaleV,'linewidth',2,'autoscale','off')
%plot(lonCoast,latCoast,'k','linewidth',2)
title('Euler WAG mean vertically integrated volume fluxes','fontsize',14)
set(gca,'fontsize',12)
axis([-5.5 -2.7 35.1 36.6])
%%
scaleV=5e5;
mvfluxvint(mvfluxvint==0)=nan;
figure; pcolor(XC,YC,mvfluxvint); shading 'flat'
axis([-6 -1 35 37.5])
caxis([-500 500])
hold on
quiver([-5;XC(index1(1:3:end))],[36.4;YC(index1(1:3:end))],[1e5;mvfluxhXint(index1(1:3:end))]/scaleV,zeros(size([1;mvfluxhXint(index1(1:3:end))])),'linewidth',2,'autoscale','off')
hold all
quiver([-5;XC(index1(1:3:end))],[36.4;YC(index1(1:3:end))],zeros(size([1;mvfluxhYint(index1(1:3:end))])),[1e5;mvfluxhYint(index1(1:3:end))]/scaleV,'linewidth',2,'autoscale','off')
%plot(lonCoast,latCoast,'k','linewidth',2)
title('Euler WAG mean vertically integrated volume fluxes','fontsize',14)
set(gca,'fontsize',12)
axis([-5.5 -2.7 35.1 36.6])
%%
load('uvwtsSigma148meanNF.mat', 'S')
figure; c1=contour(XC,YC,S(:,:,1).*double(XC<-3).*double(XC>-5.3),[36.475 36.475]);
xedge=c1(1,2:898);
yedge=c1(2,2:898);
xedge2=xedge(265:650);
yedge2=yedge(265:650);

dy=diff(yedge2); dx=diff(xedge2);
dxF=TriScatteredInterp(xedge2(1:end-1).',yedge2(1:end-1).',dx.','nearest');
dxXC=dxF(XC(index1),YC(index1));
dyF=TriScatteredInterp(xedge2(1:end-1).',yedge2(1:end-1).',dy.','nearest');
dyXC=dyF(XC(index1),YC(index1));
mags=sqrt(dxXC.^2+dyXC.^2);
vecnorm=[-dyXC.'./mags.';dxXC.'./mags.'];
magnorm=dot(vecnorm,[mvfluxhXint(index1).';mvfluxhYint(index1).']);
xnorm=magnorm.*-dyXC.'./mags.';
ynorm=magnorm.*dxXC.'./mags.';

scaleV=1e4;
mvfluxvint(mvfluxvint==0)=nan;
figure; pcolor(XC,YC,mvfluxvint); shading 'flat'
axis([-6 -1 35 37.5])
caxis([-500 500])
hold on
quiver([-5;XC(index1(1:4:end))],[36.2;YC(index1(1:4:end))],[scaleV/sqrt(2);xnorm(1:4:end).']/scaleV,[scaleV/sqrt(2);ynorm(1:4:end).']/scaleV,'linewidth',2)


