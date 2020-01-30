%% isopycnal mean(d), pman
%contour mean depth of isopycnal
%on pcolor of manifold probability function
load('isopycnalDepths2.mat')
%sigma 26
load('manifoldsIso26int8.mat')
%pman(pman==0)=nan;
figure; 
[c2,h2]=contour(xh,yh,pman,[0.0 0.0],'k');
corex=c2(1,224:246);corey=c2(2,224:246);
clf
contourf(xh,yh,pman,[0:0.1:1]); shading 'flat'; colorbar; colormap(cbrewer('seq','OrRd',9))
hold on
plot(corex,corey,'m','LineWidth',2)
[C,h]=contour(XC,YC,mean(iso26depth,3),[0 5 8 12],'k');
clabel(C,h,'fontsize',22)
%contour(xh,yh,pman,[0 0:0.1:1]);
plot(lonCoast,latCoast,'k','LineWidth',2)
caxis([0 0.9])
title('\sigma=26 manifold pdf and mean depth','fontsize',24)
xlabel('longitude','fontsize',24);ylabel('latitude','fontsize',24)
set(gca,'fontsize',22)
%% sigma 26.5
load('manifoldsIso265int8.mat')
%pman(pman==0)=nan;
figure; 
[c2,h2]=contour(xh,yh,pman,[0.0 0.0],'k');
corex=c2(1,252:296);corey=c2(2,252:296);
clf
contourf(xh,yh,pman,-0.1:0.1:1); colorbar; colormap(cbrewer('seq','OrRd',9))
hold on
plot(corex,corey,'g','LineWidth',2)
[C,h]=contour(XC,YC,mean(iso265depth,3),5:20:150,'k');
plot(lonCoast,latCoast,'k','LineWidth',2)
clabel(C,h,'fontsize',22)
caxis([0 0.9])
title('\sigma=26.5 manifold pdf and mean depth','fontsize',24)
xlabel('longitude','fontsize',24);ylabel('latitude','fontsize',24)
set(gca,'fontsize',22)
%% sigma 27
load('manifoldsIso27int8.mat')
%pman(pman==0)=nan;
figure; 
[c2,h2]=contour(xh,yh,pman,[0.0 0.0],'k');
corex=c2(1,235:285);corey=c2(2,235:285);
clf
contourf(xh,yh,pman,-0.1:0.1:1); shading 'flat'; colorbar; colormap(cbrewer('seq','OrRd',9))
hold on
plot(corex,corey,'c','LineWidth',2);
[C,h]=contour(XC,YC,mean(iso27depth,3),[20:25:150],'k');
clabel(C,h,'fontsize',22)
caxis([0 0.9])
plot(lonCoast,latCoast,'k','LineWidth',2)
title('\sigma=27 manifold pdf and mean depth','fontsize',24)
xlabel('longitude','fontsize',24);ylabel('latitude','fontsize',24)
set(gca,'fontsize',22)
%% sigma 27.5
load('manifoldsIso275int8.mat')
%pman(pman==0)=nan;
figure; 
[c2,h2]=contour(xh,yh,pman,[0.0 0.0],'k');
corex=c2(1,133:197);corey=c2(2,133:197);
clf
contourf(xh,yh,pman,-0.1:0.1:1); shading 'flat'; colorbar; colormap(cbrewer('seq','OrRd',9))
hold on
plot(corex,corey,'b','LineWidth',2);
[C,h]=contour(XC,YC,mean(iso275depth,3),[30:25:250],'k');
clabel(C,h,'fontsize',22)
caxis([0 0.9])
plot(lonCoast,latCoast,'k','LineWidth',2)
title('\sigma=27.5 manifold pdf and mean depth','fontsize',24)
xlabel('longitude','fontsize',24);ylabel('latitude','fontsize',24)
set(gca,'fontsize',22)
%% surface
%
load('analyzedSurfaceMan8day.mat', 'pman')
load('analyzedSurfaceMan8day.mat', 'xh')
load('analyzedSurfaceMan8day.mat', 'yh')

figure; 
[c2,h2]=contour(xh,yh,pman,[0.0 0.0],'k');
corex=c2(1,663:704);corey=c2(2,663:704);
clf
load('manifoldsIso275int8.mat')
load('analyzedSurfaceMan8dayPman2.mat','pman*')
contourf(xh,yh,pman,-0.1:0.1:1); shading 'flat'; colorbar; colormap(cbrewer('seq','OrRd',9))
hold on; plot(corex,corey,'r','LineWidth',2);
plot(lonCoast,latCoast,'k','LineWidth',2)
caxis([0 0.9])
title('surface manifold pdf','fontsize',24)
xlabel('longitude','fontsize',24);ylabel('latitude','fontsize',24)
set(gca,'fontsize',22)
%% domain

figure; contourf(XC,YC,d,0:200:2000); colorbar; shading 'flat'
colormap(cbrewer('seq','BuPu',10))
hold on
plot(XC(1,:),YC(1,:),'g','LineWidth',2)
plot(lonCoast,latCoast,'k','LineWidth',2)
h_ledge=legend('depth','domain','coast','Location','SE');
set(h_ledge,'fontsize',24)

plot([-7.5 -7.5],[37.15 39],'k')
plot([-2 -2],[33 35.05],'k')

plot(XC(700,:),YC(700,:),'g','LineWidth',2)
plot(XC(:,1),YC(:,1),'g','LineWidth',2)
plot(XC(:,200),YC(:,200),'g','LineWidth',2)

xlabel('longitude','fontsize',24);
ylabel('latitude','fontsize',24)
set(gca,'fontsize',22)

%% ssh
load('sshCompareAvisoFull.mat')
mA=mean(sshA(:,:,1:148),3);
mM=mean(sshMc(:,:,1:148),3);
sM=std(sshMc(:,:,1:148),0,3);
sA=std(sshA(:,:,1:148),0,3);

figure; subplot(1,2,1); contourf(long2,latg2,mM,-0.2:0.02:0.2); colormap(b2r(-0.18,0.2))
hold on; plot(lonCoast,latCoast,'k','LineWidth',2)
axis([-5.9 1.5 35 38.5])
ylabel('Mean SSH','fontsize',24)
set(gca,'fontsize',22)
subplot(1,2,2); contourf(long2,latg2,mA,-0.2:0.02:0.2); colormap(b2r(-0.18,0.2)); colorbar
hold on; plot(lonCoast,latCoast,'k','LineWidth',2)
set(gca,'fontsize',22)
axis([-5.9 1.5 35 38.5])

figure;
subplot(1,2,1); contourf(long2,latg2,sM,0:0.01:1); colormap(cbrewer('seq','Reds',10))
hold on; plot(lonCoast,latCoast,'k','LineWidth',2)
ylabel('std SSH','fontsize',24)
axis([-5.9 1.5 35 38.5])
xlabel('Model','fontsize',24)
set(gca,'fontsize',22)
subplot(1,2,2); contourf(long2,latg2,sA,0:0.01:1); colormap(cbrewer('seq','Reds',10)); colorbar
hold on; plot(lonCoast,latCoast,'k','LineWidth',2)
set(gca,'fontsize',22)
xlabel('AVISO','fontsize',24)
axis([-5.9 1.5 35 38.5])

%% center
load('manifoldsIso275int8.mat')


figure(10)
[c2,h2]=contour(xh,yh,pman,[0.0 0.0],'k');
corex=c2(1,133:197);corey=c2(2,133:197);
clf
figure(11)
plot(corex,corey,'b','LineWidth',2)
hold on

load('manifoldsIso27int8.mat')
figure(10)
[c2,h2]=contour(xh,yh,pman,[0.0 0.0],'k');
corex=c2(1,235:285);corey=c2(2,235:285);
clf
figure(11)
plot(corex,corey,'c','LineWidth',2)

load('manifoldsIso265int8.mat')
figure(10)
[c2,h2]=contour(xh,yh,pman,[0.0 0.0],'k');
corex=c2(1,252:296);corey=c2(2,252:296);
clf
figure(11)
plot(corex,corey,'g','LineWidth',2)

load('manifoldsIso26int8.mat')
figure(10)
[c2,h2]=contour(xh,yh,pman,[0.0 0.0],'k');
corex=c2(1,224:246);corey=c2(2,224:246);
clf
figure(11)
plot(corex,corey,'m','LineWidth',2)

load('analyzedSurfaceMan8day.mat', 'xh')
load('analyzedSurfaceMan8day.mat', 'yh')
load('analyzedSurfaceMan8day.mat', 'pman')
figure(10)
[c2,h2]=contour(xh,yh,pman,[0.0 0.0],'k');
corex=c2(1,663:704);corey=c2(2,663:704);
clf
figure(11)
plot(corex,corey,'r','LineWidth',2)
plot(lonCoast,latCoast,'k','LineWidth',8)
legend('center \sigma=27.5','center \sigma=27','center \sigma=26.5','center \sigma=26','surface center','coast')
axis([-5.5 -2.5 35 37])
set(gca,'fontsize',22)
xlabel('longitude','fontsize',24);
ylabel('latitude','fontsize',24)
%% center 3D
load('isopycnalDepths2.mat')
load('manifoldsIso275int8.mat')
figure(10)
[c2,h2]=contour(xh,yh,pman,[0.0 0.0],'k');
corex=c2(1,133:197);corey=c2(2,133:197);
clf
corez=griddata(XC,YC,nanmean(iso275depth,3),corex,corey);
corez(corez<125)=125;
figure(11)
surf(XC(200:550,:),YC(200:550,:),-150-d(200:550,:)/150); shading 'flat'
colormap(cbrewer('seq','BuPu',50)); colormap(flipud(colormap))
shading 'flat'
hold on; 
plot3(corex,corey,-corez,'b','LineWidth',3)


load('manifoldsIso27int8.mat')
figure(10)
[c2,h2]=contour(xh,yh,pman,[0.0 0.0],'k');
corex=c2(1,235:285);corey=c2(2,235:285);
clf
corez=griddata(XC,YC,nanmean(iso27depth,3),corex,corey);
figure(11)
plot3(corex,corey,-corez,'c','LineWidth',3)

load('manifoldsIso265int8.mat')
figure(10)
[c2,h2]=contour(xh,yh,pman,[0.0 0.0],'k');
corex=c2(1,252:296);corey=c2(2,252:296);
clf
corez=griddata(XC,YC,nanmean(iso265depth,3),corex,corey);
figure(11)
plot3(corex,corey,-corez,'g','LineWidth',3)

load('manifoldsIso26int8.mat')
figure(10)
[c2,h2]=contour(xh,yh,pman,[0.0 0.0],'k');
corex=c2(1,224:246);corey=c2(2,224:246);
clf
corez=griddata(XC,YC,nanmean(iso26depth,3),corex,corey);
figure(11)
plot3(corex,corey,-corez,'m','LineWidth',2)

load('analyzedSurfaceMan8day.mat', 'xh')
load('analyzedSurfaceMan8day.mat', 'yh')
load('analyzedSurfaceMan8day.mat', 'pman')
figure(10)
[c2,h2]=contour(xh,yh,pman,[0.0 0.0],'k');
corex=c2(1,663:704);corey=c2(2,663:704);
corez=zeros(size(corex));
clf
figure(11)
plot3(corex,corey,corez,'r','LineWidth',3)
%plot(lonCoast,latCoast,'k','LineWidth',8)
contour3(XC(200:550,:),YC(200:550,:),-150-d(200:550,:)/150,[-150 -150],'k')
legend('topography','center \sigma=27.5','center \sigma=27','center \sigma=26.5','center \sigma=26','surface center','coast')
axis tight
set(gca,'fontsize',22)
xlabel('longitude','fontsize',24);
ylabel('latitude','fontsize',24)
%% center 14 day manifolds
load('manifoldsIso275int14.mat')
figure(10)
[c2,h2]=contour(xh,yh,pman,[0.0 0.0],'k');
corex=c2(1,[207:243 207]);corey=c2(2,[207:243 207]);
clf
figure(12)
plot(corex,corey,'b','LineWidth',2)
hold on

load('manifoldsIso27int14.mat')
figure(10)
[c2,h2]=contour(xh,yh,pman,[0.0 0.0],'k');
corex=c2(1,[207:238 207]);corey=c2(2,[207:238 207]);
clf
figure(12)
plot(corex,corey,'c','LineWidth',2)

load('manifoldsIso265int14.mat')
figure(10)
[c2,h2]=contour(xh,yh,pman,[0.0 0.0],'k');
corex=c2(1,[217:232 217]);corey=c2(2,[217:232 217]);
clf
figure(12)
plot(corex,corey,'g','LineWidth',2)

load('manifoldsIso26int14.mat')
figure(10)
[c2,h2]=contour(xh,yh,pman,[0.0 0.0],'k');
corex=c2(1,[204:212 204]);corey=c2(2,[204:212 204]);
clf
figure(12)
plot(corex,corey,'m','LineWidth',2)

save('analyzedSurfaceMan14dayPman2.mat','pman*')
figure(10)
[c2,h2]=contour(xh,yh,pman,[0.0 0.0],'k');
corex=-4.15;corey=35.8;
clf
figure(12)
plot(corex,corey,'ro','LineWidth',2,'MarkerSize',10)
plot(lonCoast,latCoast,'k','LineWidth',8)
legend('center \sigma=27.5','center \sigma=27','center \sigma=26.5','center \sigma=26','surface center','coast')
axis([-5.5 -2.5 35 37])
set(gca,'fontsize',22)
xlabel('longitude','fontsize',24);
ylabel('latitude','fontsize',24)

%% center 14 day manifolds, 3D

load('isopycnalDepths2.mat')
corez=griddata(XC,YC,nanmean(iso275depth,3),corex,corey);
corez(corez<125)=125;
figure(13)
surf(XC(200:550,:),YC(200:550,:),-150-d(200:550,:)/150); shading 'flat'
colormap(cbrewer('seq','BuPu',50)); colormap(flipud(colormap))
shading 'flat'
hold on; 

load('manifoldsIso275int14.mat')
figure(10)
[c2,h2]=contour(xh,yh,pman,[0.0 0.0],'k');
corex=c2(1,[207:243 207]);corey=c2(2,[207:243 207]);
corez=griddata(XC,YC,nanmean(iso275depth,3),corex,corey);
clf
figure(13)
plot3(corex,corey,-corez,'b','LineWidth',3)
hold on

load('manifoldsIso27int14.mat')
figure(10)
[c2,h2]=contour(xh,yh,pman,[0.0 0.0],'k');
corex=c2(1,[207:238 207]);corey=c2(2,[207:238 207]);
corez=griddata(XC,YC,nanmean(iso27depth,3),corex,corey);
clf
figure(13)
plot3(corex,corey,-corez,'c','LineWidth',3)

load('manifoldsIso265int14.mat')
figure(10)
[c2,h2]=contour(xh,yh,pman,[0.0 0.0],'k');
corex=c2(1,[217:232 217]);corey=c2(2,[217:232 217]);
corez=griddata(XC,YC,nanmean(iso265depth,3),corex,corey);
clf
figure(13)
plot3(corex,corey,-corez,'g','LineWidth',3)

load('manifoldsIso26int14.mat')
figure(10)
[c2,h2]=contour(xh,yh,pman,[0.0 0.0],'k');
corex=c2(1,[204:212 204]);corey=c2(2,[204:212 204]);
corez=griddata(XC,YC,nanmean(iso26depth,3),corex,corey);
clf
figure(13)
plot3(corex,corey,-corez,'m','LineWidth',3)

save('analyzedSurfaceMan14dayPman2.mat','pman*')
figure(10)
[c2,h2]=contour(xh,yh,pman,[0.0 0.0],'k');
corex=-4.15;corey=35.8;
clf
figure(13)
plot3(corex,corey,0,'ro','LineWidth',2,'MarkerSize',10)
contour3(XC(200:550,:),YC(200:550,:),-150-d(200:550,:)/150,[-150 -150],'k')
legend('topography','center \sigma=27.5','center \sigma=27','center \sigma=26.5','center \sigma=26','surface center','coast')
axis tight
set(gca,'fontsize',22)
xlabel('longitude','fontsize',24);
ylabel('latitude','fontsize',24)
view(-14,34)
%% center changing times
load('manifoldsIso265int14.mat')
figure(10)
[c2,h2]=contour(xh,yh,pman,[0.0 0.0],'k');
corex=c2(1,217:end); corey=c2(2,217:end);
clf
figure(13)
plot(corex,corey,'g','LineWidth',2)
hold on

load('manifoldsIso26int10.mat','xh','yh','pman')
figure(10)
[c2,h2]=contour(xh,yh,pman,[0.0 0.0],'k');
corex=c2(1,198:end); corey=c2(2,198:end);
clf
figure(13)
plot(corex,corey,'m','LineWidth',2)

load('analyzedSurfaceMan8day.mat', 'xh')
load('analyzedSurfaceMan8day.mat', 'yh')
load('analyzedSurfaceMan8day.mat', 'pman')
figure(10)
[c2,h2]=contour(xh,yh,pman,[0.0 0.0],'k');
corex=c2(1,663:704);corey=c2(2,663:704);
clf
figure(13)
plot(corex,corey,'r','LineWidth',2)
plot(lonCoast,latCoast,'k','LineWidth',8)
title('WAG Center same-length manifolds')
legend('\sigma=26.5, 14d','\sigma=26, 10d','surface 8d')

%% 2 lobes
%finish redoing!!!
load('lobePropertiesFluxesSurfaceIso26.mat')
load('manifoldsIso265int8.mat')
i=21;
j=1;
nn=find(ylobeH(i-8,j,:)>0,1,'last');
xlobe=squeeze(xlobeH(i-8,j,1:nn));
ylobe=squeeze(ylobeH(i-8,j,1:nn));

%ERROR: should be using l-tr-all(21,...), not 13!
%but lobes not closed in 8day sigma 26.5...could use 14day?
figure; plot([lontrBall(13,58:182) lontrFall(13,108:130)],[lattrBall(13,58:182) lattrFall(13,108:130)],'c','LineWidth',1.5);
hold on
plot([lontrFall(24-8,360:-1:290) lontrBall(16,254:318)],[lattrFall(24-8,360:-1:290) lattrBall(16,254:318)],'m','LineWidth',1.5)
plot(xlobe,ylobe,'b','LineWidth',2)
% i=22;
% nn=find(ylobeH(i-8,j,:)>0,1,'last');
% xlobe=squeeze(xlobeH(i-8,j,1:nn));
% ylobe=squeeze(ylobeH(i-8,j,1:nn));
% hold all; plot(xlobe,ylobe,'LineWidth',2)
i=24;
nn=find(ylobeH(i-8,j,:)>0,1,'last');
xlobe=squeeze(xlobeH(i-8,j,1:nn));
ylobe=squeeze(ylobeH(i-8,j,1:nn));
hold all; plot(xlobe,ylobe,'r','LineWidth',2)
plot(lonCoast,latCoast,'k','LineWidth',2)
set(gca,'fontsize',22)
axis([-6.2 -2.4 35 37])
legend('day 21, \sigma=26.5','day 24, \sigma=26.5','day 21, surface','day24, surface')
title('Lobe 6','fontsize',26)

i=80;
j=1;
nn=find(ylobeH(i-8,j,:)>0,1,'last');
xlobe=squeeze(xlobeH(i-8,j,1:nn));
ylobe=squeeze(ylobeH(i-8,j,1:nn));
figure; plot(xlobe,ylobe,'LineWidth',2)
i=83;
nn=find(ylobeH(i-8,j,:)>0,1,'last');
xlobe=squeeze(xlobeH(i-8,j,1:nn));
ylobe=squeeze(ylobeH(i-8,j,1:nn));
hold all; plot(xlobe,ylobe,'LineWidth',2)
% i=63;
% nn=find(ylobeH(i-8,j,:)>0,1,'last');
% xlobe=squeeze(xlobeH(i-8,j,1:nn));
% ylobe=squeeze(ylobeH(i-8,j,1:nn));
% hold all; plot(xlobe,ylobe,'LineWidth',2)
plot(lonCoast,latCoast,'k','LineWidth',2)
set(gca,'fontsize',22)
axis([-6 -2 35 37])
legend('day 80','day 83')
title('Lobe 20','fontsize',24)

%% circulation sketched
load('uvwNativeGridIsoDepth285.mat')
load('uvwAve148.mat')

xmin=min(min(XC));
ymin=min(min(YC));
xcm=(XC-xmin).*111000.*cosd(YC);
ycm=(YC-ymin).*111000;
lonCM=(lonCoast-xmin).*111000.*cosd(latCoast);
latCM=(latCoast-ymin).*111000;

x=0:2000:max(xcm(:,1)); y=0:2000:max(ycm(end,:));
[xg,yg]=meshgrid(x,y);

u1=griddata(xcm,ycm,Urot(:,:,1),xg,yg);
v1=griddata(xcm,ycm,Vrot(:,:,1),xg,yg);
u2=griddata(xcm,ycm,mean(uIso,3),xg,yg);
v2=griddata(xcm,ycm,mean(vIso,3),xg,yg);

xs1=2e5:20000:8e5;
ys1=2e5:20000:4.5e5;
[xs,ys]=meshgrid(xs1,ys1);

% figure;
% h2=streamline(stream2(xg,yg,u2,v2,xs,ys,[0.1 120]));
% set(h2,'Color','g')
% hold on
% quiver(xg(1:22:end,1:12:end),yg(1:22:end,1:12:end),u2(1:22:end,1:12:end),v2(1:22:end,1:12:end),'Color',[0 1 0.5],'LineWidth',2)
% plot(lonCM,latCM,'k','LineWidth',2)
% 
% figure
% h1=streamline(stream2(xg,yg,u1,v1,xs,ys,[0.1 130]));
% set(h1,'Color','b')
% hold on
% quiver(xg(5:25:end,1:15:end),yg(5:25:end,1:15:end),u1(5:25:end,1:15:end),v1(5:25:end,1:15:end),'Color',[0.5 0 1],'LineWidth',2)
% plot(lonCM,latCM,'k','LineWidth',2)


figure;
h2=streamline(stream2(xg,yg,u2,v2,xs,ys,[0.1 150]));
set(h2,'Color','b','LineWidth',1.5)
hold on
quiver(xg(1:35:end,3:25:end),yg(1:35:end,3:25:end),u2(1:35:end,3:25:end),v2(1:35:end,3:25:end),'Color',[0 0.25 1],'LineWidth',2)
%quiver(xg(163,166),yg(163,166),u2(163,166),v2(163,166),'Color',[0 0.25 1],'LineWidth',2)
h1=streamline(stream2(xg,yg,u1,v1,xs,ys,[0.1 200]));
set(h1,'Color','r','LineWidth',1.5)
hold on
%quiver(xg(5:40:end,1:30:end),yg(5:40:end,1:30:end),u1(5:40:end,1:30:end),v1(5:40:end,1:30:end),'Color',[1 0.25 0],'LineWidth',2)
plot(lonCM,latCM,'k','LineWidth',2)
axis([2e5 8e5 2e5 4.5e5])

%% fluxes
load('lobePropertiesFluxesSurfaceIso26.mat')
fVolLobe2=fVolLobe; fVolLobe2(fVolLobe2==0)=nan;
fMassLobe2=fMassLobe; fMassLobe2(fMassLobe2==0)=nan;
fSaltLobe2=fSaltLobe; fSaltLobe2(fSaltLobe2==0)=nan;
fHeatLobe2=fHeatLobe; fHeatLobe2(fHeatLobe2==0)=nan;

figure
plot(1:148,fVolLobe2/1e6,'bo','MarkerSize',12,'LineWidth',2)
hold on
plot(1:148,fMassLobe2/1e9,'rp','MarkerFaceColor','r','MarkerSize',10)
plot(1:148,zeros(148,1),'k')
legend('Volume, 10^6 m^3 s^{-1}','Mass, 10^9 kg s^{-1}')
title('Volume and Mass Fluxes, Sv','fontsize',24)
set(gca,'fontsize',22)
xlabel('Simulation Day','fontsize',24)

figure
plot(1:148,fSaltLobe2/1e6,'bo','MarkerSize',12,'LineWidth',2,'LineStyle','none')
hold on
plot(1:148,fHeatLobe2/1e12,'rV','MarkerFaceColor','r','MarkerSize',10,'LineStyle','none')
plot(1:148,zeros(148,1),'k')
legend('Salt, 10^6 kg/s','Heat, 10^{12} J/s')
title('Salt and Heat Fluxes','fontsize',24)
set(gca,'fontsize',22)
xlabel('Simulation Day','fontsize',24)


%% mean 148 day velocity sections
load('uvwAve148.mat', 'Urot')
load('uvwAve148.mat', 'Vrot')
load('uvwAve148.mat', 'XC')
load('uvwAve148.mat', 'YC')
load('uvwAve148.mat', 'd')
load('uvwAve148.mat', 'dBin')

figure; plot(lonCoast,latCoast,'k')
hold on; plot(XC(:,110),YC(:,110),'b')
hold on; plot(XC(225,:),YC(225,:),'r')
hold on; plot(XC(350,:),YC(350,:),'m')
title('Cross-section locations')

figure; contourf(XC(:,110),-dBin,squeeze(Urot(:,110,:)).')
caxis([-1 1])
axis([-9 1 -2000 0])
hold on; plot(XC(:,110),-d(:,110),'k','LineWidth',2)
title('mean U, blue curve')

figure; contourf(XC(:,110),-dBin,squeeze(Vrot(:,110,:)).')
caxis([-1 1])
axis([-9 1 -2000 0])
hold on; plot(XC(:,110),-d(:,110),'k','LineWidth',2)
title('mean V, blue curve')

figure; contourf(YC(225,:),-dBin,squeeze(Urot(225,:,:)).')
caxis([-1 1])
hold on; plot(YC(225,:),-d(225,:),'k','LineWidth',2)
axis([35.7 36.2 -1000 0])
title('mean U, red curve')

figure; contourf(YC(225,:),-dBin,squeeze(Vrot(225,:,:)).')
caxis([-1 1])
hold on; plot(YC(225,:),-d(225,:),'k','LineWidth',2)
axis([35.7 36.2 -1000 0])
title('mean V, red curve')


figure; contourf(YC(350,:),-dBin,squeeze(Urot(350,:,:)).'); colorbar
caxis([-1 1])
hold on; plot(YC(350,:),-d(350,:),'k','LineWidth',2)
axis([35 37 -2000 0])
title('mean U, magenta curve')

figure; contourf(YC(350,:),-dBin,squeeze(Vrot(350,:,:)).'); colorbar
caxis([-1 1])
hold on; plot(YC(350,:),-d(350,:),'k','LineWidth',2)
axis([35 37 -2000 0])
title('mean V, magenta curve')