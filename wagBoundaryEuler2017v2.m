% %Eulerian WAG boundaries for first 148 days of the simulation
% %Nov through most of Mar
% %Boundaries by limits in time-averaged T,S
% %Fluxes from time-varying U,V
% 
% %addpath('../mStuff')
% %addpath('/nobackup1/gbrett/addforce')
% load('steady148tswRho.mat');
% 
% 
% t2d=Temp(:,:,1)>18.39;
% s2d=Sal(:,:,1)<36.48;
% figure; [c2,h2]=contour(XC,YC,t2d+s2d,[2 2]);
% px=c2(1,870:1376); py=c2(2,870:1376);
% ts2d=inpolygon(XC,YC,px,py);
% ts3d=repmat(ts2d,[1 1 46]);
% sigma3d=Sigma<27.5;
% 
% water=zeros([700 200 46]);
% for i=1:46
%     water(:,:,i)=d>dInterface(i);
% end
% 
% inWag=sigma3d.*ts3d.*water;
% 
% clear W CT Rho
% %%
% % inWag=wagS.*wagT; clear wagS wagT
% 
% % 
% openW=zeros(size(inWag));
% openE=zeros(size(inWag));
% openS=zeros(size(inWag));
% openN=zeros(size(inWag));
% isedge=zeros(size(inWag));
% 
% for zi=1:15
%     for xi=240:450
%         for yi=1:200
%             if inWag(xi,yi,zi)==1
%                 if inWag(xi-1,yi,zi)==0
%                     openW(xi,yi,zi)=1;
%                     isedge(xi,yi,zi)=1;
%                     %disp('W')
%                 end
%                 if inWag(xi+1,yi,zi)==0
%                     openE(xi,yi,zi)=1;
%                     isedge(xi,yi,zi)=1;
%                     %disp('E')
%                 end
%                 if inWag(xi,yi-1,zi)==0
%                     openS(xi,yi,zi)=1;
%                     isedge(xi,yi,zi)=1;
%                     %disp('S')
%                 end
%                 if inWag(xi,yi+1,zi)==0
%                     openN(xi,yi,zi)=1;
%                     isedge(xi,yi,zi)=1;
%                     %disp('N')
%                 end
%             end
%         end
%     end
% end
% % 
% % 
% % figure
% % plot(lonCoast,latCoast,'k')
% % hold all
% % for i=1:14
% % plot(XC(logical(isedge(:,:,i))),YC(logical(isedge(:,:,i))),'.')
% % end
% % %dBin(1:14)
% % legend('Coast','2.5m','8m','14m','20.5m','27.5m','35m','43.5m','53m','63.5m','75m','87.5m','101.5m','117m','134m')
% % title('WAG Edge with depth, Defined by T-S contours')
% % axis ([-6 -1.5 34.75 37.25])
% % %%
% %vertical differences to get vertical fluxes
% openZd=zeros(size(inWag));
% for i=1:15
%    openZd(:,:,i)=(inWag(:,:,i)==1 & inWag(:,:,i+1)==0); 
% end
% 
% openZu=zeros(size(inWag));
% for i=2:16
%    openZu(:,:,i)=(inWag(:,:,i)==1 & inWag(:,:,i-1)==0); 
% end
% save('edgesWAGeuler2017.mat','-v7.3')
% disp('done section 1')
% %% load and initialize for volume flux
% load('edgesWAGeuler2017.mat')
% load('uva148levels15.mat','U');
% load('uva148levels15.mat','V');
% load('wa148levels16.mat','W');
% load('eraDataAndFall2007meansTry2.mat', 'EminusPms')
% load('eraDataAndFall2007meansTry2.mat', 'hoursSince1900')
% load('eraDataAndFall2007meansTry2.mat', 'hflux2')
% load('ssh148.mat', 'ssh148')
% load('distancesAreas.mat','DXG')
% load('distancesAreas.mat','DYG')
% load('distancesAreas.mat','RAC')
% load('distancesAreas.mat','hFacW')
% load('distancesAreas.mat','hFacS')
% %xcM=(XC-min(min(XC))*ones(size(XC))).*111000.*cosd(YC); 
% %ycM=(YC-min(min(YC))*ones(size(YC))).*111000;
% %xuM=(XU-min(min(XU))*ones(size(XU))).*111000.*cosd(YU); 
% %yvM=(YV-min(min(YV))*ones(size(YV))).*111000;
% %[distx,~]=distance(XU(1:end-1,:),YU(1:end-1,:),XU(2:end,:),YU(2:end,:),[6371000 0]);
% %[disty,~]=distance(XV(:,1:end-1),YV(:,1:end-1),XV(:,2:end),YV(:,2:end),[6371000 0]);
% distx=DXG(1:end-1,:);%distance(XG(1:end-1,:),YG(1:end-1,:),XG(2:end,:),YG(2:end,:),[6371000 0]);
% disty=DYG(:,1:end-1);%distance(XG(:,1:end-1),YG(:,1:end-1),XG(:,2:end),YG(:,2:end),[6371000 0]);
% areasRect=RAC(1:end-1,1:end-1);%(0.5.*distx(:,1:end-1)+0.5.*distx(:,2:end)).*(0.5.*disty(1:end-1,:)+0.5.*disty(2:end,:));
% 
% %horizontal volume flux
% vfluxhE=zeros([700 200 16 148]);
% vfluxhW=zeros([700 200 16 148]);
% vfluxhN=zeros([700 200 16 148]);
% vfluxhS=zeros([700 200 16 148]);
% vfluxvD=zeros([700 200 16 148]);
% vfluxvU=zeros([700 200 16 148]);
% 
% %% do calculations for volume flux, save
% disp('start zi=2:15')
% for zi=15%2:15
%       vfluxhW(1:end-1,1:end-1,zi,:)=(dInterface(zi+1)-dInterface(zi)).*repmat((disty(1:end-1,:))...
%           .*openW(1:end-1,1:end-1,zi).*hFacW(1:end-1,1:end-1,zi),[1 1 1 148]).*U(1:end-1,1:end-1,zi,:);
%       vfluxhE(1:end-1,1:end-1,zi,:)=(dInterface(zi+1)-dInterface(zi)).*repmat((disty(2:end,:))...
%           .*openE(1:end-1,1:end-1,zi).*hFacW(2:end,1:end-1,zi),[1 1 1 148]).*U(2:end,1:end-1,zi,:);
%       vfluxhS(1:end-1,1:end-1,zi,:)=(dInterface(zi+1)-dInterface(zi)).*repmat((distx(:,1:end-1))...
%           .*openS(1:end-1,1:end-1,zi).*hFacS(1:end-1,1:end-1,zi),[1 1 1 148]).*V(1:end-1,1:end-1,zi,:);
%       vfluxhN(1:end-1,1:end-1,zi,:)=(dInterface(zi+1)-dInterface(zi)).*repmat((distx(:,2:end))...
%           .*openN(1:end-1,1:end-1,zi).*hFacS(1:end-1,2:end,zi),[1 1 1 148]).*V(1:end-1,2:end,zi,:);
%       vfluxvD(1:end-1,1:end-1,zi,:)=W(1:end-1,1:end-1,zi+1,:).*...
%           repmat(openZd(1:end-1,1:end-1,zi).*...
%         areasRect,[1 1 1 148]);
%       vfluxvU(1:end-1,1:end-1,zi,:)=W(1:end-1,1:end-1,zi,:).*...
%           repmat(openZu(1:end-1,1:end-1,zi).*...
%         areasRect,[1 1 1 148]);
% end
% disp('done zi=2:15')
% max(max(max(max(abs(vfluxhW)))))
% max(max(max(max(abs(vfluxhE)))))
% max(max(max(max(abs(vfluxhN)))))
% max(max(max(max(abs(vfluxhS)))))
% max(max(max(max(abs(vfluxvD)))))
% max(max(max(max(abs(vfluxvU)))))
% zi=1
% vfluxhW(1:end-1,1:end-1,zi,:)=repmat((disty(1:end-1,:))...
%           .*squeeze(openW(1:end-1,1:end-1,zi)).*hFacW(1:end-1,1:end-1,zi),[1 1 148]).*squeeze(U(1:end-1,1:end-1,zi,:)).*(ssh148(1:end-1,1:end-1,:)+5);
% vfluxhE(1:end-1,1:end-1,zi,:)=repmat((disty(2:end,:))...
%           .*squeeze(openE(1:end-1,1:end-1,zi)).*hFacW(2:end,1:end-1,zi),[1 1 148]).*squeeze(U(2:end,1:end-1,zi,:)).*(ssh148(1:end-1,1:end-1,:)+5);
% vfluxhS(1:end-1,1:end-1,zi,:)=repmat((distx(:,1:end-1))...
%           .*squeeze(openS(1:end-1,1:end-1,zi)).*hFacS(1:end-1,1:end-1,zi),[1 1 148]).*squeeze(V(1:end-1,1:end-1,zi,:)).*(ssh148(1:end-1,1:end-1,:)+5);
% vfluxhN(1:end-1,1:end-1,zi,:)=repmat((distx(:,2:end))...
%           .*squeeze(openN(1:end-1,1:end-1,zi)).*hFacS(1:end-1,2:end,zi),[1 1 148]).*squeeze(V(1:end-1,2:end,zi,:)).*(ssh148(1:end-1,1:end-1,:)+5);
% vfluxvD(1:end-1,1:end-1,zi,:)=W(1:end-1,1:end-1,zi+1,:).*...
%           repmat(openZd(1:end-1,1:end-1,zi).*...
%         (areasRect),[1 1 1 148]);
% vfluxvU(1:end-1,1:end-1,zi,:)=W(1:end-1,1:end-1,zi,:).*...
%           repmat(openZu(1:end-1,1:end-1,zi).*...
%         (areasRect),[1 1 1 148]);
% disp('done calc vfluxh* vfluxv*')
% max(max(max(abs(vfluxvD(:,:,1,:)))))
% max(max(max(abs(vfluxvU(:,:,1,:)))))
% max(max(max(abs(vfluxhE(:,:,1,:)))))
% max(max(max(abs(vfluxhW(:,:,1,:)))))
% max(max(max(abs(vfluxhN(:,:,1,:)))))
% max(max(max(abs(vfluxhS(:,:,1,:)))))
% clear U V W open* isopen
% 
% % % vfluxv1=zeros([700 200 15]);
% % % for i=1:15
% % %     vfluxv1(1:end-1,1:end-1,i)=openZ(1:end-1,1:end-1,i).*W(1:end-1,1:end-1,i).*...
% % %         (yvM(1:end-1,2:end)-yvM(1:end-1,1:end-1)).*...
% % %         (xuM(2:end,1:end-1)-xuM(1:end-1,1:end-1));
% % % end
% % % vfluxvA=sum(sum(sum(vfluxv1)));
% 
% vfluxh=squeeze(sum(sum(sum(vfluxhW-vfluxhE-vfluxhN+vfluxhS))));
% vfluxv=squeeze(sum(sum(sum(vfluxvD-vfluxvU))));
% vfluxvExchange=squeeze(sum(sum(sum(abs(vfluxvU+vfluxvD)))));
% vfluxhExchange=squeeze(sum(sum(sum(abs(vfluxhW)+abs(vfluxhE)+abs(vfluxhN)+abs(vfluxhS)))));
% disp('done calc exchanges')
% 
% daysPrecip=(double(hoursSince1900(608:end))-double(hoursSince1900(608)))./24;
% meanPrecipWag=-EminusPms(7,10,608:end)./3 -EminusPms(7,9,608:end)./3 -EminusPms(8,9,608:end)./3;
% areaPrecip=inWag(1:end-1,1:end-1,1).*areasRect;%(1:end-1,1:end-1);
%    %(yvM(1:end-1,2:end)-yvM(1:end-1,1:end-1)).*(xuM(2:end,1:end-1)-xuM(1:end-1,1:end-1));
% vfluxs=squeeze(sum(sum(areaPrecip))*meanPrecipWag(1:295));
% 
% volumeSSH=squeeze(sum(sum(repmat(inWag(1:end-1,1:end-1,1).*areasRect,[1 1 148]).*ssh148(1:end-1,1:end-1,:))));
% 
% figure; plot(0:147,vfluxh./(1e6))
% hold all
% plot(0:147,vfluxv./1e6)
% plot(daysPrecip(1:295),vfluxs./1e6)
% hold on; plot(0:147,(vfluxh./1e6)+(vfluxv./1e6)+(vfluxs(1:2:end)./1e6),'k','LineWidth',2)
% hold all; plot(0:147,(volumeSSH-mean(volumeSSH))./(1e6*86400))
% title('Wag Eulerian Boundary Volume Flux in Sverdrups')
% legend('horizontal','vertical','precipitation','total','WAG volume change from SSH')
% xlabel('Time in Days')
% %save2pdf('eulerWAGvolumeflux.pdf')
% 
% figure; plot(0:147,vfluxhExchange./1e6,0:147,vfluxvExchange./1e6)
% title('Wag Eulerian Boundary Exchange Volume Flux in Sverdrups')
% legend('horizontal','vertical')
% xlabel('Time in Days')
% %save2pdf('eulerWAGvolumeExchange.pdf')
% 
% 
% save('volumesfluxesWAGeuler2017.mat','-v7.3')
% disp('done volume fluxes')
% %%
% %load('varying148ts16levelsRho.mat','Rho')
% disp('calculating mass in WAG')
% 
% wagR1=zeros([700 200 15 148]);
% dBin=0.5*dInterface(2:end)+0.5*dInterface(1:end-1);
% depths1(1,1,1:15,1)=dBin(1:15);
% wagR1(1:end-1,1:end-1,2:15,:)=repmat(inWag(1:end-1,1:end-1,2:15),[1 1 1 148]).*repmat(areasRect,[1 1 14 148]).*repmat(depths1(:,:,2:15,:),[699 199 1 148]).*Rho(1:end-1,1:end-1,2:15,:);
% maxWRA=max(max(max(max(wagR1))))
% wagR1(1:end-1,1:end-1,1,:)=repmat(squeeze(inWag(1:end-1,1:end-1,1)),[1 1 148]).*repmat(areasRect,[1 1 148]).*(5+ssh148(1:end-1,1:end-1,:)).*squeeze(Rho(1:end-1,1:end-1,1,:));
% maxWR1=max(max(max(wagR1(:,:,1,:))))
% wagMass=squeeze(sum(sum(sum(wagR1))));
% %%
% %load('varying148ts16levelsRho.mat','Rho')
% %load('volumesfluxesWAGeuler2017.mat')
% rfluxhE=zeros([700 200 15 148]);
% rfluxhW=zeros([700 200 15 148]);
% rfluxhN=zeros([700 200 15 148]);
% rfluxhS=zeros([700 200 15 148]);
% 
% rfluxv1=zeros([700 200 15 148]);
% 
% %% REDO TO HAVE NSEW FLUXES for mass!
% 
% rfluxhE(2:end-1,2:end-1,:,:)=0.5.*Rho(3:end,2:end-1,1:15,:).*(vfluxhE(2:end-1,2:end-1,1:15,:)-abs(vfluxhE(2:end-1,2:end-1,1:15,:)))...
%             +0.5.*Rho(2:end-1,2:end-1,1:15,:).*(vfluxhE(2:end-1,2:end-1,1:15,:)+abs(vfluxhE(2:end-1,2:end-1,1:15,:)));
% rEmax=max(max(max(max(abs(rfluxhE)))))
% disp('rE done')
% rfluxhW(2:end-1,2:end-1,:,:)=0.5.*Rho(2:end-1,2:end-1,1:15,:).*(vfluxhW(2:end-1,2:end-1,1:15,:)-abs(vfluxhW(2:end-1,2:end-1,1:15,:)))...
%     +0.5.*Rho(1:end-2,2:end-1,1:15,:).*(vfluxhW(2:end-1,2:end-1,1:15,:)+abs(vfluxhW(2:end-1,2:end-1,1:15,:)));
% disp('rW done')
% rWmax=max(max(max(max(abs(rfluxhW)))))
% rfluxhS(2:end-1,2:end-1,:,:)=0.5.*Rho(2:end-1,2:end-1,1:15,:).*(vfluxhS(2:end-1,2:end-1,1:15,:)-abs(vfluxhS(2:end-1,2:end-1,1:15,:)))...
%    +0.5.*Rho(2:end-1,1:end-2,1:15,:).*(vfluxhS(2:end-1,2:end-1,1:15,:)+abs(vfluxhS(2:end-1,2:end-1,1:15,:)));
% disp('rS done')
% rSmax=max(max(max(max(abs(rfluxhS)))))
% rfluxhN(2:end-1,2:end-1,:,:)=0.5.*Rho(2:end-1,2:end-1,1:15,:).*(vfluxhN(2:end-1,2:end-1,1:15,:)+abs(vfluxhN(2:end-1,2:end-1,1:15,:)))...
%     +0.5.*Rho(2:end-1,3:end,1:15,:).*(vfluxhN(2:end-1,2:end-1,1:15,:)-abs(vfluxhN(2:end-1,2:end-1,1:15,:)));
% disp('rN done')
% rNmax=max(max(max(max(abs(rfluxhN)))))
% sizevU=size(vfluxvU)
% sizevD=size(vfluxvD)
% sizeRho=size(Rho)
% rfluxv1(:,:,2:15,:)=Rho(:,:,2:15,:).*0.5.*(vfluxvD(:,:,2:15,:)-abs(vfluxvD(:,:,2:15,:))+...
% 		vfluxvU(:,:,2:15,:)+abs(vfluxvU(:,:,2:15,:)))...
%             +Rho(:,:,3:16,:).*0.5.*(vfluxvD(:,:,2:15,:)+abs(vfluxvD(:,:,2:15,:)))+...
% 		Rho(:,:,1:14,:).*0.5.*(vfluxvU(:,:,2:15,:)-abs(vfluxvU(:,:,2:15,:)));
% rv2max=max(max(max(max(rfluxv1))))
% rfluxv1(:,:,1,:)=Rho(:,:,1,:).*0.5.*(vfluxvD(:,:,1,:)-abs(vfluxvD(:,:,1,:)))...
% 		+Rho(:,:,2,:).*0.5.*(vfluxvD(:,:,1,:)+abs(vfluxvD(:,:,1,:)));
% rv1max=max(max(max(rfluxv1(:,:,1,:))))
% rfluxv=squeeze(sum(sum(sum(rfluxv1))));
% rfluxh=squeeze(sum(sum(sum(rfluxhW+rfluxhS-rfluxhE-rfluxhN))));
% rfluxs=vfluxs.*1000;
% disp('r done')
% rfluxvExchange=squeeze(sum(sum(sum(abs(rfluxv1)))));
% rfluxhExchange=squeeze(sum(sum(sum(abs(rfluxhW)+abs(rfluxhE)+abs(rfluxhN)+abs(rfluxhS)))));
% disp('r exchange done')
% clear vflux*
% 
% figure; plot(0:147,rfluxh./(1e9),'LineWidth',2)
% hold all
% plot(0:147,rfluxv./1e9,'LineWidth',2)
% plot(daysPrecip(1:295),rfluxs./1e9,'LineWidth',2)
% hold on; plot(0:147,(rfluxh./1e9)+(rfluxv./1e9)+(rfluxs(1:2:end)./1e9),'k','LineWidth',4)
% hold on; plot(0:146,diff(wagMass)./(86400.*1e9),'m','LineWidth',2)
% title('WAG Eulerian Boundary Mass Flux in Sverdrups','fontsize',18)
% legend('horizontal','vertical','precipitation','total','changes in mass within WAG')
% xlabel('Time in Days','fontsize',14)
% ylabel('flux in Sv, 10^9 kg/s, mass in kg','fontsize',14)
% set(gca,'fontsize',14)
% %save2pdf('eulerWAGmassflux148.pdf')
% 
% figure; plot(0:147,rfluxhExchange./1e9,0:147,rfluxvExchange./1e9,'linewidth',2)
% title('WAG Eulerian Boundary Exchange Mass Flux in Sverdrups','fontsize',18)
% legend('horizontal','vertical')
% xlabel('Time in Days','fontsize',14)
% ylabel('flux in Sv, 10^9 kg/s','fontsize',14)
% set(gca,'fontsize',14)
% %save2pdf('eulerWAGmassExchange148.pdf')
% disp('saving')
% save('massfluxesWAGeuler.mat','-v7.3')
% disp('done mass fluxes')
% %%
% %S fluxes; h is horizontal, v is vertical, s is surface (atmosphere)
% sfluxh1=zeros([700 200 15 148]);
% sfluxv1=zeros([700 200 15 148]);
% 
% %sfluxs=zeros([700 200 15 148]);
% 
% 
% load('varying148ts16levelsRho.mat','S');
% 
% %REDO USING MASS FLUXES!!!
% 
% sfluxh1(2:end-1,2:end-1,:,:)=0.5.*S(2:end-1,2:end-1,1:15,:).*(rfluxhE(2:end-1,2:end-1,:,:)-abs(rfluxhE(2:end-1,2:end-1,:,:))...
%     +rfluxhS(2:end-1,2:end-1,:,:)-abs(rfluxhS(2:end-1,2:end-1,:,:))-rfluxhW(2:end-1,2:end-1,:,:)-abs(rfluxhW(2:end-1,2:end-1,:,:))...
%     -rfluxhN(2:end-1,2:end-1,:,:)-abs(rfluxhN(2:end-1,2:end-1,:,:)))...
%     +0.5.*S(1:end-2,2:end-1,1:15,:).*(rfluxhE(2:end-1,2:end-1,:,:)+abs(rfluxhE(2:end-1,2:end-1,:,:)))...
%     +0.5.*S(2:end-1,1:end-2,1:15,:).*(rfluxhS(2:end-1,2:end-1,:,:)+abs(rfluxhS(2:end-1,2:end-1,:,:)))...
%     -0.5.*S(3:end,2:end-1,1:15,:).*(rfluxhW(2:end-1,2:end-1,:,:)-abs(rfluxhW(2:end-1,2:end-1,:,:)))...
%     -0.5.*S(2:end-1,3:end,1:15,:).*(rfluxhN(2:end-1,2:end-1,:,:)-abs(rfluxhN(2:end-1,2:end-1,:,:)));
% 
% sfluxv1(:,:,2:15,:)=S(:,:,1:14,:).*0.5.*(rfluxv1(:,:,2:15,:)+abs(rfluxv1(:,:,2:15,:)))...
%             +S(:,:,2:15,:).*0.5.*(rfluxv1(:,:,2:15,:)-abs(rfluxv1(:,:,2:15,:)));
% sfluxv1(:,:,1,:)=S(:,:,1,:).*rfluxv1(:,:,1,:);
% clear rflux*
% sfluxh=squeeze(sum(sum(sum(sfluxh1))));
% sfluxv=squeeze(sum(sum(sum(sfluxv1))));
% save('saltfluxesWAGeuler2017.mat','sflux*','-v7.3')
% disp('done salt fluxes')
% 
% figure; plot(0:147,sfluxh./1000)%(3e10))
% hold all
% plot(0:147,sfluxv./1000)%3e10)
% legend('horizontal','vertical')
% title('Advective Salt fluxes','fontsize',18)
% xlabel('Time in days','fontsize',14)
% ylabel('Salt fluxes, kg/s of salt','fontsize',14)
% set(gca,'fontsize',14)
% 
% clear sflux*1
% %% diffusion of salt
% % addpath('/nobackup1/gbrett/mStuff')
% % load('edgesWAGeuler2017.mat','dBin','open*','inWag')
% % load('varying148ts16levelsRho.mat','S');
% % load('geometrySpinupSteady.mat') %for running on server
% 
% saltDiffHFlux=zeros([148,1]);
% saltDiffVFlux=zeros([148,1]);
% load('varying148ts16levelsRho.mat','Rho');
% load('edgesWAGeuler2017.mat');
% load('distancesAreas.mat','DXG')
% load('distancesAreas.mat','DYG')
% load('distancesAreas.mat','RAC')
% load('distancesAreas.mat','hFacW')
% load('distancesAreas.mat','hFacS')
% kh=1e4;
% kz=0.02;
% nh=1;
% 
% dBin3(1,1,1:16)=dBin(1:16);
% dBin3=repmat(dBin3,[700 200 1]);
% DXG3=repmat(DXG,[1 1 16]);
% DYG3=repmat(DYG,[1 1 16]);
% RAC3=repmat(RAC,[1 1 16]);
% for k=1:148
%     k
%     csal=S(:,:,:,k).*Rho(:,:,:,k);
%     csal(csal<1)=NaN;
%     %ch=tS(:,:,k).*cpS(:,:,k).*rhoS(:,:,k);
%     %dZ=layerdepth(:,:,k);
%     %inWAG=inWAG2(:,:,k);
%     
%     saltDiffHFlux(k)=nansum(nansum(nansum(openW(1+nh:end,:,1:16).*(-kh*(csal(1+nh:end,:,:)-csal(1:end-nh,:,:)).*dBin3(1+nh:end,:,:).*DYG3(1+nh:end,:,:)./(DXG3(1:end-nh,:,:)+DXG3(1+nh:end,:,:))))))...
%         +nansum(nansum(nansum(openE(1:end-nh,:,1:16).*(kh*(csal(1+nh:end,:,:)-csal(1:end-nh,:,:)).*dBin3(1:end-nh,:,:).*DYG3(1:end-nh,:,:)./(DXG3(1:end-nh,:,:)+DXG3(1+nh:end,:,:))))))...
%         +nansum(nansum(nansum(openS(:,1+nh:end,1:16).*(-kh*(csal(:,1+nh:end,:)-csal(:,1:end-nh,:)).*dBin3(:,1+nh:end,:).*DXG3(:,1+nh:end,:)./(DYG3(:,1:end-nh,:)+DYG3(:,1+nh:end,:))))))...
%         +nansum(nansum(nansum(openN(:,1:end-nh,1:16).*(kh*(csal(:,1+nh:end,:)-csal(:,1:end-nh,:)).*dBin3(:,1:end-nh,:).*DXG3(:,1:end-nh,:)./(DYG3(:,1:end-nh,:)+DYG3(:,1+nh:end,:))))));
%     %!!!! redo this! need to divide by dBin!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%     saltDiffVFlux(k)=kz*nansum(nansum(nansum(openZd(:,:,1:15).*(csal(:,:,2:16)-csal(:,:,1:15)).*RAC3(:,:,1:15)./diff(dBin3,1,3)))); %g salt/s    
% end
% 
% save('saltWAGeulerDiffusion2017corrected.mat','saltD*');
% 
% 
% figure; plot(0:147,sfluxh./1000,'linewidth',2)%(3e10))
% hold all
% plot(0:147,sfluxv./1000,'linewidth',2)%3e10)
% plot(0:147,saltDiffHFlux./1000,'linewidth',2);
% plot(0:147,saltDiffVFlux./1000,'linewidth',2);
% legend('horizontal advection','vertical advection','horizontal diffusion','vertical diffusion')
% title('Salt fluxes','fontsize',18)
% xlabel('Time in days','fontsize',14)
% ylabel('Salt fluxes, kg/s of salt','fontsize',14)
% set(gca,'fontsize',14)

%% heat advection, diffusion, surface forcing, plot budget

load('varying148ts16levelsRho.mat','T','Rho');
load('cp16levels148.mat', 'cp')
%load('massfluxesWAGeuler.mat','rflux*')

% %% advection
% tfluxh1=0.5.*T(2:end-1,2:end-1,1:15,:).*cp(2:end-1,2:end-1,1:15,:).*(rfluxhE(2:end-1,2:end-1,:,:)-abs(rfluxhE(2:end-1,2:end-1,:,:))...
%     +rfluxhS(2:end-1,2:end-1,:,:)-abs(rfluxhS(2:end-1,2:end-1,:,:))-rfluxhW(2:end-1,2:end-1,:,:)-abs(rfluxhW(2:end-1,2:end-1,:,:))...
%     -rfluxhN(2:end-1,2:end-1,:,:)-abs(rfluxhN(2:end-1,2:end-1,:,:)))...
%     +0.5.*T(1:end-2,2:end-1,1:15,:).*cp(1:end-2,2:end-1,1:15,:).*(rfluxhE(2:end-1,2:end-1,:,:)+abs(rfluxhE(2:end-1,2:end-1,:,:)))...
%     +0.5.*T(2:end-1,1:end-2,1:15,:).*cp(2:end-1,1:end-2,1:15,:).*(rfluxhS(2:end-1,2:end-1,:,:)+abs(rfluxhS(2:end-1,2:end-1,:,:)))...
%     -0.5.*T(3:end,2:end-1,1:15,:).*cp(3:end,2:end-1,1:15,:).*(rfluxhW(2:end-1,2:end-1,:,:)-abs(rfluxhW(2:end-1,2:end-1,:,:)))...
%     -0.5.*T(2:end-1,3:end,1:15,:).*cp(2:end-1,3:end,1:15,:).*(rfluxhN(2:end-1,2:end-1,:,:)-abs(rfluxhN(2:end-1,2:end-1,:,:)));
% disp('line 1 done')
% tfluxv1(:,:,2:15,:)=T(:,:,1:14,:).*cp(:,:,1:14,:).*0.5.*(rfluxv1(:,:,2:15,:)+abs(rfluxv1(:,:,2:15,:)))...
%             +T(:,:,2:15,:).*cp(:,:,2:15,:).*0.5.*(rfluxv1(:,:,2:15,:)-abs(rfluxv1(:,:,2:15,:)));
% tfluxv1(:,:,1,:)=T(:,:,1,:).*cp(:,:,1,:).*rfluxv1(:,:,1,:);
% clear rflux*
%  %advection totals
% tfluxh=squeeze(nansum(nansum(nansum(tfluxh1))));
% tfluxv=squeeze(nansum(nansum(nansum(tfluxv1))));
% disp('heat advection done')
% clear tfluxv1 tfluxh1
%% diffusion
% addpath('/nobackup1/gbrett/mStuff')
% load('edgesWAGeuler2017.mat','dBin','open*','inWag')
% load('distancesAreas.mat','DXG')
% load('distancesAreas.mat','DYG')
% load('distancesAreas.mat','RAC')
% load('distancesAreas.mat','hFacW')
% load('distancesAreas.mat','hFacS')
% dBin3(1,1,1:16)=dBin(1:16);
% dBin3=repmat(dBin3,[700 200 1]);
% DXG3=repmat(DXG,[1 1 16]);
% DYG3=repmat(DYG,[1 1 16]);
% RAC3=repmat(RAC,[1 1 16]);
% nh=1
% kh=1e4;
% kz=0.02;
% load('geometrySpinupSteady.mat')
% for k=1:148
%     k
%     csal=T(:,:,:,k).*cp(:,:,:,k).*Rho(:,:,:,k);
%     csal(csal<1)=NaN;
%     %ch=tS(:,:,k).*cpS(:,:,k).*rhoS(:,:,k);
%     %dZ=layerdepth(:,:,k);
%     %inWAG=inWAG2(:,:,k);
%     
%     heatDiffHFlux(k)=nansum(nansum(nansum(openW(1+nh:end,:,1:16).*(-kh*(csal(1+nh:end,:,:)-csal(1:end-nh,:,:)).*dBin3(1+nh:end,:,:).*DYG3(1+nh:end,:,:)./(DXG3(1:end-nh,:,:)+DXG3(1+nh:end,:,:))))))...
%         +nansum(nansum(nansum(openE(1:end-nh,:,1:16).*(kh*(csal(1+nh:end,:,:)-csal(1:end-nh,:,:)).*dBin3(1:end-nh,:,:).*DYG3(1:end-nh,:,:)./(DXG3(1:end-nh,:,:)+DXG3(1+nh:end,:,:))))))...
%         +nansum(nansum(nansum(openS(:,1+nh:end,1:16).*(-kh*(csal(:,1+nh:end,:)-csal(:,1:end-nh,:)).*dBin3(:,1+nh:end,:).*DXG3(:,1+nh:end,:)./(DYG3(:,1:end-nh,:)+DYG3(:,1+nh:end,:))))))...
%         +nansum(nansum(nansum(openN(:,1:end-nh,1:16).*(kh*(csal(:,1+nh:end,:)-csal(:,1:end-nh,:)).*dBin3(:,1:end-nh,:).*DXG3(:,1:end-nh,:)./(DYG3(:,1:end-nh,:)+DYG3(:,1+nh:end,:))))));
%     heatDiffVFlux(k)=kz*nansum(nansum(nansum(openZd(:,:,1:15).*(csal(:,:,2:16)-csal(:,:,1:15)).*RAC3(:,:,1:15)./diff(dBin3,1,3))));    
% end
% clear open* T cp Rho
% disp('heat diffusion done')
% 
% %%
% %surface heat flux
% %load('C:\Users\JayB\Desktop\climatologyIfremar\eraInterim200720082009.mat')
% load('eraInterim200720082009.mat')
% load('eraDataAndFall2007meansTry2.mat')
% iNov12007=608;
% i148=904;%(march 28 2008)
% h148=hflux2(:,:,iNov12007:i148);
% inWag2=inWag(:,:,1);
% for i=1:148
%    hNow=griddata(double(lonEg).',double(latEg).',mean(h148(:,:,(2*i-1):2*i),3),XC,YC);
%    eulerWAGheatInt(i)=-nansum(nansum(hNow.*RAC.*inWag2));
% end
% clear inWag lon* lat* h148 hflux2
% disp('surface heat done, saving')
% 
% load('heatAdvectionWAGeuler2017.mat')
% save('heatWAGeuler2017.mat','tflux*','heat*','euler*','-v7.3')
% 
% figure; plot(0:147,tfluxh,'linewidth',2)%(3e10))
% hold all
% plot(0:147,tfluxv,'linewidth',2)%3e10)
% plot(0:147,heatDiffHFlux,'linewidth',2);
% plot(0:147,heatDiffVFlux,'linewidth',2);
% plot(0:147,eulerWAGheatInt,'linewidth',2)
% plot(0:147,tfluxh(:)+tfluxv(:)+heatDiffHFlux(:)+heatDiffVFlux(:)+eulerWAGheatInt(:),'k--','linewidth',4)
% legend('horizontal advection','vertical advection','horizontal diffusion','vertical diffusion','surface forcing','total')
% title('Heat fluxes','fontsize',18)
% xlabel('Time in days','fontsize',14)
% ylabel('Heat fluxes, Watts','fontsize',14)
% set(gca,'fontsize',14)
% save2pdf('heatBudgetEuler2017.pdf')


%% quantities in eulerian wag over time

addpath('/nobackup1/gbrett/mStuff')
load('geometrySpinupSteady.mat')
load('edgesWAGeuler2017.mat','dBin','inWag')
load('distancesAreas.mat','RAC')
load('ssh148.mat','ssh148');
RAC3=repmat(RAC,[1 1 16]);
dBin3(1,1,1:16)=dBin(1:16);
dBin3=repmat(dBin3,[700 200 1]);
inWag=inWag(:,:,1:16);

volWagEuler1=sum(sum(sum(RAC3.*inWag.*dBin3)));
volWagEuler148=volWagEuler1*ones(148,1)+squeeze(sum(sum(ssh148.*repmat(inWag(:,:,1).*RAC,[1 1 148]))));

load('varying148ts16levelsRho.mat','Rho');

massWagEuler=squeeze(sum(sum(sum(repmat(RAC3.*inWag.*dBin3,[1 1 1 148]).*Rho(:,:,1:16,:)))))+squeeze(sum(sum(ssh148.*squeeze(Rho(:,:,1,:)).*repmat(inWag(:,:,1).*RAC,[1 1 148]))));

load('varying148ts16levelsRho.mat','S');
saltWagEuler=squeeze(sum(sum(sum(repmat(RAC3.*inWag.*dBin3,[1 1 1 148]).*Rho(:,:,1:16,:).*S(:,:,1:16,:)))))+squeeze(sum(sum(ssh148.*squeeze(Rho(:,:,1,:).*S(:,:,1,:)).*repmat(inWag(:,:,1).*RAC,[1 1 148]))));
clear S
load('cp16levels148.mat', 'cp')
load('varying148ts16levelsRho.mat','T');
heatWagEuler=squeeze(sum(sum(sum(repmat(RAC3.*inWag.*dBin3,[1 1 1 148]).*Rho(:,:,1:16,:).*cp(:,:,1:16,:).*T(:,:,1:16,:)))))+squeeze(sum(sum(ssh148.*squeeze(Rho(:,:,1,:).*cp(:,:,1,:).*T(:,:,1,:)).*repmat(inWag(:,:,1).*RAC,[1 1 148]))));

clear T Rho cp *3 

save('contentsWAGeuler2017.mat')
