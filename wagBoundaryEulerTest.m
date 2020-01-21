%Eulerian WAG boundaries for first 148 days of the simulation
%Nov through most of Mar
%Boundaries by limits in time-averaged T,S
%Fluxes from time-varying U,V
%addpath('../mStuff')
% load('steady148tswRho.mat');
% 
% wagS=zeros([700 200 16]);
% for i=1:10
%     wagS(:,:,i)=Sal(:,:,i)<36.5;
% end
% for i=11:15
%     wagS(:,:,i)=Sal(:,:,i)<(36.5+0.1*(i-10));
% end
% 
% wagS(1:249,:,:)=0;
% wagS(427:end,:,:)=0;
% 
% wagT=zeros(size(wagS));
% for i=1:5
%    wagT(:,:,i)=Temp(:,:,i)>(18.1-0.1*i); 
% end
% for i=6:10
%     wagT(:,:,i)=Temp(:,:,i)>(18.1-0.15*i);
% end
% for i=11:15
%     wagT(:,:,i)=Temp(:,:,i)>(18.1-0.2*i);
% end
% 
% wagT(1:249,:,:)=0;
% wagT(427:end,:,:)=0;
% 
% wagR=zeros(size(wagS));
% for i=1:16
%     wagR(:,:,i)=Sigma(:,:,i)<26.6;
% end
% wagR(1:249,:,:)=0;
% wagR(427:end,:,:)=0;
% 
% 
% clear Sal Temp W
% %%
% inWag=wagS.*wagT; clear wagS wagT
% openW=zeros(size(inWag));
% openE=zeros(size(inWag));
% openS=zeros(size(inWag));
% openN=zeros(size(inWag));
% isedge=zeros(size(inWag));
% 
% for zi=1:14
%     for xi=249:428
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
% 
% figure
% plot(lonCoast,latCoast,'k')
% hold all
% for i=1:14
% plot(XC(logical(isedge(:,:,i))),YC(logical(isedge(:,:,i))),'.')
% end
% %dBin(1:14)
% legend('Coast','2.5m','8m','14m','20.5m','27.5m','35m','43.5m','53m','63.5m','75m','87.5m','101.5m','117m','134m')
% title('WAG Edge with depth, Defined by T-S contours')
% axis ([-6 -1.5 34.75 37.25])
% %%
% %vertical differences to get vertical fluxes
% openZ=zeros(size(inWag));
% for i=1:15
%    openZ(:,:,i)=(inWag(:,:,i)==1 & inWag(:,:,i+1)==0); 
% end
% save('edgesWAGeuler.mat','-v7.3')
% disp('done section 1')
%%
addpath('../mStuff')
load('edgesWAGeuler.mat')
load('uva148levels14.mat','U');
load('uva148levels14.mat','V');
load('wa148levels16.mat','W');
load('eraDataAndFall2007meansTry2.mat', 'EminusPms')
load('eraDataAndFall2007meansTry2.mat', 'hoursSince1900')
load('eraDataAndFall2007meansTry2.mat', 'hflux2')
load('ssh148.mat', 'ssh148')
load('distancesAreas.mat','DXG')
load('distancesAreas.mat','DYG')
load('distancesAreas.mat','RAC')
%xcM=(XC-min(min(XC))*ones(size(XC))).*111000.*cosd(YC); 
%ycM=(YC-min(min(YC))*ones(size(YC))).*111000;
%xuM=(XU-min(min(XU))*ones(size(XU))).*111000.*cosd(YU); 
%yvM=(YV-min(min(YV))*ones(size(YV))).*111000;
distx=DXG(1:end-1,:);%distance(XG(1:end-1,:),YG(1:end-1,:),XG(2:end,:),YG(2:end,:),[6371000 0]);
disty=DYG(:,1:end-1);%distance(XG(:,1:end-1),YG(:,1:end-1),XG(:,2:end),YG(:,2:end),[6371000 0]);
areasRectApprox=RAC(1:end-1,1:end-1);%(0.5.*distx(:,1:end-1)+0.5.*distx(:,2:end)).*(0.5.*disty(1:end-1,:)+0.5.*disty(2:end,:));

openW=zeros([700 200 16]);
openW(350,75,5)=1;
openE=openW;
openN=openW;
openS=openW;
openZd=openW;
openZu=openW;
inWag=openW;

% openZd=openZ;
% openZu=zeros(size(inWag));
% for i=2:16
%    openZu(:,:,i)=(inWag(:,:,i)==1 & inWag(:,:,i-1)==0); 
% end

%horizontal volume flux
vfluxhE=zeros([700 200 15 148]);
vfluxhW=zeros([700 200 15 148]);
vfluxhN=zeros([700 200 15 148]);
vfluxhS=zeros([700 200 15 148]);
vfluxvD=zeros([700 200 15 148]);
vfluxvU=zeros([700 200 15 148]);

disp('start zi=2:15')
for zi=2:15
      vfluxhW(1:end-1,1:end-1,zi,:)=(dInterface(zi+1)-dInterface(zi)).*repmat((disty(1:end-1,:))...
          .*openW(1:end-1,1:end-1,zi),[1 1 1 148]).*U(1:end-1,1:end-1,zi,:);
      vfluxhE(1:end-1,1:end-1,zi,:)=(dInterface(zi+1)-dInterface(zi)).*repmat((disty(2:end,:))...
          .*openE(1:end-1,1:end-1,zi),[1 1 1 148]).*U(2:end,1:end-1,zi,:);
      vfluxhS(1:end-1,1:end-1,zi,:)=(dInterface(zi+1)-dInterface(zi)).*repmat((distx(:,1:end-1))...
          .*openS(1:end-1,1:end-1,zi),[1 1 1 148]).*V(1:end-1,1:end-1,zi,:);
      vfluxhN(1:end-1,1:end-1,zi,:)=(dInterface(zi+1)-dInterface(zi)).*repmat((distx(:,2:end))...
          .*openN(1:end-1,1:end-1,zi),[1 1 1 148]).*V(1:end-1,2:end,zi,:);
      vfluxvD(1:end-1,1:end-1,zi,:)=W(1:end-1,1:end-1,zi+1,:).*...
          repmat(openZd(1:end-1,1:end-1,zi).*...
        areasRectApprox,[1 1 1 148]);
      vfluxvU(1:end-1,1:end-1,zi,:)=W(1:end-1,1:end-1,zi,:).*...
          repmat(openZu(1:end-1,1:end-1,zi).*...
        areasRectApprox,[1 1 1 148]);
end
disp('done zi=2:15')
max(max(max(max(abs(vfluxhW)))))
max(max(max(max(abs(vfluxhE)))))
max(max(max(max(abs(vfluxhN)))))
max(max(max(max(abs(vfluxhS)))))
max(max(max(max(abs(vfluxvD)))))
max(max(max(max(abs(vfluxvU)))))
zi=1
vfluxhW(1:end-1,1:end-1,zi,:)=repmat((disty(1:end-1,:))...
          .*squeeze(openW(1:end-1,1:end-1,zi)),[1 1 148]).*squeeze(U(1:end-1,1:end-1,zi,:)).*(ssh148(1:end-1,1:end-1,:)+5);
vfluxhE(1:end-1,1:end-1,zi,:)=repmat((disty(2:end,:))...
          .*squeeze(openE(1:end-1,1:end-1,zi)),[1 1 148]).*squeeze(U(2:end,1:end-1,zi,:)).*(ssh148(1:end-1,1:end-1,:)+5);
vfluxhS(1:end-1,1:end-1,zi,:)=repmat((distx(:,1:end-1))...
          .*squeeze(openS(1:end-1,1:end-1,zi)),[1 1 148]).*squeeze(V(1:end-1,1:end-1,zi,:)).*(ssh148(1:end-1,1:end-1,:)+5);
vfluxhN(1:end-1,1:end-1,zi,:)=repmat((distx(:,2:end))...
          .*squeeze(openN(1:end-1,1:end-1,zi)),[1 1 148]).*squeeze(V(1:end-1,2:end,zi,:)).*(ssh148(1:end-1,1:end-1,:)+5);
vfluxvD(1:end-1,1:end-1,zi,:)=W(1:end-1,1:end-1,zi+1,:).*...
          repmat(openZd(1:end-1,1:end-1,zi).*...
        (disty(1:end-1,:)).*(distx(:,1:end-1)),[1 1 1 148]);
vfluxvU(1:end-1,1:end-1,zi,:)=W(1:end-1,1:end-1,zi,:).*...
          repmat(openZu(1:end-1,1:end-1,zi).*...
        (disty(1:end-1,:)).*(distx(:,1:end-1)),[1 1 1 148]);
disp('done calc vfluxh* vfluxv1')
max(max(max(abs(vfluxvD(:,:,1,:)))))
max(max(max(abs(vfluxvU(:,:,1,:)))))
max(max(max(abs(vfluxhE(:,:,1,:)))))
max(max(max(abs(vfluxhW(:,:,1,:)))))
max(max(max(abs(vfluxhN(:,:,1,:)))))
max(max(max(abs(vfluxhS(:,:,1,:)))))
clear U V W open* isopen

% % vfluxv1=zeros([700 200 15]);
% % for i=1:15
% %     vfluxv1(1:end-1,1:end-1,i)=openZ(1:end-1,1:end-1,i).*W(1:end-1,1:end-1,i).*...
% %         (yvM(1:end-1,2:end)-yvM(1:end-1,1:end-1)).*...
% %         (xuM(2:end,1:end-1)-xuM(1:end-1,1:end-1));
% % end
% % vfluxvA=sum(sum(sum(vfluxv1)));

vfluxh=squeeze(sum(sum(sum(vfluxhW-vfluxhE-vfluxhN+vfluxhS))));
vfluxv=squeeze(sum(sum(sum(vfluxvD-vfluxvU))));
vfluxvExchange=squeeze(sum(sum(sum(abs(vfluxvU+vfluxvD)))));
vfluxhExchange=squeeze(sum(sum(sum(abs(vfluxhW)+abs(vfluxhE)+abs(vfluxhN)+abs(vfluxhS)))));
disp('done calc exchanges')

daysPrecip=(double(hoursSince1900(608:end))-double(hoursSince1900(608)))./24;
meanPrecipWag=-EminusPms(7,10,608:end)./3 -EminusPms(7,9,608:end)./3 -EminusPms(8,9,608:end)./3;
areaPrecip=0;%inWag(1:end-1,1:end-1,1).*...
   %(yvM(1:end-1,2:end)-yvM(1:end-1,1:end-1)).*(xuM(2:end,1:end-1)-xuM(1:end-1,1:end-1));
vfluxs=squeeze(sum(sum(areaPrecip))*meanPrecipWag(1:295));

%volumeSSH=squeeze(sum(sum(repmat(inWag(1:end-1,1:end-1,1).*distx(:,1:end-1).*disty(1:end-1,:),[1 1 148]).*ssh148(1:end-1,1:end-1,:))));

figure; plot(0:147,vfluxh./(1e6))
hold all
plot(0:147,vfluxv./1e6)
plot(daysPrecip(1:295),vfluxs./1e6)
hold on; plot(0:147,(vfluxh./1e6)+(vfluxv./1e6)+(vfluxs(1:2:end)./1e6),'k','LineWidth',2)
%hold all; plot(0:147,(volumeSSH-mean(volumeSSH))./(1e6*86400))
title('Wag Eulerian Boundary Volume Flux in Sverdrups')
legend('horizontal','vertical','precipitation','total','WAG volume change from SSH')
xlabel('Time in Days')
save2pdf('eulerWAGvolumefluxTestOne.pdf')

figure; plot(0:147,vfluxhExchange./1e6,0:147,vfluxvExchange./1e6)
title('Wag Eulerian Boundary Exchange Volume Flux in Sverdrups')
legend('horizontal','vertical')
xlabel('Time in Days')
save2pdf('eulerWAGvolumeExchangeTestOne.pdf')


save('volumesfluxesWAGeulerTestOne.mat','-v7.3')
disp('done volume fluxes')

%%
% load('varying148ts16levelsRho.mat','Rho')
% load('volumesfluxesWAGeuler.mat')
% rfluxhE=zeros([700 200 15 148]);
% rfluxhW=zeros([700 200 15 148]);
% rfluxhN=zeros([700 200 15 148]);
% rfluxhS=zeros([700 200 15 148]);
% 
% rfluxv1=zeros([700 200 15 148]);
% 
% %%REDO TO HAVE NSEW FLUXES for mass!
% 
% rfluxhE(2:end-1,2:end-1,:,:)=0.5.*Rho(2:end-1,2:end-1,1:15,:).*(vfluxhE(2:end-1,2:end-1,:,:)-abs(vfluxhE(2:end-1,2:end-1,:,:)))...
%             +0.5.*Rho(1:end-2,2:end-1,1:15,:).*(vfluxhE(2:end-1,2:end-1,:,:)+abs(vfluxhE(2:end-1,2:end-1,:,:)));
% rEmax=max(max(max(max(abs(rfluxhE)))))
% disp('rE done')
% rfluxhW(2:end-1,2:end-1,:,:)=0.5.*Rho(3:end,2:end-1,1:15,:).*(vfluxhW(2:end-1,2:end-1,:,:)-abs(vfluxhW(2:end-1,2:end-1,:,:)))...
%     +0.5.*Rho(2:end-1,2:end-1,1:15,:).*(vfluxhW(2:end-1,2:end-1,:,:)-abs(vfluxhW(2:end-1,2:end-1,:,:)));
% disp('rW done')
% rWmax=max(max(max(max(abs(rfluxhW)))))
% rfluxhS(2:end-1,2:end-1,:,:)=0.5.*Rho(2:end-1,2:end-1,1:15,:).*(vfluxhS(2:end-1,2:end-1,:,:)-abs(vfluxhS(2:end-1,2:end-1,:,:)))...
%    +0.5.*Rho(2:end-1,1:end-2,1:15,:).*(vfluxhS(2:end-1,2:end-1,:,:)+abs(vfluxhS(2:end-1,2:end-1,:,:)));
% disp('rS done')
% rSmax=max(max(max(max(abs(rfluxhS)))))
% rfluxhN(2:end-1,2:end-1,:,:)=0.5.*Rho(2:end-1,2:end-1,1:15,:).*(vfluxhN(2:end-1,2:end-1,:,:)-abs(vfluxhN(2:end-1,2:end-1,:,:)))...
%     +0.5.*Rho(2:end-1,3:end,1:15,:).*(vfluxhN(2:end-1,2:end-1,:,:)-abs(vfluxhN(2:end-1,2:end-1,:,:)));
% disp('rN done')
% rNmax=max(max(max(max(abs(rfluxhN)))))
% rfluxv1(:,:,2:15,:)=Rho(:,:,1:14,:).*0.5.*(vfluxv1(:,:,2:15,:)+abs(vfluxv1(:,:,2:15,:)))...
%             +Rho(:,:,2:15,:).*0.5.*(vfluxv1(:,:,2:15,:)-abs(vfluxv1(:,:,2:15,:)));
% rv2max=max(max(max(max(rfluxv1))))
% rfluxv1(:,:,1,:)=Rho(:,:,1,:).*vfluxv1(:,:,1,:);
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
% figure; plot(0:147,rfluxh./(1e9))
% hold all
% plot(0:147,rfluxv./1e9)
% plot(daysPrecip(1:295),rfluxs./1e9)
% hold on; plot(0:147,(rfluxh./1e9)+(rfluxv./1e9)+(rfluxs(1:2:end)./1e9),'k','LineWidth',2)
% title('WAG Eulerian Boundary Mass Flux in Sverdrups')
% legend('horizontal','vertical','precipitation','total')
% xlabel('Time in Days')
% ylabel('flux in 10^9 kg')
% save2pdf('eulerWAGmassflux148.pdf')
% 
% figure; plot(0:147,rfluxhExchange./1e9,0:147,rfluxvExchange./1e9)
% title('WAG Eulerian Boundary Exchange Mass Flux in Sverdrups')
% legend('horizontal','vertical')
% xlabel('Time in Days')
% ylabel('flux in 10^9 kg')
% save2pdf('eulerWAGmassExchange148.pdf')
% 
% save('massfluxesWAGeuler.mat','-v7.3')
% disp('done mass fluxes')
% %%
% %T,S fluxes; h is horizontal, v is vertical, s is surface (atmosphere)
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
% sfluxv1(2:end-1,2:end-1,2:15,:)=S(:,:,1:14,:).*0.5.*(rfluxv1(:,:,2:15,:)+abs(rfluxv1(:,:,2:15,:)))...
%             +S(:,:,2:15,:).*0.5.*(rfluxv1(:,:,2:15,:)-abs(rfluxv1(:,:,2:15,:)));
% sfluxv1(:,:,1,:)=S(:,:,1,:).*rfluxv1(:,:,1,:);
% clear rflux*
% save('saltfluxesWAGeuler.mat','-v7.3')
% disp('done salt fluxes')
% sfluxh=squeeze(sum(sum(sum(sfluxh1))));
% sfluxv=squeeze(sum(sum(sum(sfluxv1))));
% figure; plot(0:147,sfluxh./(3e10))
% hold all
% plot(0:147,sfluxv./3e10)
% legend('horizontal','vertical')
% title('Salt fluxes')
% xlabel('Days')
% ylabel('Sverdrup equivalents, 3*10^7 kg/s of salt')
%%
% load('varying148ts16levelsRho.mat','T');
% load('massfluxesWAGeuler.mat');
% tfluxv1=zeros([700 200 15 148]);
% tfluxh1=zeros([700 200 15 148]);
% %C_P!!!
% tfluxh1(2:end-1,2:end-1,:,:)=0.5.*S(2:end-1,2:end-1,1:15,:).*(rfluxhE(2:end-1,2:end-1,:,:)-abs(rfluxhE(2:end-1,2:end-1,:,:))...
%     +rfluxhS(2:end-1,2:end-1,:,:)-abs(rfluxhS(2:end-1,2:end-1,:,:))-rfluxhW(2:end-1,2:end-1,:,:)-abs(rfluxhW(2:end-1,2:end-1,:,:))...
%     -rfluxhN(2:end-1,2:end-1,:,:)-abs(rfluxhN(2:end-1,2:end-1,:,:)))...
%     +0.5.*S(1:end-2,2:end-1,1:15,:).*(rfluxhE(2:end-1,2:end-1,:,:)+abs(rfluxhE(2:end-1,2:end-1,:,:)))...
%     +0.5.*S(2:end-1,1:end-2,1:15,:).*(rfluxhS(2:end-1,2:end-1,:,:)+abs(rfluxhS(2:end-1,2:end-1,:,:)))...
%     -0.5.*S(3:end,2:end-1,1:15,:).*(rfluxhW(2:end-1,2:end-1,:,:)-abs(rfluxhW(2:end-1,2:end-1,:,:)))...
%     -0.5.*S(2:end-1,3:end,1:15,:).*(rfluxhN(2:end-1,2:end-1,:,:)-abs(rfluxhN(2:end-1,2:end-1,:,:)));
% tfluxv1(2:end-1,2:end-1,zi,ti)=T(:,:,1:14,:).*0.5.*(rfluxv1(:,:,2:15,:)+abs(rfluxv1(:,:,2:15,:)))...
%             +T(:,:,2:15,:).*0.5.*(rfluxv1(:,:,2:15,:)-abs(rfluxv1(:,:,2:15,:)));
% tfluxv1(:,:,1,:)=T(:,:,1,:).*rfluxv1(:,:,1,:);
% clear rflux*
% tfluxh=squeeze(sum(sum(sum(tfluxh1))));
% tfluxv=squeeze(sum(sum(sum(tfluxv1))));
% save('heatfluxesWAGeuler.mat','-v7.3')
% disp('done heat fluxes')
% 
% figure; plot(0:147,tfluxh)
% hold all
% plot(0:147,tfluxv)
% legend('horizontal','vertical')
% title('Heat fluxes')
% xlabel('Days')
%ylabel('Sverdrup equivalents, 3*10^7 kg/s of salt')


% tfluxs=zeros([700 200 15 148]);

% 
% tfluxh1(2:end-1,2:end-1,:,:)=0.5.*T(2:end-1,2:end-1,1:15,:).*(vfluxhE(2:end-1,2:end-1,:,:)-abs(vfluxhE(2:end-1,2:end-1,:,:))...
%     +vfluxhS(2:end-1,2:end-1,:,:)-abs(vfluxhS(2:end-1,2:end-1,:,:))-vfluxhW(2:end-1,2:end-1,:,:)-abs(vfluxhW(2:end-1,2:end-1,:,:))...
%     -vfluxhN(2:end-1,2:end-1,:,:)-abs(vfluxhN(2:end-1,2:end-1,:,:)))...
%     +0.5.*T(1:end-2,2:end-1,1:15,:).*(vfluxhE(2:end-1,2:end-1,:,:)+abs(vfluxhE(2:end-1,2:end-1,:,:)))...
%     +0.5.*T(2:end-1,1:end-2,1:15,:).*(vfluxhS(2:end-1,2:end-1,:,:)+abs(vfluxhS(2:end-1,2:end-1,:,:)))...
%     -0.5.*T(3:end,2:end-1,1:15,:).*(vfluxhW(2:end-1,2:end-1,:,:)-abs(vfluxhW(2:end-1,2:end-1,:,:)))...
%     -0.5.*T(2:end-1,3:end,1:15,:).*(vfluxhN(2:end-1,2:end-1,:,:)-abs(vfluxhN(2:end-1,2:end-1,:,:)));
% 
% 
% tfluxv1(:,:,2:15,:)=T(:,:,1:14,:).*0.5.*(vfluxv1(:,:,2:15,:)+abs(vfluxv1(:,:,2:15,:)))...
%             +T(:,:,2:15,:).*0.5.*(vfluxv1(:,:,2:15,:)-abs(vfluxv1(:,:,2:15,:)));
% tfluxv1(:,:,1,:)=T(:,:,1,:).*vfluxv1(:,:,1,:);
% 
% 
% %Swag=squeeze(sum(sum(sum(S(:,:,1:15,:).*repmat(inWag,[1 1 1 148]))))./);
% %Twag=squeeze(sum(sum(sum(T(:,:,1:15,:).*repmat(inWag,[1 1 1 148]))))./sum(sum(sum(inWag))));
% areas=(xcM(2:end,2:end)-xcM(1:end-1,2:end)).*(ycM(2:end,2:end)-ycM(2:end,1:end-1));
% wagVolume=zeros(15,1);
% Swag=zeros(15,148);
% Twag=Swag;
% for i=1:15
%     wagVolume(i)=sum(sum(inWag(2:end,2:end,i).*areas.*(dInterface(i+1)-dInterface(i))));
%     Swag(i,:)=sum(sum(repmat(inWag(2:end,2:end,i).*areas.*(dInterface(i+1)-dInterface(i)),[1 1 1 148]).*S(2:end,2:end,i,:)));
%     Twag(i,:)=sum(sum(repmat(inWag(2:end,2:end,i).*areas.*(dInterface(i+1)-dInterface(i)),[1 1 1 148]).*T(2:end,2:end,i,:)));
% end
% 
% 
% sfluxv=squeeze(sum(sum(sum(sfluxv1))));
% sfluxh=squeeze(sum(sum(sum(sfluxh1))));
% tfluxv=squeeze(sum(sum(sum(tfluxv1))));
% tfluxh=squeeze(sum(sum(sum(tfluxh1))));
% 
% meanHfluxWag=-hflux2(7,10,608:end)./3 -hflux2(7,9,608:end)./3 -hflux2(8,9,608:end)./3;
% hfluxs=squeeze(sum(sum(areaPrecip))*meanHfluxWag(1:295));
% tfluxs=hfluxs/(3990*1021);  %cp 3990 J/(kg degrees), density 1021 kg/m3; hflux J/s; tflux m^3 degrees/s
% 
% 
% figure; plot(squeeze(sum(Swag))-sum(Swag(:,1))); hold all; plot(86400.*sfluxh); plot(86400.*sfluxv);
% legend('changes in daily WAG Salinity*Volume','horizontal salinity flux, m^3/day','vertical salinity flux, m^3/day')
% xlabel('Time in days')
% save2pdf('salinityfluxWAGeuler.pdf')
% figure; plot(squeeze(sum(Twag))-sum(Twag(:,1))); hold all; plot(86400*tfluxh); plot(86400*tfluxv); 
% legend('changes in daily WAG Temperature*Volume','horizontal temperature flux, m^3/day','vertical temperature flux, m^3/day')
% xlabel('Time in days')
% save2pdf('temperaturefluxWAGeuler.pdf')
% save('tsfluxWAGeuler.mat','-v7.3')
% disp('done ts fluxes, done all')

