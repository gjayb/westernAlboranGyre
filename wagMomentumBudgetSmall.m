%%load variables
 times=8640:8640:(14*8640)
 Utend=rdmds('Utend',times);
 sizeUtend=size(Utend)
 Vtend=rdmds('Vtend',times);
 UDiss=rdmds('UDiss',times);
 VDiss=rdmds('VDiss',times);
 AdvU=rdmds('AdvU',times);
 AdvV=rdmds('AdvV',times);
 UCori=rdmds('UCori',times);
 VCori=rdmds('VCori',times);
 UdPdx=rdmds('UdPdx',times);
 VdPdy=rdmds('VdPdy',times);
 Uext=rdmds('Uext',times);
 Vext=rdmds('Vext',times);
 VisZU=rdmds('ViscZiU',times);
 VisZV=rdmds('ViscZiV',times);
 
 SSH=rdmds('SSHave',times);
 AbU=rdmds('ABU',times);
 AbV=rdmds('ABV',times);
 
 %AdvXU=rdmds('AdvXU',times);%zero
 %AdvYU=rdmds('AdvYU',times);%zero
 %AdvZU=rdmds('AdvZU',times);%zero
 AdvVortU=rdmds('AdvVortU',times);
 AdvReU=rdmds('AdvReU',times);%nonzero!
 %AdvXV=rdmds('AdvXV',times);%zero
 %AdvYV=rdmds('AdvYV',times);%zero
 %AdvZV=rdmds('AdvZV',times);%zero
 AdvVortV=rdmds('AdvVortV',times);%nonzero!
 AdvReV=rdmds('AdvReV',times);%wasn't saved whoops
 
 KppK=rdmds('KppVis',times);
 Usd=rdmds('Usd',times);
 Vsd=rdmds('Vsd',times);
 Ubd=rdmds('Ubd',times);
 Vbd=rdmds('Vbd',times);
 disp('saving diagnostics')
 save('momentumDiagnostics14daySplit.mat','-v7.3')
%%
%load('momentumDiagnostics11day.mat')

sizeUtend=size(Utend)
load('edgesWAGeuler2017NF275S.mat','inWag','open*','dInterface','XC','YC','*Coast')
sizeInWag=size(inWag)
inWag1=inWag(:,:,1);
inWag(:,:,1:15)=repmat(inWag1,[1 1 15]);
%inWag=zeros(size(inWag));
%inWag(300:301,100:101,1:15)=1;
% %inWag(300,100,1:16)=1;
load('distancesAreas.mat','ra*','hFac*','dxc','dyc')
DXC=reshape(dxc,[700 200]); DYC=reshape(dyc,[700 200]);
rAw=reshape(raw,[700 200]); rAs=reshape(ras,[700 200]);
load('geometrySpinupSteady.mat','XC','YC','dInterface')
xmin=min(min(XC));
ymin=min(min(YC));
xcm=(XC-xmin).*111000.*cosd(YC);
ycm=(YC-ymin).*111000;
% % load('varying148ts16levelsRho.mat', 'Rho')
% % Rho=squeeze(Rho(:,:,1,1:58));
% % load('cp16levels148.mat', 'cp')
% % cp=squeeze(cp(:,:,1,1:58));
g=9.81;

nt=14;

dZ(1,1,1:46)=diff(dInterface);
cellVol=repmat(rAw,[1 1 46]).*hFacW.*repmat(dZ,[700 200 1]);
cellVolU=repmat(cellVol,[1 1 1 nt]);
cellVol=repmat(rAs,[1 1 46]).*hFacS.*repmat(dZ,[700 200 1]);
cellVolV=repmat(cellVol,[1 1 1 nt]);

%%
dUdt=squeeze(nansum(nansum(nansum(Utend.*repmat(inWag,[1 1 1 nt])./86400))));
dVdt=squeeze(nansum(nansum(nansum(Vtend.*repmat(inWag,[1 1 1 nt])./86400))));


%need to offset press1 by press2 before inwag and adding up
dhdx=zeros([700 200 1 nt]);
dhdy=dhdx;
% for i=1:nt
% [dhdy1,dhdx1]=gradient(SSH(:,:,i),ycm,xcm);
% dhdx(:,:,1,i)=-g*dhdx1;
% dhdy(:,:,1,i)=-g*dhdy1;
% end
dhdx(2:end,:,1,:)=-g*(SSH(2:end,:,:)-SSH(1:end-1,:,:))./repmat(DXC(1:end-1,:),[1 1 nt]);
UP=UdPdx+repmat(dhdx,[1 1 46 1]);
dhdy(:,2:end,1,:)=-g*(SSH(:,2:end,:)-SSH(:,1:end-1,:))./repmat(DYC(:,1:end-1),[1 1 nt]);
VP=VdPdy+repmat(dhdy,[1 1 46 1]);
UPress1=squeeze(nansum(nansum(nansum(UP.*repmat(inWag,[1 1 1 nt])))));
VPress1=squeeze(nansum(nansum(nansum(VP.*repmat(inWag,[1 1 1 nt])))));
%UPress2=squeeze(nansum(nansum(-g*repmat(inWag(1:end-1,:,1),[1 1 nt]).*(SSH(2:end,:,:)-SSH(1:end-1,:,:))./repmat(DXC(1:end-1,:),[1 1 nt]))));
%VPress2=squeeze(nansum(nansum(-g*repmat(inWag(:,1:end-1,1),[1 1 nt]).*(SSH(:,2:end,:)-SSH(:,1:end-1,:))./repmat(DYC(:,1:end-1),[1 1 nt]))));

UAdvec=squeeze(nansum(nansum(nansum(AdvU.*repmat(inWag,[1 1 1 nt])))));
VAdvec=squeeze(nansum(nansum(nansum(AdvV.*repmat(inWag,[1 1 1 nt])))));

Ufv=squeeze(nansum(nansum(nansum(UCori.*repmat(inWag,[1 1 1 nt])))));
Vfu=squeeze(nansum(nansum(nansum(VCori.*repmat(inWag,[1 1 1 nt])))));

uSurf=squeeze(nansum(nansum(nansum(Uext.*repmat(inWag,[1 1 1 nt])))));
vSurf=squeeze(nansum(nansum(nansum(Vext.*repmat(inWag,[1 1 1 nt])))));

UDif1=squeeze(nansum(nansum(nansum(UDiss.*repmat(inWag,[1 1 1 nt])))));
VDif1=squeeze(nansum(nansum(nansum(VDiss.*repmat(inWag,[1 1 1 nt])))));
UDif2a=(VisZU(:,:,2:end,:)-VisZU(:,:,1:end-1,:))./cellVolU(:,:,1:end-1,:);
VDif2a=(VisZV(:,:,2:end,:)-VisZV(:,:,1:end-1,:))./cellVolV(:,:,1:end-1,:);
UDif2=squeeze(nansum(nansum(nansum(repmat(inWag(:,:,1:end-1,:),[1 1 1 nt]).*UDif2a))));
VDif2=squeeze(nansum(nansum(nansum(repmat(inWag(:,:,1:end-1,:),[1 1 1 nt]).*VDif2a))));

Utot1=UP(:,:,1:end-1,:)+AdvU(:,:,1:end-1,:)+Uext(:,:,1:end-1,:)+UDiss(:,:,1:end-1,:)+UDif2a+AbU(:,:,1:end-1,:);
Vtot1=VP(:,:,1:end-1,:)+AdvV(:,:,1:end-1,:)+Vext(:,:,1:end-1,:)+VDiss(:,:,1:end-1,:)+VDif2a+AbV(:,:,1:end-1,:);
Utot=squeeze(nansum(nansum(nansum(Utot1.*repmat(inWag(:,:,1:end-1,:),[1 1 1 nt])))));
Vtot=squeeze(nansum(nansum(nansum(Vtot1.*repmat(inWag(:,:,1:end-1,:),[1 1 1 nt])))));

uAB=squeeze(nansum(nansum(nansum(AbU.*repmat(inWag,[1 1 1 nt])))));
vAB=squeeze(nansum(nansum(nansum(AbV.*repmat(inWag,[1 1 1 nt])))));

perrorU=100*(Utot-dUdt)./dUdt;
perrorV=100*(Vtot-dVdt)./dVdt;

Uside=squeeze(nansum(nansum(nansum(Usd.*repmat(inWag,[1 1 1 nt])))));
Vside=squeeze(nansum(nansum(nansum(Vsd.*repmat(inWag,[1 1 1 nt])))));
Ubot=squeeze(nansum(nansum(nansum(Ubd.*repmat(inWag,[1 1 1 nt])))));
Vbot=squeeze(nansum(nansum(nansum(Vbd.*repmat(inWag,[1 1 1 nt])))));

%advUh=squeeze(nansum(nansum(nansum((AdvXU(1:end-1,1:end-1,:,:)-AdvXU(2:end,1:end-1,:,:)+AdvYU(1:end-1,1:end-1,:,:)...
%    -AdvYU(1:end-1,2:end,:,:)).*repmat(inWag(1:end-1,1:end-1,:),[1 1 1 nt])./cellVolU(1:end-1,1:end-1,:,:)))));
%advUz=squeeze(nansum(nansum(nansum((AdvZU(:,:,2:46,:)-AdvZU(:,:,1:45,:)).*repmat(inWag(:,:,1:45),[1 1 1 nt])./cellVolU(:,:,1:45,:)))));

%advVh=squeeze(nansum(nansum(nansum((AdvXV(1:end-1,1:end-1,:,:)-AdvXV(2:end,1:end-1,:,:)+AdvYV(1:end-1,1:end-1,:,:)...
%    -AdvYV(1:end-1,2:end,:,:)).*repmat(inWag(1:end-1,1:end-1,:),[1 1 1 nt])./cellVolV(1:end-1,1:end-1,:,:)))));
%advVz=squeeze(nansum(nansum(nansum((AdvZV(:,:,2:46,:)-AdvZV(:,:,1:45,:)).*repmat(inWag(:,:,1:45),[1 1 1 nt])./cellVolV(:,:,1:45,:)))));


UAdvecVort=squeeze(nansum(nansum(nansum(AdvVortU.*repmat(inWag,[1 1 1 nt])))));
VAdvecVort=squeeze(nansum(nansum(nansum(AdvVortV.*repmat(inWag,[1 1 1 nt])))));
UAdvecZ=squeeze(nansum(nansum(nansum(AdvReU.*repmat(inWag,[1 1 1 nt])))));
VAdvecZ=squeeze(nansum(nansum(nansum(AdvReV.*repmat(inWag,[1 1 1 nt])))));

% need native grid u and v for this!
% dZ2(1,1,1:45)=0.5*dZ(1:45)+0.5*dZ(2:46);
% for i=1:nt
%     dudz=diff(U(:,:,:,i),1,3)./repmat(dZ2,[700 200 1]);
%     UDifK(i)=squeeze(nansum(nansum(nansum((KppV(:,:,1:44,i).*dudz(:,:,1:44)-KppV(:,:,2:45,i).*dudz(:,:,2:45)).*inWag(:,:,1:44)./repmat(dZ(1,1,1:44),[700 200 1])))));
% 
%     dvdz=diff(V(:,:,:,i),1,3)./repmat(dZ2,[700 200 1]);
%     VDifK(i)=squeeze(nansum(nansum(nansum((KppV(:,:,1:44,i).*dvdz(:,:,1:44)-KppV(:,:,2:45,i).*dvdz(:,:,2:45)).*inWag(:,:,1:44)./repmat(dZ(1,1,1:44),[700 200 1])))));
% 
% end


save('momentumBudget14daySplit.mat')
%save('momentumBudget58dayTest2.mat')
%% vertical plot day 1 cells at 300,100
% load('momentumBudget11dayTest.mat','AdvU')
% AdvU=squeeze(AdvU(300,100,:,:));
% load('momentumBudget11dayTest.mat','UP')
% UP=squeeze(UP(300,100,:,:));
% load('momentumBudget11dayTest.mat','Uext')
% Uext=squeeze(Uext(300,100,:,:));
% load('momentumBudget11dayTest.mat','UDiss')
% UDiss=squeeze(UDiss(300,100,:,:));
% load('momentumBudget11dayTest.mat','UDif2a')
% UDif2a=squeeze(UDif2a(300,100,:,:));
% load('momentumBudget11dayTest.mat','Utend')
% Utend=squeeze(Utend(300,100,:,:));
% load('geometrySpinupSteady.mat','dInterface')
% dBin=0.5.*dInterface(1:end-1)+0.5.*dInterface(2:end);
% for i=1:3:11
% figure;
% plot(AdvU(:,i),-dBin)
% hold all
% plot(UP(:,i),-dBin);
% plot(Uext(:,i),-dBin);
% plot(UDiss(:,i),-dBin);
% plot(UDif2a(:,i),-dBin(1:end-1));
% plot(AdvU(1:end-1,i)+UP(1:end-1,i)+Uext(1:end-1,i)+UDiss(1:end-1,i)+UDif2a(:,i),-dBin(1:end-1),'k--','linewidth',2);
% plot(Utend(:,i)./86400,-dBin)
% legend('advection','pressure','wind','dissipation','diffusion','total','dUdt')
% title(strcat('Umom Day ',num2str(i)))
% ylabel('depth (m)')
% xlabel('momentum flux m/s^2')
% ylim([-1000 0])
% end
% 
% i=1;
% figure; 
% plot(nansum(AdvU,i))
% hold all
% plot(nansum(UP,i));
% plot(nansum(Uext,i));
% plot(nansum(UDiss,i));
% plot(nansum(UDif2a,i));
% plot(nansum(AdvU(1:end-1,:)+UP(1:end-1,:)+Uext(1:end-1,:)+UDiss(1:end-1,:)+UDif2a(:,:),i),'k--','linewidth',2);
% plot(nansum(Utend./86400,i))
% legend('advection','pressure','wind','dissipation','diffusion','total','dUdt')
% title('Umom')
% xlabel('day')
% ylabel('momentum flux m/s^2')
% 
% 
% load('momentumBudget11dayTest.mat','AdvV')
% AdvV=squeeze(AdvV(300,100,:,:));
% load('momentumBudget11dayTest.mat','VP')
% VP=squeeze(VP(300,100,:,:));
% load('momentumBudget11dayTest.mat','Vext')
% Vext=squeeze(Vext(300,100,:,:));
% load('momentumBudget11dayTest.mat','VDiss')
% VDiss=squeeze(VDiss(300,100,:,:));
% load('momentumBudget11dayTest.mat','VDif2a')
% VDif2a=squeeze(VDif2a(300,100,:,:));
% load('momentumBudget11dayTest.mat','Vtend')
% Vtend=squeeze(Vtend(300,100,:,:));
% load('geometrySpinupSteady.mat','dInterface')
% dBin=0.5.*dInterface(1:end-1)+0.5.*dInterface(2:end);
% for i=1:3:11
% figure;
% plot(AdvV(:,i),-dBin)
% hold all
% plot(VP(:,i),-dBin);
% plot(Vext(:,i),-dBin);
% plot(VDiss(:,i),-dBin);
% plot(VDif2a(:,i),-dBin(1:end-1));
% plot(AdvV(1:end-1,i)+VP(1:end-1,i)+Vext(1:end-1,i)+VDiss(1:end-1,i)+VDif2a(:,i),-dBin(1:end-1),'k--','linewidth',2);
% plot(Vtend(:,i)./86400,-dBin)
% legend('advection','pressure','wind','dissipation','diffusion','total','dUdt')
% title(strcat('Vmom Day ',num2str(i)))
% ylabel('depth (m)')
% xlabel('momentum flux m/s^2')
% ylim([-1000 0])
% end
% 
% %% plots
% load('momentumBudget11dayTest.mat','*P*','*Surf','*Dif*','*Advec*','d*dt','*tot*','*AB')
% 
% figure; plot(-dUdt,'linewidth',2); hold all
% plot(UPress1); plot(UAdvec); %plot(Ufv);
% plot(uSurf); plot(UDif1+UDif2);
% plot(UAdvecVort); plot(UAdvecZ)
% plot(uAB);
% % plot(-dUdt+UDif1+UDif2+UPress1+UPress2+uSurf+UAdvec-Ufv,'r--','linewidth',2)
% %plot(-dUdt+UPress1+UPress2+uSurf+UDif1+UDif2+UAdvecVort+UAdvecZ,'m--','linewidth',2)
% %plot(-dUdt+UPress1+uSurf+UDif1+UDif2,'k--','linewidth',2)
% %plot(dUdt-UPress1-uSurf-UDif1-UDif2,'g--','linewidth',2)
% plot(-dUdt+(UPress1+uSurf+UDif1+UDif2+UAdvec+uAB),'r--','linewidth',2)
% legend('-dUdt','Pressure term','Advection term','Surface','Dissipation-Diffusion','Advection by Vorticity','Vertical Advection','timestepping','total')
% title('U momentum')
% 
% 
% figure; plot(-dVdt,'linewidth',2); hold all
% plot(VPress1); plot(VAdvec); %plot(Vfu);
% plot(vSurf); plot(VDif1+VDif2);
% plot(VAdvecVort); plot(VAdvecZ)
% plot(vAB);
% plot(-dVdt+VDif1+VDif2+VPress1+vSurf+VAdvecVort+VAdvecZ,'r--','linewidth',2)
% plot(-dVdt+VDif1+VDif2+VPress1+vSurf+VAdvec,'k--','linewidth',2)
% legend('-dVdt','Pressure term','Advection','Surface','Dissipation-Diffusion','Advection by vorticity','vertical advection','timestepping','error using vort+z advec','error using advection term')
% title('V momentum')
%% plots for me
% figure; plot(dUdt,'linewidth',2); hold all
% plot(UPress1);plot(UPress2); %plot(UAdvec/10); %plot(Ufv);
% plot(uSurf); plot(UDif1); plot(UDif2);
% plot(UDif1+UDif2+UPress1+UPress2+uSurf,'r--','linewidth',2)
% plot(dUdt-(UDif1+UDif2+UPress1+UPress2+uSurf),'m','linewidth',2)
% legend('dUdt','Pressure term 1','pressure 2','Surface','Dissipation','Diffusion','sumflux','proper adv')
% 
% figure; plot(dVdt,'linewidth',2); hold all
% plot(VPress1); plot(VPress2); %plot(VAdvec/10); %plot(Vfu);
% plot(vSurf); plot(VDif1);plot(VDif2);
% plot(VDif1+VDif2+VPress1+VPress2+vSurf,'r--','linewidth',2)
% plot(dVdt-(VDif1+VDif2+VPress1+VPress2+vSurf),'m','linewidth',2)
% legend('dVdt','Pressure term 1','pressure 2','Surface','Dissipation','Diffusion','sumflux','proper adv')
% 
% advU1=dUdt-(UDif1+UDif2+UPress1+UPress2+uSurf);
% advV1=dVdt-(VDif1+VDif2+VPress1+VPress2+vSurf);
% 
% figure; plot(UAdvec); hold all; plot(advU1); plot(VAdvec); plot(advV1)
% legend('output u','balance u','output v','balance v')

