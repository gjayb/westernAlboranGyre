%%load variables
 times=8640:8640:(148*8640)
 Utend=rdmds('Utend',times);
 sizeUtend=size(Utend)
 Vtend=rdmds('Vtend',times);
 UDiss=rdmds('UDiss',times);
 VDiss=rdmds('VDiss',times);
 AdvU=rdmds('AdvU',times);
 AdvV=rdmds('AdvV',times);
 disp('saving diagnostics')
 save('momentumDiagnostics148dayNF1.mat','-v7.3')
 clear
 
 times=8640:8640:(148*8640)
 UdPdx=rdmds('UdPdx',times);
 VdPdy=rdmds('VdPdy',times);
 Uext=rdmds('Uext',times);
 Vext=rdmds('Vext',times);
 VisZU=rdmds('ViscZiU',times);
 VisZV=rdmds('ViscZiV',times);
 disp('saving diagnostics')
 save('momentumDiagnostics148dayNF2.mat','-v7.3')
 clear
 
 UCori=rdmds('UCori',times);
 VCori=rdmds('VCori',times);
% 
 times=8640:8640:(148*8640)
 SSH=rdmds('SSHave',times);
 AbU=rdmds('ABU',times);
 AbV=rdmds('ABV',times);
 
% %AdvXU=rdmds('AdvXU',times);%zero
% %AdvYU=rdmds('AdvYU',times);%zero
% %AdvZU=rdmds('AdvZU',times);%zero
% %AdvVortU=rdmds('AdvVortU',times);
% %AdvReU=rdmds('AdvReU',times);%nonzero!
% %AdvXV=rdmds('AdvXV',times);%zero
% %AdvYV=rdmds('AdvYV',times);%zero
% %AdvZV=rdmds('AdvZV',times);%zero
% %AdvVortV=rdmds('AdvVortV',times);%nonzero!
% %AdvReV=rdmds('AdvReV',times);%wasn't saved whoops
% 
  KppK=rdmds('KppVis',times);
 disp('saving diagnostics')
 save('momentumDiagnostics148dayNF3.mat','-v7.3')
 clear
%%

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

nt=148;

dZ(1,1,1:46)=diff(dInterface);
cellVolU=repmat(rAw,[1 1 46]).*hFacW.*repmat(dZ,[700 200 1]);
%cellVolU=repmat(cellVol,[1 1 1 nt]);
cellVolV=repmat(rAs,[1 1 46]).*hFacS.*repmat(dZ,[700 200 1]);
%cellVolV=repmat(cellVol,[1 1 1 nt]);

load('momentumDiagnostics148dayNF1.mat')

sizeUtend=size(Utend)
load('edgesWAGeuler2017NF.mat','inWag','open*','XC','YC')
sizeInWag=size(inWag)
%inWag1=inWag(:,:,1);
%inWag(:,:,1:15)=repmat(inWag1,[1 1 15]);

%removed one nansum to get depth-dependent version of this

dUdt=squeeze(nansum(nansum((Utend.*repmat(inWag.*cellVolU,[1 1 1 nt])./86400))));
dVdt=squeeze(nansum(nansum((Vtend.*repmat(inWag.*cellVolV,[1 1 1 nt])./86400))));
clear Utend Vtend
UDif1=squeeze(nansum(nansum((UDiss.*repmat(inWag.*cellVolU,[1 1 1 nt])))));
VDif1=squeeze(nansum(nansum((VDiss.*repmat(inWag.*cellVolV,[1 1 1 nt])))));
UAdvec=squeeze(nansum(nansum((AdvU.*repmat(inWag.*cellVolU,[1 1 1 nt])))));
VAdvec=squeeze(nansum(nansum((AdvV.*repmat(inWag.*cellVolV,[1 1 1 nt])))));
%Utot1=AdvU(:,:,1:end-1,:)+UDiss(:,:,1:end-1,:);
%Vtot1=AdvV(:,:,1:end-1,:)+VDiss(:,:,1:end-1,:);
clear AdvU AdvV UDiss VDiss
%%


load('momentumDiagnostics148dayNF2.mat','Vis*')
UDif2a=(VisZU(:,:,2:end,:)-VisZU(:,:,1:end-1,:));%./cellVolU(:,:,1:end-1,:);
VDif2a=(VisZV(:,:,2:end,:)-VisZV(:,:,1:end-1,:));%./cellVolV(:,:,1:end-1,:);
clear VisZU VisZV
load('momentumDiagnostics148dayNF2.mat','*ext')
uSurf=squeeze(nansum(nansum((Uext.*repmat(inWag.*cellVolU,[1 1 1 nt])))));
vSurf=squeeze(nansum(nansum((Vext.*repmat(inWag.*cellVolV,[1 1 1 nt])))));
UDif2=squeeze(nansum(nansum((repmat(inWag(:,:,1:end-1,:),[1 1 1 nt]).*UDif2a))));
VDif2=squeeze(nansum(nansum((repmat(inWag(:,:,1:end-1,:),[1 1 1 nt]).*VDif2a))));
%Utot1=Utot1+Uext(:,:,1:end-1,:)+UDif2a./repmat(cellVolU(:,:,1:end-1),[1 1 1 nt]);
%Vtot1=Vtot1+Vext(:,:,1:end-1,:)+VDif2a./repmat(cellVolV(:,:,1:end-1),[1 1 1 nt]);
clear Uext Vext UDif2a VDif2a

load('momentumDiagnostics148dayNF2.mat','*dP*')
load('momentumDiagnostics148dayNF3.mat','SSH')
%need to offset press1 by press2 before inwag and adding up
dhdx=zeros([700 200 1 nt]);
dhdy=dhdx;

% for i=1:nt
% [dhdy1,dhdx1]=gradient(SSH(:,:,i),ycm,xcm);
% dhdx(:,:,1,i)=-g*dhdx1;
% dhdy(:,:,1,i)=-g*dhdy1;
% end
dhdx(2:end,:,1,:)=-g*(SSH(2:end,:,:)-SSH(1:end-1,:,:))./repmat(DXC(1:end-1,:),[1 1 nt]);
UP2=UdPdx+repmat(dhdx,[1 1 46 1]);
sizeUP2=size(UP2)
%dhdx2=dhdx; 
dhdx(2:end-1,:,1,:)=-g*(SSH(3:end,:,:)-SSH(1:end-2,:,:))./repmat(0.5*DXC(1:end-2,:)+DXC(2:end-1,:)+0.5*DXC(3:end,:),[1 1 nt]);
UP=UdPdx+repmat(dhdx,[1 1 46 1]);
sizeUP=size(UP)
%UP2=UdPdx+repmat(dhdx2,[1 1 46 1]);
clear dhdx dhdx2
clear UdPdx 
dhdy(:,2:end,1,:)=-g*(SSH(:,2:end,:)-SSH(:,1:end-1,:))./repmat(DYC(:,1:end-1),[1 1 nt]);
VP2=VdPdy+repmat(dhdy,[1 1 46 1]);
sizeVP2=size(VP2)
%dhdy2=dhdy;
dhdy(:,2:end-1,1,:)=-g*(SSH(:,3:end,:)-SSH(:,1:end-2,:))./repmat(0.5*DYC(:,1:end-2)+DYC(:,2:end-1)+0.5*DYC(:,3:end),[1 1 nt]);
VP=VdPdy+repmat(dhdy,[1 1 46 1]);
sizeVP=size(VP)
clear dhdy
clear VdPdy 
UPress1=squeeze(nansum(nansum((UP.*repmat(inWag.*cellVolU,[1 1 1 nt])))));
VPress1=squeeze(nansum(nansum((VP.*repmat(inWag.*cellVolV,[1 1 1 nt])))));
UPress2=squeeze(nansum(nansum((UP2.*repmat(inWag.*cellVolU,[1 1 1 nt])))));
VPress2=squeeze(nansum(nansum((VP2.*repmat(inWag.*cellVolV,[1 1 1 nt])))));
%UPress2=squeeze(nansum(nansum(-g*repmat(inWag(1:end-1,:,1),[1 1 nt]).*(SSH(2:end,:,:)-SSH(1:end-1,:,:))./repmat(DXC(1:end-1,:),[1 1 nt]))));
%VPress2=squeeze(nansum(nansum(-g*repmat(inWag(:,1:end-1,1),[1 1 nt]).*(SSH(:,2:end,:)-SSH(:,1:end-1,:))./repmat(DYC(:,1:end-1),[1 1 nt]))));
%Ufv=squeeze(nansum(nansum(nansum(UCori.*repmat(inWag,[1 1 1 nt])))));
%Vfu=squeeze(nansum(nansum(nansum(VCori.*repmat(inWag,[1 1 1 nt])))));

load('momentumDiagnostics148dayNF3.mat','Ab*')
%Utot1=Utot1+UP(:,:,1:end-1,:)+AbU(:,:,1:end-1,:);
%Vtot1=Vtot1+VP(:,:,1:end-1,:)+AbV(:,:,1:end-1,:);
%Utot=squeeze(nansum(nansum(nansum(Utot1.*repmat(inWag(:,:,1:end-1,:).*cellVolU(:,:,1:end-1,:),[1 1 1 nt])))));
%Vtot=squeeze(nansum(nansum(nansum(Vtot1.*repmat(inWag(:,:,1:end-1,:).*cellVolV(:,:,1:end-1,:),[1 1 1 nt])))));
uAB=squeeze(nansum(nansum((AbU.*repmat(inWag.*cellVolU,[1 1 1 nt])))));
vAB=squeeze(nansum(nansum((AbV.*repmat(inWag.*cellVolV,[1 1 1 nt])))));
clear AbU AbV %Utot1 Vtot1
sizeUP=size(UP)

perrorU=100*(Utot-dUdt)./dUdt;
perrorV=100*(Vtot-dVdt)./dVdt;



%advUh=squeeze(nansum(nansum(nansum((AdvXU(1:end-1,1:end-1,:,:)-AdvXU(2:end,1:end-1,:,:)+AdvYU(1:end-1,1:end-1,:,:)...
%    -AdvYU(1:end-1,2:end,:,:)).*repmat(inWag(1:end-1,1:end-1,:),[1 1 1 nt])./cellVolU(1:end-1,1:end-1,:,:)))));
%advUz=squeeze(nansum(nansum(nansum((AdvZU(:,:,2:46,:)-AdvZU(:,:,1:45,:)).*repmat(inWag(:,:,1:45),[1 1 1 nt])./cellVolU(:,:,1:45,:)))));

%advVh=squeeze(nansum(nansum(nansum((AdvXV(1:end-1,1:end-1,:,:)-AdvXV(2:end,1:end-1,:,:)+AdvYV(1:end-1,1:end-1,:,:)...
%    -AdvYV(1:end-1,2:end,:,:)).*repmat(inWag(1:end-1,1:end-1,:),[1 1 1 nt])./cellVolV(1:end-1,1:end-1,:,:)))));
%advVz=squeeze(nansum(nansum(nansum((AdvZV(:,:,2:46,:)-AdvZV(:,:,1:45,:)).*repmat(inWag(:,:,1:45),[1 1 1 nt])./cellVolV(:,:,1:45,:)))));


% UAdvecVort=squeeze(nansum(nansum(nansum(AdvVortU.*repmat(inWag,[1 1 1 nt])))));
% VAdvecVort=squeeze(nansum(nansum(nansum(AdvVortV.*repmat(inWag,[1 1 1 nt])))));
% UAdvecZ=squeeze(nansum(nansum(nansum(AdvReU.*repmat(inWag,[1 1 1 nt])))));
% VAdvecZ=squeeze(nansum(nansum(nansum(AdvReV.*repmat(inWag,[1 1 1 nt])))));

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
load('momentumCori');
uCori=squeeze(nansum(nansum((UCori.*repmat(inWag.*cellVolU,[1 1 1 nt])))));
vCori=squeeze(nansum(nansum((VCori.*repmat(inWag.*cellVolV,[1 1 1 nt])))));

VPress3=dVdt-VAdvec-VDif1-VDif2-vAB-vSurf;
UPress3=dUdt-UAdvec-UDif1-UDif2-uAB-uSurf;

save('momentumBudget148dayNFdepth.mat')
%save('momentumBudget58dayTest2.mat')
%% vertical plot day 1 cells at 300,100
sizeUP=size(UP)
runthis=false;
if runthis
xa=350; ya=70;
UP=squeeze(UP(xa,ya,:,:));
VP=squeeze(VP(xa,ya,:,:));
UP2=squeeze(UP2(xa,ya,:,:));
VP2=squeeze(VP2(xa,ya,:,:));
load('momentumDiagnostics148dayNF1.mat','AdvU')
AdvU=squeeze(AdvU(xa,ya,:,:));
load('momentumDiagnostics148dayNF2.mat','Uext')
Uext=squeeze(Uext(xa,ya,:,:));
load('momentumDiagnostics148dayNF1.mat','UDiss')
UDiss=squeeze(UDiss(xa,ya,:,:));
load('momentumDiagnostics148dayNF2.mat','VisZU')
UDif2a=(VisZU(:,:,2:end,:)-VisZU(:,:,1:end-1,:))./repmat(cellVolU(:,:,1:end-1),[1 1 1 nt]);
UDif2a=squeeze(UDif2a(xa,ya,:,:)); clear VisZU
load('momentumDiagnostics148dayNF1.mat','Utend')
Utend=squeeze(Utend(xa,ya,:,:));
load('geometrySpinupSteady.mat','dInterface')
dBin=0.5.*dInterface(1:end-1)+0.5.*dInterface(2:end);
% for i=1:3:11
% figure;
% plot(AdvU(:,i),-dBin)
% hold all
% plot(UP(:,i),-dBin);
% plot(Uext(:,i),-dBin);
% plot(UDiss(:,i),-dBin);
% plot(UDif2a(:,i),-dBin(1:end-1));
% plot(-Utend(:,i)./86400,-dBin);
% plot(-Utend(1:end-1,i)./86400+AdvU(1:end-1,i)+UP(1:end-1,i)+Uext(1:end-1,i)+UDiss(1:end-1,i)+UDif2a(:,i),-dBin(1:end-1),'k--','linewidth',2);
% legend('advection','pressure','wind','dissipation','diffusion','-dUdt','total')
% title(strcat('Umom Day ',num2str(i)))
% ylabel('depth (m)')
% xlabel('momentum flux m/s^2')
% ylim([-1000 0])
% end

% i=1;
% figure; 
% plot(nansum(AdvU,i))
% hold all
% plot(nansum(UP,i));
% plot(nansum(Uext,i));
% plot(nansum(UDiss,i));
% plot(nansum(UDif2a,i));
% plot(nansum(-Utend./86400,i))
% plot(nansum(-Utend(1:end-1,:)./86400+AdvU(1:end-1,:)+UP(1:end-1,:)+Uext(1:end-1,:)+UDiss(1:end-1,:)+UDif2a(:,:),i),'k--','linewidth',2);
% legend('advection','pressure','wind','dissipation','diffusion','-dUdt','total')
% title('Umom')
% xlabel('day')
% ylabel('momentum flux m/s^2')


load('momentumDiagnostics148dayNF1.mat','AdvV')
AdvV=squeeze(AdvV(xa,ya,:,:));
load('momentumDiagnostics148dayNF2.mat','Vext')
Vext=squeeze(Vext(xa,ya,:,:));
load('momentumDiagnostics148dayNF2.mat','VisZV')
VDif2a=(VisZV(:,:,2:end,:)-VisZV(:,:,1:end-1,:))./repmat(cellVolV(:,:,1:end-1),[1 1 1 nt]);
VDif2a=squeeze(VDif2a(xa,ya,:,:)); clear VisZV cellVolV
load('momentumDiagnostics148dayNF1.mat','VDiss')
VDiss=squeeze(VDiss(xa,ya,:,:));
load('momentumDiagnostics148dayNF1.mat','Vtend')
Vtend=squeeze(Vtend(xa,ya,:,:));
load('geometrySpinupSteady.mat','dInterface')
dBin=0.5.*dInterface(1:end-1)+0.5.*dInterface(2:end);
% for i=1:3:11
% figure;
% plot(AdvV(:,i),-dBin)
% hold all
% plot(VP(:,i),-dBin);
% plot(Vext(:,i),-dBin);
% plot(VDiss(:,i),-dBin);
% plot(VDif2a(:,i),-dBin(1:end-1));
% plot(-Vtend(:,i)./86400,-dBin);
% plot(-Vtend(1:end-1,i)./86400+AdvV(1:end-1,i)+VP(1:end-1,i)+Vext(1:end-1,i)+VDiss(1:end-1,i)+VDif2a(:,i),-dBin(1:end-1),'k--','linewidth',2);
% legend('advection','pressure','wind','dissipation','diffusion','-dUdt','total')
% title(strcat('Vmom Day ',num2str(i)))
% ylabel('depth (m)')
% xlabel('momentum flux m/s^2')
% ylim([-1000 0])
% end
% 
% i=1;
% figure; 
% plot(nansum(AdvV,i))
% hold all
% plot(nansum(VP,i));
% plot(nansum(Vext,i));
% plot(nansum(VDiss,i));
% plot(nansum(VDif2a,i));
% plot(nansum(-Vtend./86400,i))
% plot(nansum(-Vtend(1:end-1,:)./86400+AdvV(1:end-1,:)+VP(1:end-1,:)+Vext(1:end-1,:)+VDiss(1:end-1,:)+VDif2a(:,:),i),'k--','linewidth',2);
% legend('advection','pressure','wind','dissipation','diffusion','-dVdt','total')
% title('Vmom')
% xlabel('day')
% ylabel('momentum flux m/s^2')
save('momentumBudgetOneColumnNF.mat','*V*','*U*','dBin')
end
%% plots
% load('momentumBudget11dayEulerWag.mat','*P*','*Surf','*Dif*','*Advec*','d*dt','*tot*','*AB')

% figure; plot(-dUdt,'linewidth',2); hold all
% plot(UPress1); plot(UAdvec); %plot(Ufv);
% plot(uSurf); plot(UDif1+UDif2);
% %plot(UAdvecVort); plot(UAdvecZ)
% plot(uAB);
% % plot(-dUdt+UDif1+UDif2+UPress1+UPress2+uSurf+UAdvec-Ufv,'r--','linewidth',2)
% %plot(-dUdt+UPress1+UPress2+uSurf+UDif1+UDif2+UAdvecVort+UAdvecZ,'m--','linewidth',2)
% %plot(-dUdt+UPress1+uSurf+UDif1+UDif2,'k--','linewidth',2)
% %plot(dUdt-UPress1-uSurf-UDif1-UDif2,'g--','linewidth',2)
% plot(-dUdt+(UPress1+uSurf+UDif1+UDif2+UAdvec+uAB),'r--','linewidth',2)
% plot(Utot-dUdt,'k--','linewidth',2)
% legend('-dUdt','Pressure term','Advection term','Surface','Dissipation-Diffusion','Advection by Vorticity','Vertical Advection','timestepping','total','total')
% title('U momentum')
% %Utot=nansum(nansum(UP(:,:,1:end-1,:)+AdvU(:,:,1:end-1,:)+Uext(:,:,1:end-1,:)+UDiss(:,:,1:end-1,:)+UDif2a+AbU(:,:,1:end-1,:)));
% 
% figure; plot(-dVdt,'linewidth',2); hold all
% plot(VPress1); plot(VAdvec); %plot(Vfu);
% plot(vSurf); plot(VDif1+VDif2);
% %plot(VAdvecVort); plot(VAdvecZ)
% plot(vAB);
% %plot(-dVdt+VDif1+VDif2+VPress1+vSurf+VAdvecVort+VAdvecZ,'r--','linewidth',2)
% %plot(-dVdt+VDif1+VDif2+VPress1+vSurf+VAdvec,'k--','linewidth',2)
% plot(Vtot-dVdt,'k--','linewidth',2)
% legend('-dVdt','Pressure term','Advection','Surface','Dissipation-Diffusion','Advection by vorticity','vertical advection','timestepping','total')
% title('V momentum')
%% plots for me
% figure; plot(dUdt,'linewidth',2); hold all
% plot(UAdvec+UPress1);
% %plot(UPress0); %plot(UPress2); %plot(UAdvec/10); %plot(Ufv);
% plot(uSurf); plot(UDif1); plot(UDif2); plot(uAB)
% plot(UAdvec+UDif1+UDif2+UPress1+uSurf+uAB,'r--','linewidth',2)
% plot(-dUdt+(UAdvec+UDif1+UDif2+UPress1+uAB+uSurf),'m','linewidth',2)
% legend('dUdt','advec+Pressure term','Surface','Dissipation','Diffusion','timestep','sumflux','err')
% 
% 
% verr=-dVdt+(VAdvec+VDif1+VDif2+VPress1+vAB+vSurf);
% figure; plot(dVdt,'linewidth',2); hold all
% plot(VAdvec+VPress1); %plot(VPress2); %plot(VAdvec/10); %plot(Vfu);
% plot(vSurf); plot(VDif1);plot(VDif2);plot(vAB);
% plot(VAdvec+VDif1+VDif2+VPress1+vAB+vSurf,'r--','linewidth',2)
% plot(-dVdt+VAdvec+VDif1+VDif2+VPress1+vAB+vSurf,'m','linewidth',2)
% legend('dVdt','advec+Pressure term 1','Surface','Dissipation','Diffusion','timestep','sumflux','err')
% 
% % advU1=dUdt-(UDif1+UDif2+UPress1+UPress2+uSurf);
% % advV1=dVdt-(VDif1+VDif2+VPress1+VPress2+vSurf);
% 
% %figure; plot(UAdvec); hold all; plot(advU1); plot(VAdvec); plot(advV1)
% %legend('output u','balance u','output v','balance v')

%%
% load('momentumCori');
% for i=1:148
% uCori(i)=squeeze(nansum(nansum(nansum(UCori(:,:,:,i).*double(inWag).*cellVolU))));
% end
% clear UCori
% 
% for i=1:148
% vCori(i)=squeeze(nansum(nansum(nansum(VCori(:,:,:,i).*double(inWag).*cellVolV))));
% end
% clear VCori
%%
% %VPress2-dVdt+VAdvec+VDif1+VDif2+vAB+vSurf=0
% VPress2=dVdt-VAdvec-VDif1-VDif2-vAB-vSurf;
% figure; plot(-dUdt,'linewidth',2); hold all
% plot(UAdvec-uCori.','linewidth',2); plot(uCori,'linewidth',2); plot(UPress1,'linewidth',2);
% plot(uSurf,'linewidth',2); plot(UDif1,'linewidth',2); plot(UDif2,'linewidth',2); plot(uAB,'linewidth',2)
% %plot(UAdvec+UDif1+UDif2+UPress1+uSurf+uAB,'r--','linewidth',2)
% plot(-dUdt+(UAdvec+UDif1+UDif2+UPress1+uAB+uSurf),'m--','linewidth',2)
% legend('-dU/dt','Advection','Coriolis','Pressure term','Surface','Dissipation','Diffusion','Timestep','Total')
% title('U momentum','fontsize',14)
% set(gca,'fontsize',12)
% xlabel('simulation day')
% ylabel('volume-integrated momentum budget terms, m^4/s^2')
% axis tight
% 
% figure; plot(-dVdt,'linewidth',2); hold all
% plot(VAdvec-vCori.','linewidth',2); plot(vCori,'linewidth',2); plot(VPress2,'linewidth',2); 
% plot(vSurf,'linewidth',2); plot(VDif1,'linewidth',2);plot(VDif2,'linewidth',2);plot(vAB,'linewidth',2);
% %plot(VAdvec+VDif1+VDif2+VPress1+vAB+vSurf,'r--','linewidth',2)
% plot(-dVdt+VAdvec+VDif1+VDif2+VPress2+vAB+vSurf,'m--','linewidth',2)
% legend('dV/dt','Advection','Coriolis','Pressure term 1','Surface','Dissipation','Diffusion','Timestep','Total')
% title('V momentum','fontsize',14)
% set(gca,'fontsize',12)
% xlabel('simulation day')
% ylabel('volume-integrated momentum budget terms, m^4/s^2')
% axis tight
% 
% UPress2=dUdt-UAdvec-UDif1-UDif2-uAB-uSurf;
% figure; plot(dUdt,'linewidth',2); hold all
% plot(uCori.'+UPress2,'linewidth',2); plot(UAdvec-uCori.','linewidth',2);
% %plot(UPress0); %plot(UPress2); %plot(UAdvec/10); %plot(Ufv);
% plot(uSurf,'linewidth',2); plot(UDif1,'linewidth',2); %plot(UDif2); plot(uAB)
% %plot(UAdvec+UDif1+UDif2+UPress1+uSurf+uAB,'r--','linewidth',2)
% %plot(-dUdt+(UAdvec+UDif1+UDif2+UPress1+uAB+uSurf),'m--','linewidth',2)
% legend('dU/dt','Ageostrophy','Advection','Surface','Dissipation')%,'Diffusion','Timestep','Total')
% title('U momentum','fontsize',14)
% set(gca,'fontsize',12)
% xlabel('simulation day')
% ylabel('volume-integrated momentum budget terms, m^4/s^2')
% axis tight
% 
% figure; plot(dVdt,'linewidth',2); hold all
% plot(vCori.'+VPress2,'linewidth',2); %plot(VPress2); %plot(VAdvec/10); %plot(Vfu);
% plot(VAdvec-vCori.','linewidth',2);
% plot(vSurf,'linewidth',2); plot(VDif1,'linewidth',2);%plot(VDif2);plot(vAB);
% %plot(VAdvec+VDif1+VDif2+VPress1+vAB+vSurf,'r--','linewidth',2)
% %plot(-dVdt+VAdvec+VDif1+VDif2+VPress2+vAB+vSurf,'m--','linewidth',2)
% legend('dV/dt','Ageostrophy','Advection','Surface','Dissipation')
% title('V momentum','fontsize',14)
% set(gca,'fontsize',12)
% xlabel('simulation day')
% ylabel('volume-integrated momentum budget terms, m^4/s^2')
% axis tight



