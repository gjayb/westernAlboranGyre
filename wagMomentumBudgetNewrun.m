% %load variables
% times=8640:8640:(58*8640)
% Utend=rdmds('Utend',times);
% Vtend=rdmds('Vtend',times);
% UDiss=rdmds('UDiss',times);
% VDiss=rdmds('VDiss',times);
% AdvU=rdmds('AdvU',times);
% AdvV=rdmds('AdvV',times);
% UCori=rdmds('UCori',times);
% VCori=rdmds('VCori',times);
% UdPdx=rdmds('UdPdx',times);
% VdPdy=rdmds('VdPdy',times);
% Uext=rdmds('Uext',times);
% Vext=rdmds('Vext',times);
% VisZU=rdmds('ViscZiU',times);
% VisZV=rdmds('ViscZiV',times);
% 
% SSH=rdmds('SSHave',times);
% 
% save('momentumDiagnostics58day.mat','-v7.3')

%%
load('momentumDiagnostics58day.mat')
%SSH=rdmds('SSHave',times);
%sizeSSH=size(SSH)
%save('momentumDiagnostics58day.mat','-v7.3')
sizeUtend=size(Utend)
load('edgesWAGeuler2017b.mat','inWag','open*','dInterface','XC','YC','*Coast')
sizeInWag=size(inWag)
inWag=zeros(size(inWag));
%inWag(300,100,5)=1;
inWag(300,100,1:16)=1;
load('distancesAreas.mat','ra*','hFac*','dxc','dyc')
DXC=reshape(dxc,[700 200]); DYC=reshape(dyc,[700 200]);
rAw=reshape(raw,[700 200]); rAs=reshape(ras,[700 200]);
load('varying148ts16levelsRho.mat', 'Rho')
Rho=squeeze(Rho(:,:,1,1:58));
load('cp16levels148.mat', 'cp')
cp=squeeze(cp(:,:,1,1:58));
g=9.81;

nt=58;

dZ(1,1,1:46)=diff(dInterface);
cellVol=repmat(rAw,[1 1 46]).*hFacW.*repmat(dZ,[700 200 1]);
cellVolU=repmat(cellVol,[1 1 1 nt]);
cellVol=repmat(rAs,[1 1 46]).*hFacS.*repmat(dZ,[700 200 1]);
cellVolV=repmat(cellVol,[1 1 1 nt]);

%%
dUdt=squeeze(nansum(nansum(nansum(Utend.*repmat(inWag,[1 1 1 nt])./86400))));
dVdt=squeeze(nansum(nansum(nansum(Vtend.*repmat(inWag,[1 1 1 nt])./86400))));

UPress1=squeeze(nansum(nansum(nansum(UdPdx.*repmat(inWag,[1 1 1 nt])))));
VPress1=squeeze(nansum(nansum(nansum(VdPdy.*repmat(inWag,[1 1 1 nt])))));
UPress2=squeeze(nansum(nansum(-g*repmat(inWag(1:end-1,:,1),[1 1 nt]).*(SSH(2:end,:,:)-SSH(1:end-1,:,:))./repmat(DXC(1:end-1,:),[1 1 nt]))));
VPress2=squeeze(nansum(nansum(-g*repmat(inWag(:,1:end-1,1),[1 1 nt]).*(SSH(:,2:end,:)-SSH(:,1:end-1,:))./repmat(DYC(:,1:end-1),[1 1 nt]))));

UAdvec=squeeze(nansum(nansum(nansum(AdvU.*repmat(inWag,[1 1 1 nt])))));
VAdvec=squeeze(nansum(nansum(nansum(AdvV.*repmat(inWag,[1 1 1 nt])))));

Ufv=squeeze(nansum(nansum(nansum(UCori.*repmat(inWag,[1 1 1 nt])))));
Vfu=squeeze(nansum(nansum(nansum(VCori.*repmat(inWag,[1 1 1 nt])))));

uSurf=squeeze(nansum(nansum(nansum(Uext.*repmat(inWag,[1 1 1 nt])))));
vSurf=squeeze(nansum(nansum(nansum(Vext.*repmat(inWag,[1 1 1 nt])))));

UDif1=squeeze(nansum(nansum(nansum(UDiss.*repmat(inWag,[1 1 1 nt])))));
VDif1=squeeze(nansum(nansum(nansum(VDiss.*repmat(inWag,[1 1 1 nt])))));
UDif2=squeeze(nansum(nansum(nansum(repmat(inWag(:,:,1:end-1,:),[1 1 1 nt]).*(VisZU(:,:,2:end,:)-VisZU(:,:,1:end-1,:))./cellVolU(:,:,1:end-1,:)))));
VDif2=squeeze(nansum(nansum(nansum(repmat(inWag(:,:,1:end-1,:),[1 1 1 nt]).*(VisZV(:,:,2:end,:)-VisZV(:,:,1:end-1,:))./cellVolV(:,:,1:end-1,:)))));

Utot=UPress1+UPress2+UAdvec+Ufv+uSurf+UDif1+UDif2;
Vtot=VPress1+VPress2+VAdvec+Vfu+vSurf+VDif1+VDif2;

perrorU=100*(Utot-dUdt)./dUdt;
perrorV=100*(Vtot-dVdt)./dVdt;

%save('momentumBudget58dayTest.mat')
save('momentumBudget58dayTest2.mat')
%% plots
% figure; plot(-dUdt,'linewidth',2); hold all
% plot(UPress1+UPress2); plot(UAdvec); plot(Ufv);
% plot(uSurf); plot(UDif1+UDif2);
% plot(-dUdt+UDif1+UDif2+UPress1+UPress2+uSurf+UAdvec-Ufv,'r--','linewidth',2)
% legend('-dUdt','Pressure term','Advection','Coriolis','Surface','Dissipation-Diffusion','error')
% 
% figure; plot(-dVdt,'linewidth',2); hold all
% plot(VPress1+VPress2); plot(VAdvec); plot(Vfu);
% plot(vSurf); plot(VDif1+VDif2);
% plot(-dVdt+VDif1+VDif2+VPress1+VPress2+vSurf+VAdvec-Vfu,'r--','linewidth',2)
% legend('-dVdt','Pressure term','Advection','Coriolis','Surface','Dissipation-Diffusion','error')
