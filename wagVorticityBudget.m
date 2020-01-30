%based on wagMomentumBudget, eulerian vorticity budget, circulation
%integral- circulation is int(u*ds)=int(nabla cross u)dA=int(vort)dA
%surface only
%% load
% times=(110*8640):8640:(148*8640)
% nt=length(times)
% Utend=rdmds('Utend',times);
% sizeUtend=size(Utend)
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
% AbU=rdmds('ABU',times);
% AbV=rdmds('ABV',times);
% Usnap=rdmds('U',times);
% Vsnap=rdmds('V',times);
% times2=1334880:8640:1395360; %vort
% vortSnap=rdmds('VORT',times2);
% TauEa=rdmds('TAUXave',times);
% TauNa=rdmds('TAUYave',times);
% save('vorticityDiagnostics8day.mat')
%%
%load('vorticityDiagnostics8day.mat')
%vort=rdmds('VortAve',times);

%sizeUtend=size(Utend)
load('edgesWAGeuler2017b.mat','inWag','dInterface','XC','YC','*Coast')
inWag1=inWag(:,:,1);


load('distancesAreas.mat','ra*','hFac*','dxc','dyc')
DXC=reshape(dxc,[700 200]); DYC=reshape(dyc,[700 200]);
rAw=reshape(raw,[700 200]); rAs=reshape(ras,[700 200]);
load('geometrySpinupSteady.mat','XC','YC','dInterface')
xmin=min(min(XC));
ymin=min(min(YC));
xcm=(XC-xmin).*111000.*cosd(YC);
ycm=(YC-ymin).*111000;
g=9.81;

nt=8;

dZ(1,1,1:46)=diff(dInterface);
cellVol=repmat(rAw,[1 1 46]).*hFacW.*repmat(dZ,[700 200 1]);
cellVolU=repmat(cellVol,[1 1 1 nt]);
cellVol=repmat(rAs,[1 1 46]).*hFacS.*repmat(dZ,[700 200 1]);
cellVolV=repmat(cellVol,[1 1 1 nt]);
%% prep variables
%curl of gradient is zero, so don't add surface gradient in here for
%vorticity terms, only for momentum terms
dhdx=zeros([700 200 1 nt]);
dhdy=dhdx;
dhdx(2:end,:,1,:)=-g*(SSH(2:end,:,:)-SSH(1:end-1,:,:))./repmat(DXC(1:end-1,:),[1 1 nt]);
UP=UdPdx+repmat(dhdx,[1 1 46 1]);
dhdy(:,2:end,1,:)=-g*(SSH(:,2:end,:)-SSH(:,1:end-1,:))./repmat(DYC(:,1:end-1),[1 1 nt]);
VP=VdPdy+repmat(dhdy,[1 1 46 1]);
for i=1:8
curlgradH(:,:,i)=squeeze(DYC(2:end,2:end).*(dhdy(2:end,2:end,1,i)-dhdy(1:end-1,2:end,1,i))-DXC(2:end,2:end).*(dhdx(2:end,2:end,1,i)-dhdx(2:end,1:end-1,1,i)))./rAz(2:end,2:end);
end

UDif2a=(VisZU(:,:,2:end,:)-VisZU(:,:,1:end-1,:))./cellVolU(:,:,1:end-1,:);
VDif2a=(VisZV(:,:,2:end,:)-VisZV(:,:,1:end-1,:))./cellVolV(:,:,1:end-1,:);

Utot1=UP(:,:,1:end-1,:)+AdvU(:,:,1:end-1,:)+Uext(:,:,1:end-1,:)+UDiss(:,:,1:end-1,:)+UDif2a+AbU(:,:,1:end-1,:);
Vtot1=VP(:,:,1:end-1,:)+AdvV(:,:,1:end-1,:)+Vext(:,:,1:end-1,:)+VDiss(:,:,1:end-1,:)+VDif2a+AbV(:,:,1:end-1,:);

%advection AdvU AdvV, coriolis in advection but have UCori VCori
%dissipation UDiss VDiss, wind Uext Vext, timestep AbU AbV
%pressure UP VP, diffusion UDifa VDifa

AdvU1=AdvU-UCori;
AdvV1=AdvV-VCori;
Utend=Utend./86400; 
Vtend=Vtend./86400;

clear Vis* cellVol*
%% taking the curl mitgcm style
load('distancesAreas.mat','raz','dxc','dyc')
rAz=reshape(raz,[700 200]);
DXC=reshape(dxc,[700 200]);
DYC=reshape(dyc,[700 200]);
for i=1:nt
hold1=(DYC(2:end,2:end).*(AdvV(2:end,2:end,1,i)-AdvV(1:end-1,2:end,1,i))-DXC(2:end,2:end).*(AdvU(2:end,2:end,1,i)-AdvU(2:end,1:end-1,1,i)))./rAz(2:end,2:end);
cAdv(:,:,i)=hold1;
hold1=(DYC(2:end,2:end).*(AdvV1(2:end,2:end,1,i)-AdvV1(1:end-1,2:end,1,i))-DXC(2:end,2:end).*(AdvU1(2:end,2:end,1,i)-AdvU1(2:end,1:end-1,1,i)))./rAz(2:end,2:end);
cAdvOnly(:,:,i)=hold1;
hold1=(DYC(2:end,2:end).*(VCori(2:end,2:end,1,i)-VCori(1:end-1,2:end,1,i))-DXC(2:end,2:end).*(UCori(2:end,2:end,1,i)-UCori(2:end,1:end-1,1,i)))./rAz(2:end,2:end);
cCori(:,:,i)=hold1;

hold1=(DYC(2:end,2:end).*(VDiss(2:end,2:end,1,i)-VDiss(1:end-1,2:end,1,i))-DXC(2:end,2:end).*(UDiss(2:end,2:end,1,i)-UDiss(2:end,1:end-1,1,i)))./rAz(2:end,2:end);
cDiss(:,:,i)=hold1;
hold1=(DYC(2:end,2:end).*(VDif2a(2:end,2:end,1,i)-VDif2a(1:end-1,2:end,1,i))-DXC(2:end,2:end).*(UDif2a(2:end,2:end,1,i)-UDif2a(2:end,1:end-1,1,i)))./rAz(2:end,2:end);
cDiff(:,:,i)=hold1;
hold1=(DYC(2:end,2:end).*(VdPdy(2:end,2:end,1,i)-VdPdy(1:end-1,2:end,1,i))-DXC(2:end,2:end).*(UdPdx(2:end,2:end,1,i)-UdPdx(2:end,1:end-1,1,i)))./rAz(2:end,2:end);
cPclinic(:,:,i)=hold1;

hold1=(DYC(2:end,2:end).*(VP(2:end,2:end,1,i)-VP(1:end-1,2:end,1,i))-DXC(2:end,2:end).*(UP(2:end,2:end,1,i)-UP(2:end,1:end-1,1,i)))./rAz(2:end,2:end);
cP(:,:,i)=hold1;


hold1=(DYC(2:end,2:end).*(Vext(2:end,2:end,1,i)-Vext(1:end-1,2:end,1,i))-DXC(2:end,2:end).*(Uext(2:end,2:end,1,i)-Uext(2:end,1:end-1,1,i)))./rAz(2:end,2:end);
cWind(:,:,i)=hold1;
hold1=(DYC(2:end,2:end).*(AbV(2:end,2:end,1,i)-AbV(1:end-1,2:end,1,i))-DXC(2:end,2:end).*(AbU(2:end,2:end,1,i)-AbU(2:end,1:end-1,1,i)))./rAz(2:end,2:end);
cAb(:,:,i)=hold1;
hold1=(DYC(2:end,2:end).*(Vtend(2:end,2:end,1,i)-Vtend(1:end-1,2:end,1,i))-DXC(2:end,2:end).*(Utend(2:end,2:end,1,i)-Utend(2:end,1:end-1,1,i)))./rAz(2:end,2:end);
cTend(:,:,i)=hold1;
hold1=(DYC(2:end,2:end).*(Vtot1(2:end,2:end,1,i)-Vtot1(1:end-1,2:end,1,i))-DXC(2:end,2:end).*(Utot1(2:end,2:end,1,i)-Utot1(2:end,1:end-1,1,i)))./rAz(2:end,2:end);
cTot1(:,:,i)=hold1;

hold1=(DYC(2:end,2:end).*(Vsnap(2:end,2:end,1,i)-Vsnap(1:end-1,2:end,1,i))-DXC(2:end,2:end).*(Usnap(2:end,2:end,1,i)-Usnap(2:end,1:end-1,1,i)))./rAz(2:end,2:end);
cVort(:,:,i)=hold1;

if i>1
vortTend(:,:,i)=(vortSnap(:,:,1,i)-vortSnap(:,:,1,i-1))./86400;
end
end

cPtropic=-cPclinic+cP;
%vortSnap=squeeze(vortSnap(:,:,1,:));
%tauN=squeeze(TauNa(:,:,1,:));
%tauE=squeeze(TauEa(:,:,1,:));
save('surfaceVorticity8day.mat','c*','inWag1','XC','x*','YC','y*','*Coast','DXC','DYC','rAz','vortTend','vortSnap','TauNa','TauEa')
%% plot 1
load('surfaceVorticity8day.mat')
i=350; j=70;
figure; plot(squeeze(cAdv(i,j,:))); hold all
plot(squeeze(cDiss(i,j,:))); plot(squeeze(cWind(i,j,:))); plot(squeeze(cDiff(i,j,:))); plot(squeeze(cP(i,j,:))); 
plot(squeeze(cAb(i,j,:)))
plot(squeeze(cTend(i,j,:)),'linewidth',2); plot(squeeze(cTot1(i,j,:)),'--','linewidth',2); plot(squeeze(vortTend(i,j,:)),'--')
plot(squeeze(cAdv(i,j,:)+cDiss(i,j,:)+cWind(i,j,:)+cDiff(i,j,:)+cPtropic(i,j,:)+cPclinic(i,j,:)+cAb(i,j,:)-cTend(i,j,:)),'linewidth',2)
plot(squeeze(cPtropic(i,j,:))); plot(squeeze(cPclinic(i,j,:)));
legend('adv','diss','wind','diff','timestep','tend','tot1','vortTend','tot2','Tropic','Clinic')
title('1 point')

%these should probably be multiplied by rAz, right?
figure; plot(squeeze(nansum(nansum(cAdv.*repmat(rAz(2:end,2:end).*inWag1(2:end,2:end),[1 1 nt]))))); hold all
plot(squeeze(nansum(nansum(cDiss.*repmat(rAz(2:end,2:end).*inWag1(2:end,2:end),[1 1 nt]))))); plot(squeeze(nansum(nansum(cWind.*repmat(rAz(2:end,2:end).*inWag1(2:end,2:end),[1 1 nt]))))); 
plot(squeeze(nansum(nansum(cDiff.*repmat(rAz(2:end,2:end).*inWag1(2:end,2:end),[1 1 nt]))))); plot(squeeze(nansum(nansum(cP.*repmat(rAz(2:end,2:end).*inWag1(2:end,2:end),[1 1 nt]))))); plot(squeeze(nansum(nansum(cAb.*repmat(rAz(2:end,2:end).*inWag1(2:end,2:end),[1 1 nt])))));
plot(squeeze(nansum(nansum(cTend.*repmat(rAz(2:end,2:end).*inWag1(2:end,2:end),[1 1 nt])))),'linewidth',2); plot(squeeze(nansum(nansum(cTot1.*repmat(rAz(2:end,2:end).*inWag1(2:end,2:end),[1 1 nt])))),'--','linewidth',2); plot(squeeze(nansum(nansum(vortTend.*repmat(rAz.*inWag1,[1 1 nt])))),'--')
plot(squeeze(nansum(nansum((cAdv+cDiss+cWind+cDiff+cP+cAb-cTend).*repmat(rAz(2:end,2:end).*inWag1(2:end,2:end),[1 1 nt])))),'linewidth',2)
legend('adv','diss','wind','diff','Press','timestep','tend','tot1','vortTend','tot2')
title('surface wag vorticity budget by area integral')

figure; plot(squeeze(nansum(nansum(cP.*repmat(rAz(2:end,2:end).*inWag1(2:end,2:end),[1 1 nt])))));
hold all
plot(squeeze(nansum(nansum(cPtropic.*repmat(rAz(2:end,2:end).*inWag1(2:end,2:end),[1 1 nt])))));
plot(squeeze(nansum(nansum(cPclinic.*repmat(rAz(2:end,2:end).*inWag1(2:end,2:end),[1 1 nt])))));
legend('pressure','barotropic part only','baroclinic part only')
title('vorticity budget pressure term')
%% circulation integral
% figure; [cWag,h]=contour(xcm,ycm,inWag1,[1 1]);
% xWag=cWag(1,2:494); yWag=cWag(2,2:494);
% vecWag=[xWag(:).';yWag(:).'];
% vecDS=diff(vecWag,1,2);
% for i=1:8
%     uWag(:,:,i)=griddata(xcm,ycm,Usnap(:,:,i),xWag,yWag);
%     vWag(:,:,i)=griddata(xcm,ycm,Vsnap(:,:,i),xWag,yWag);
%     tauEwag(:,:,i)=griddata(xcm,ycm,TauEa(:,:,i),xWag,yWag);
%     tauNwag(:,:,i)=griddata(xcm,ycm,TauNa(:,:,i),xWag,yWag);
%     hold1=(DYC(2:end,2:end).*(TauNa(2:end,2:end,i)-TauNa(1:end-1,2:end,i))-DXC(2:end,2:end).*(TauEa(2:end,2:end,i)-TauEa(2:end,1:end-1,i)))./rAz(2:end,2:end);
%     curlWind(:,:,i)=hold1;
% end
% vecTau=[reshape(tauEwag(1,1:end-1,:),1,492,8);reshape(tauNwag(1,1:end-1,:),1,492,8)];
% vecU=[reshape(uWag(1,1:end-1,:),1,492,8);reshape(vWag(1,1:end-1,:),1,492,8)];
% for i=1:8
% intUds(i)=sum(dot(vecU(:,:,i),vecDS));
% intTauDS(i)=sum(dot(vecTau(:,:,i),vecDS));
% intDelTauDA(i)=squeeze(nansum(nansum(curlWind(:,:,i).*rAz(2:end,2:end).*inWag1(2:end,2:end))));
% intVortDA(i)=squeeze(nansum(nansum(cVort(:,:,i).*rAz(2:end,2:end).*inWag1(2:end,2:end))));
% end

%% old
% %matlab curl requires meshgrid-like coordinates
% x=0:2000:max(max(xcm));
% y=0:2000:max(max(ycm));
% [xg,yg]=meshgrid(x,y);
% for i=1:nt
%     hold1=griddata(xcm,ycm,UCori(:,:,1,i),xg,yg);
%     UCoriM(:,:,i)=hold1;
%     hold1=griddata(xcm,ycm,VCori(:,:,1,i),xg,yg);
%     VCoriM(:,:,i)=hold1;
%     hold1=griddata(xcm,ycm,AdvU1(:,:,1,i),xg,yg);
%     AdvU1M(:,:,i)=hold1;
%     hold1=griddata(xcm,ycm,AdvV1(:,:,1,i),xg,yg);
%     AdvV1M(:,:,i)=hold1;
%     hold1=griddata(xcm,ycm,UDiss(:,:,1,i),xg,yg);
%     UDissM(:,:,i)=hold1;
%     hold1=griddata(xcm,ycm,VDiss(:,:,1,i),xg,yg);
%     VDissM(:,:,i)=hold1;
%     hold1=griddata(xcm,ycm,Uext(:,:,1,i),xg,yg);
%     UextM(:,:,i)=hold1;
%     hold1=griddata(xcm,ycm,Vext(:,:,1,i),xg,yg);
%     VextM(:,:,i)=hold1;
%     hold1=griddata(xcm,ycm,AbU(:,:,1,i),xg,yg);
%     AbUM(:,:,i)=hold1;
%     hold1=griddata(xcm,ycm,AbV(:,:,1,i),xg,yg);
%     AbVM(:,:,i)=hold1;
%     hold1=griddata(xcm,ycm,UDif2a(:,:,1,i),xg,yg);
%     UDifaM(:,:,i)=hold1;
%     hold1=griddata(xcm,ycm,VDif2a(:,:,1,i),xg,yg);
%     VDifaM(:,:,i)=hold1;
%     hold1=griddata(xcm,ycm,Utot1(:,:,1,i),xg,yg);
%     Utot1M(:,:,i)=hold1;
%     hold1=griddata(xcm,ycm,Vtot1(:,:,1,i),xg,yg);
%     Vtot1M(:,:,i)=hold1;
%     hold1=griddata(xcm,ycm,Utend(:,:,1,i)./86400,xg,yg);
%     UtendM(:,:,i)=hold1;
%     hold1=griddata(xcm,ycm,Vtend(:,:,1,i)./86400,xg,yg);
%     VtendM(:,:,i)=hold1;
%     hold1=griddata(xcm,ycm,UP(:,:,1,i),xg,yg);
%     UPM(:,:,i)=hold1;
%     hold1=griddata(xcm,ycm,VP(:,:,1,i),xg,yg);
%     VPM(:,:,i)=hold1;
% end
% 
% 
% %% surface layer vorticity budget
% disp('starting curl')
% for i=1:nt
% [betaVi,~]=curl(xg,yg,UCoriM(:,:,i),VCoriM(:,:,i));
% betaV(:,:,i)=betaVi;
% [advVi,~]=curl(xg,yg,AdvU1M(:,:,i),AdvV1M(:,:,i));
% advV(:,:,i)=advVi;
% [dissVi,~]=curl(xg,yg,UDissM(:,:,i),VDissM(:,:,i));
% dissV(:,:,i)=dissVi;
% [windVi,~]=curl(xg,yg,UextM(:,:,i),VextM(:,:,i));
% windV(:,:,i)=windVi;
% [pressVi,~]=curl(xg,yg,UPM(:,:,i),VPM(:,:,i));
% pressV(:,:,i)=pressVi;
% [timeVi,~]=curl(xg,yg,AbUM(:,:,i),AbVM(:,:,i));
% timeV(:,:,i)=timeVi;
% [diffVi,~]=curl(xg,yg,UDifaM(:,:,i),VDifaM(:,:,i));
% diffV(:,:,i)=diffVi;
% [totVi,~]=curl(xg,yg,Utot1M(:,:,i),Vtot1M(:,:,i));
% totV(:,:,i)=totVi;
% [tendVi,~]=curl(xg,yg,UtendM(:,:,i),VtendM(:,:,i));
% tendV(:,:,i)=tendVi;
% end
% 
% tot2V=diffV+timeV+pressV+windV+dissV+advV+betaV;
% clear AbV
% save('surfaceVorticity1.mat','*V','vort','inWag1','XC','x*','YC','y*','*Coast')
% 
% %% plot/check results of surface layer budget
% load('surfaceVorticity1.mat')
% figure; pcolor(xg,yg,sum(totV-tot2V,3)); shading 'flat'
% 
% inWag2=griddata(xcm,ycm,inWag1,xg,yg);
% figure; pcolor(xg,yg,inWag2); shading 'flat'; hold on
% contour(xcm,ycm,inWag1,[1 1])
% 
% inWag3=double(inWag2>0.8);
% figure; pcolor(xg,yg,inWag3); shading 'flat'; hold on
% contour(xcm,ycm,inWag1,[1 1],'k','linewidth',2)
% %%
% %inWag3=zeros(size(xg)); inWag3(150,250)=1;
% %nt=11;
% % bv=squeeze(nansum(nansum(betaV.*repmat(inWag3,[1 1 nt]))));
% % advect=squeeze(nansum(nansum(advV.*repmat(inWag3,[1 1 nt]))));
% % wind=squeeze(nansum(nansum(windV.*repmat(inWag3,[1 1 nt]))));
% % dissDiff=squeeze(nansum(nansum((dissV+diffV).*repmat(inWag3,[1 1 nt]))));
% % diss=squeeze(nansum(nansum((dissV).*repmat(inWag3,[1 1 nt]))));
% % diff=squeeze(nansum(nansum((diffV).*repmat(inWag3,[1 1 nt]))));
% % time=squeeze(nansum(nansum(timeV.*repmat(inWag3,[1 1 nt]))));
% % tend=squeeze(nansum(nansum(tendV.*repmat(inWag3,[1 1 nt])./86400)));
% % tot=squeeze(nansum(nansum(totV.*repmat(inWag3,[1 1 nt]))));
% 
% bv=squeeze(betaV(150,250,:));
% advect=squeeze(advV(150,250,:));
% wind=squeeze(windV(150,250,:));
% dissDiff=squeeze(dissV(150,250,:)+diffV(150,250,:));
% diss=squeeze(dissV(150,250,:));
% diff=squeeze(diffV(150,250,:));
% time=squeeze(timeV(150,250,:));
% tend=squeeze(tendV(150,250,:))./86400;
% tot=squeeze(totV(150,250,:));
% 
% figure; plot(advect); hold all; plot(bv); plot(wind);
% plot(dissDiff); plot(time);  plot(tot,'k--'); 
% plot(-advect+bv+wind+dissDiff+time,'linewidth',2)
% plot(tend);
% %plot(advect+bv); plot(diss+diff+wind)
% legend('advect','bv','wind','diss+diff','ab','total','sum','tend')%'advect+bv','dissdiff+wind')
% %% circulation integral at surface
% %can't do this with existing data- wasn't saving out velocity snapshots
% 
% %circulation=sum(dot([u v],[dx dy]));
% %wind=sum(dot([taux tauy],[dx dy]));
% %disf=nansum(nansum((dissV+diffV).*repmat(inWag3,[1 1 nt])));
% %no pressure gradient in closed loop at surface
% %advect=dot(u,grad(vort)); %??
% 
