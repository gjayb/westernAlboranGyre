%how isopycnals move due to vertical diffusion of T,S
%no horizontal because it is orders of magnitude smaller
%this is to get an approximation of the cross-isopycnal flux with the
%assumption that advection is along isopycnals

load('geometrySpinupSteady','dInterface','d','nLayers','XC','YC')
load('distancesAreas','hFacC','RAC')
dZ(1,1,1:46)=diff(dInterface);
cellVol=repmat(RAC,[1 1 46]).*hFacC.*repmat(dZ,[700 200 1]);
cellVol4=repmat(cellVol,[1 1 1 148]);

% load('fluxesDifS148day.mat','DifEs','DifIs')
% load('fluxesOtherS148day.mat','KpS')
% difSz2=(DifEs(:,:,2:46,:)+DifIs(:,:,2:46,:)+KpS(:,:,2:46,:)-DifEs(:,:,1:45,:)-DifIs(:,:,1:45,:)-KpS(:,:,1:45,:))./cellVol4(:,:,1:45,:);
% clear DifEs DifIs KpS
% 
% load('tsPress162NF.mat','PractSave')
% PractSave=PractSave(:,:,:,1:148);
% S2=PractSave(:,:,1:45,:)+86400.*difSz2;
% clear PractSave
% disp('S2 done')
% 
% load('fluxesDifT148day.mat','DifzEt','DifIt')
% load('fluxesOtherT148day.mat','KpT')
% difTz2=(DifzEt(:,:,2:46,:)+DifIt(:,:,2:46,:)+KpT(:,:,2:46,:)-DifzEt(:,:,1:45,:)-DifIt(:,:,1:45,:)-KpT(:,:,1:45,:))./cellVol4(:,:,1:45,:);
% clear DifEt DifIt KpT
% 
% load('tsPress162NF.mat','PotTave')
% PotTave=PotTave(:,:,:,1:148);%fix by adding this semicolon if rerunning on server!
% T2=PotTave(:,:,1:45,:)+86400.*difTz2;
% disp('T2 done')
% %%
% load('tsPress162NF.mat','PressDbar')
% Sigma=zeros([700 200 45 148]);
% for k=1:148
%     k
%     for di=1:45
%         SA=gsw_SA_from_SP(S2(:,:,di,k),PressDbar(:,:,di,k),XC,YC);
%         CT=gsw_CT_from_pt(SA,T2(:,:,di,k));
%         Sigma(:,:,di,k)=gsw_sigma0(SA,CT);
%         
%     end
% end
% disp('Sigma done, saving')
% save('diffusiveTSsigma.mat','S2','T2','Sigma','-v7.3')
%%
load('diffusiveTSsigma.mat')
dBin=0.5.*(dInterface(2:end)+dInterface(1:end-1));
isopycs=[26.3 26.5 26.75 27 27.5 28];
isoFSs=[263 265 2675 27 275 28];
nLayers(nLayers==46)=45;
Sigma(isnan(Sigma))=40;%all land points, making sigma large so isopycnals defined in water
for isoi=1:6
    isopyc=isopycs(isoi)
    isoFS=isoFSs(isoi);
 for i=1:148
    i
        for j=1:700
            for k=1:200
                if d(j,k)>0
                depths1=dBin(1:nLayers(j,k));
                if nLayers(j,k)>1
                isoDepth(j,k,i)=interp1(squeeze(Sigma(j,k,1:nLayers(j,k),i)),depths1,isopyc);
                end
                else
                    isoDepth(j,k,i)=0;
                end
            end
        end
    
 end
    fn=strcat('diffusiveMovedIsopycnal',num2str(isoFS),'.mat')
    save(fn,'isoDepth')
end
disp('done')
%%
isoFS=275;
fn=strcat('diffusiveMovedIsopycnal',num2str(isoFS),'.mat')
load(fn)
isoNew=isoDepth;
load('iso275depthNFrev')
nt=148;
 snapDiff=isoNew-isoDepth(:,:,1:nt);%positive difference is deepening
 volDiff275in27=inwag27areas.*(isoNew(:,:,1:146)-isoDepth(:,:,1:146));%positive is deepening, so volume change is into gyre
 dvdtDiff275in27=squeeze(nansum(nansum(volDiff275in27)))/86400;
 dvdtDiff275in275=squeeze(nansum(nansum(inwag275areas.*(isoNew(:,:,1:146)-isoDepth(:,:,1:146)))))./86400;%positive is deepening, so volume change is into gyre

%% diffusive motion of isopycnals at 'steps'
load('diffusiveMovedIsopycnal263.mat')
isoNew=isoDepth;
load('iso263depthNFrev')
dvdt263StepTd=squeeze(nansum(nansum( (inwagSareas1-inwag263areas).*(isoNew(:,:,1:146)-isoDepth(:,:,1:146)) )))./86400;
dvdt263StepBu=squeeze(nansum(nansum( (inwag263areas-inwag265areas).*(isoNew(:,:,1:146)-isoDepth(:,:,1:146)) )))./86400;

load('diffusiveMovedIsopycnal265.mat')
isoNew=isoDepth;
load('iso265depthNFrev')
dvdt265StepTd=squeeze(nansum(nansum( (inwag263areas-inwag265areas).*(isoNew(:,:,1:146)-isoDepth(:,:,1:146)) )))./86400;
dvdt265StepBu=squeeze(nansum(nansum( (inwag265areas-inwag2675areas).*(isoNew(:,:,1:146)-isoDepth(:,:,1:146)) )))./86400;

load('diffusiveMovedIsopycnal2675.mat')
isoNew=isoDepth;
load('iso2675depthNFrev')
dvdt2675StepTd=squeeze(nansum(nansum( (inwag265areas-inwag2675areas).*(isoNew(:,:,1:146)-isoDepth(:,:,1:146)) )))./86400;
dvdt2675StepBu=squeeze(nansum(nansum( (inwag2675areas-inwag27areas).*(isoNew(:,:,1:146)-isoDepth(:,:,1:146)) )))./86400;

load('diffusiveMovedIsopycnal27.mat')
isoNew=isoDepth;
load('iso27depthNFrev')
dvdt27StepTd=squeeze(nansum(nansum( (inwag2675areas-inwag27areas).*(isoNew(:,:,1:146)-isoDepth(:,:,1:146)) )))./86400;
dvdt27StepBu=squeeze(nansum(nansum( (inwag27areas-inwag275areas).*(isoNew(:,:,1:146)-isoDepth(:,:,1:146)) )))./86400;
%% SSH component of volume, dV/dt
% load('uvwSSHDailyDepth1rotated148F.mat', 'SSHa')
% load('sshSnap.mat', 'sshSnap')
volssh1=squeeze(nansum(nansum(SSHa(:,:,1:146).*inwagSareas1)));
volssh1s=squeeze(nansum(nansum(sshSnap(:,:,1:146).*inwagSareas1)));
volssh2=squeeze(nansum(nansum(SSHa(:,:,1:146).*inwag263areas)));
volssh2s=squeeze(nansum(nansum(sshSnap(:,:,1:146).*inwag263areas)));
%% 3 day running mean gate flux and dVdt (incl ssh)
dVdt11=dVdt1+diff(volssh1).'/86400; dVdt11(13)=0;
dVdt3day1=[0 dVdt11(1:end-2)+dVdt11(2:end-1)+dVdt11(3:end) 0]./3;
dVdt1s1=dVdt1s+diff(volssh1s).'/86400; dVdt1s1(13)=0;
dVdt3day1s=[0 dVdt1s1(1:end-2)+dVdt1s1(2:end-1)+dVdt1s1(3:end) 0]./3;

dVdt21=dVdt2+diff(volssh2).'/86400; dVdt12(13)=0;
dVdt3day2=[0 dVdt12(1:end-2)+dVdt12(2:end-1)+dVdt12(3:end) 0]./3;
dVdt2s1=dVdt2s+diff(volssh2s).'/86400; dVdt2s1(13)=0;
dVdt3day2s=[0 dVdt2s1(1:end-2)+dVdt2s1(2:end-1)+dVdt2s1(3:end) 0]./3;

gateFluxTot3day=[0 gateFluxTot(1:end-2)+gateFluxTot(2:end-1)+gateFluxTot(3:end) 0]./3;
gateFluxTot3dayi=[0 gateFluxToti(1:end-2)+gateFluxToti(2:end-1)+gateFluxToti(3:end) 0]./3;
gateFluxTot3day2=[0 gateFluxTot2(1:end-2)+gateFluxTot2(2:end-1)+gateFluxTot2(3:end) 0]./3;
gateFluxTot3day2i=[0 gateFluxTot2i(1:end-2)+gateFluxTot2i(2:end-1)+gateFluxTot2i(3:end) 0]./3;
%% plot gate-dVdt, small fluxes
%%
 load('sigma275sigmagradientNFrev.mat')
 [nxg,nyg,ntg]=size(gsxI)
 gsxI(700,200,ntg)=0; gsyI(700,200,ntg)=0; gszI(700,200,ntg)=0; dZi(700,200,ntg)=0;
mag1=sqrt(gsxI(:,:,1:nt).^2+gsyI(:,:,1:nt).^2+gszI(:,:,1:nt).^2);%magnitude of grad(sigma)
vecGI2=[reshape(gsxI(:,:,1:nt)./mag1(:,:,1:nt),1,[]);reshape(gsyI(:,:,1:nt)./mag1(:,:,1:nt),1,[]);reshape(gszI(:,:,1:nt)./mag1(:,:,1:nt),1,[])];
clear gsxI gsyI

 load('distancesAreas','RAC','dxg','dyg')
DXG=reshape(dxg,[700 200]);
DYG=reshape(dyg,[700 200]);
 surfaceXN=dZi(:,:,1:148).*repmat(DYG,[1 1 148]);
 surfaceYN=dZi(:,:,1:148).*repmat(DXG,[1 1 148]);
 load('geometrySpinupSteady','Angle*','*C*')
 %%
 surfaceX=surfaceXN.*repmat(AngleCS,[1 1 148]) - surfaceYN.*repmat(AngleSN,[1 1 148]);  
 surfaceY=surfaceXN.*repmat(AngleSN,[1 1 148]) + surfaceYN.*repmat(AngleCS,[1 1 148]); 
sizesY=size(surfaceY) 
 surfaceZ=repmat(RAC,[1 1 nt]);

vecSurf=[reshape(surfaceX,1,[]);reshape(surfaceY,1,[]);reshape(surfaceZ,1,[])];%all positive

 snapVel=snapDiff.*gszI(:,:,1:148)./(mag1(:,:,1:148).*86400);%positive down
surfaces=reshape(dot(vecSurf,vecGI2),[700 200 148]);%vecGI usually negative

  vertVol=snapVel.*surfaces;%this is positive down*negative surface areas, giving Positive UP
