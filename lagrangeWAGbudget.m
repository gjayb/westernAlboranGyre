% %% explanation
% %diffusion for lagrangian WAG based on Eulerian diagnostic fluxes
% 
% %% location
% %load('wagAreaAndFluxGate.mat', 'inWAG')%adjust by adding 26, 265, or 27 to filename
% load('geometrySpinupSteady.mat','XC','YC')
% load('wagAreaAndFluxGate27.mat', '*Closed')
% figure;
% for i=1:131
%     c=[];
%     if latWagClosed(i,1)>0
%         np=find(latWagClosed(i,:)>0,1,'last');
%         inWAG=inpolygon(XC,YC,lonWagClosed(i,1:np),latWagClosed(i,1:np));
%         [c,h]=contour(XC,YC,double(inWAG),[1 1]);
%     end
%     if ~isempty(c)
%     npointsi=find(c(2,:)==floor(c(2,:)));
%     varhold=max(c(2,npointsi));
%     i1=1+find(c(2,:)==varhold);
%     iend=varhold+i1-1;
%     inWAG27(:,:,i)=inpolygon(XC,YC,c(1,i1:iend),c(2,i1:iend));
%     end
% end
% clear inWAG
% 
% load('wagAreaAndFluxGate265.mat', '*Closed')
% figure;
% for i=1:131
%     c=[];
%     if latWagClosed(i,1)>0
%         np=find(latWagClosed(i,:)>0,1,'last');
%         inWAG=inpolygon(XC,YC,lonWagClosed(i,1:np),latWagClosed(i,1:np));
%         [c,h]=contour(XC,YC,double(inWAG),[1 1]);
%     end
%     if ~isempty(c)
%     npointsi=find(c(2,:)==floor(c(2,:)));
%     varhold=max(c(2,npointsi));
%     i1=1+find(c(2,:)==varhold);
%     iend=varhold+i1-1;
%     inWAG265(:,:,i)=inpolygon(XC,YC,c(1,i1:iend),c(2,i1:iend));
%     end
% end
% clear inWAG
% 
% load('wagAreaAndFluxGate26.mat', '*Closed')
% figure;
% for i=1:138
%     c=[];
%     if latWagClosed(i,1)>0
%         np=find(latWagClosed(i,:)>0,1,'last');
%         inWAG=inpolygon(XC,YC,lonWagClosed(i,1:np),latWagClosed(i,1:np));
%         [c,h]=contour(XC,YC,double(inWAG),[1 1]);
%     end
%     if ~isempty(c)
%     npointsi=find(c(2,:)==floor(c(2,:)));
%     varhold=max(c(2,npointsi));
%     i1=1+find(c(2,:)==varhold);
%     iend=varhold+i1-1;
%     inWAG26(:,:,i)=inpolygon(XC,YC,c(1,i1:iend),c(2,i1:iend));
%     end
% end
% clear inWAG
% 
% load('wagAreaAndFluxGate.mat', '*Closed')
% figure;
% for i=9:141
%     c=[];
%     if latWagClosed(i-8,1)>0
%         np=find(latWagClosed(i-8,:)>0,1,'last');
%         inWAG=inpolygon(XC,YC,lonWagClosed(i-8,1:np),latWagClosed(i-8,1:np));
%         [c,h]=contour(XC,YC,double(inWAG),[1 1]);
%     end
%     if ~isempty(c)
%     npointsi=find(c(2,:)==floor(c(2,:)));
%     varhold=max(c(2,npointsi));
%     i1=1+find(c(2,:)==varhold);
%     iend=varhold+i1-1;
%     inWAGS(:,:,i)=inpolygon(XC,YC,c(1,i1:iend),c(2,i1:iend));
%     end
% end
% clear inWAG
% 
% inWag3=zeros([700 200 16 131]);
% load('varying148ts16levelsRho.mat','Sigma')
% %%
% Sigma=Sigma(:,:,1:16,1:131);
% for k=1:16
% inWag3(:,:,k,:)=double(squeeze(Sigma(:,:,k,:))<26).*inWAGS(:,:,1:131)...
%     +double(squeeze(Sigma(:,:,k,:))<26.5).*double(squeeze(Sigma(:,:,k,:))>=26).*inWAG26(:,:,1:131)...
%     +double(squeeze(Sigma(:,:,k,:))<27).*double(squeeze(Sigma(:,:,k,:))>=26.5).*inWAG265...
%     +double(squeeze(Sigma(:,:,k,:))<27.5).*double(squeeze(Sigma(:,:,k,:))>=27).*inWAG27;
% end
% 
% save('lagrangeWAGboundary.mat','inW*','-v7.3')
%% potential temperature and salinity contents
load('lagrangeWAGboundary.mat','inW*')
%load('varying148ts16levelsRho.mat','T');
load('fluxesT58day.mat','Ttend')

for i=2:58
    T1=T(:,:,:,i);
    Ttend1=Ttend(:,:,:,i);
    
    stayIn=(inWag3(:,:,:,i)==1)&(inWag3(:,:,:,i-1)==1);
    out1=(inWag3(:,:,:,i)==0)&(inWag3(:,:,:,i-1)==1);
    in1=(inWag3(:,:,:,i)==1)&(inWag3(:,:,:,i-1)==0);
    inNow=inWag3(:,:,:,i)==1;
    
    dTdt(i)=(sum(Ttend1(stayIn))+sum(T1(in1))-sum(T1(out1)))/86400;
    dTdt2(i)=squeeze(nansum(nansum(nansum(Ttend1(inNow)))))/86400;
end
clear Ttend T

load('varying148ts16levelsRho.mat','S');
load('fluxesS58day.mat','Stend')

for i=2:58
    S1=S(:,:,:,i);
    Stend1=Stend(:,:,:,i);
    
    stayIn=(inWag3(:,:,:,i)==1)&(inWag3(:,:,:,i-1)==1);
    out1=(inWag3(:,:,:,i)==0)&(inWag3(:,:,:,i-1)==1);
    in1=(inWag3(:,:,:,i)==1)&(inWag3(:,:,:,i-1)==0);
    inNow=inWag3(:,:,:,i)==1;
    
    dSdt(i)=(sum(Stend1(stayIn))+sum(S1(in1))-sum(S1(out1)))/86400;
    dSdt2(i)=squeeze(nansum(nansum(nansum(Stend1(inNow)))))/86400;
end
clear S Stend

save('lagrangeWAGcontentsTS58.mat','dTdt*','dSdt*')
%% T fluxes the Eulerian way, from diagnostics- load 
load('distancesAreas.mat','RAC','hFacC','dxc','dyc')
DXC=reshape(dxc,[700 200]); DYC=reshape(dyc,[700 200]);
load('geometrySpinupSteady.mat','dInterface')
nt=58;

dZ(1,1,1:46)=diff(dInterface);
cellVol=repmat(RAC,[1 1 46]).*hFacC.*repmat(dZ,[700 200 1]);
cellVol4=repmat(cellVol,[1 1 1 nt]);

load('fluxesT58day.mat')
load('lagrangeWAGboundary.mat','inWag3')
hold1=zeros([700 200 46 nt]);
[a1,a2,a3,a4]=size(inWag3);
hold1(1:a1,1:a2,1:a3,1:nt)=inWag3(:,:,:,1:nt);
inWag3=hold1; clear hold1
%% T fluxes the Eulerian way, from diagnostics
advTh=squeeze(nansum(nansum(nansum((AdvUt(1:end-1,1:end-1,:,:)-AdvUt(2:end,1:end-1,:,:)+AdvVt(1:end-1,1:end-1,:,:)...
    -AdvVt(1:end-1,2:end,:,:)).*inWag3(1:end-1,1:end-1,:,1:nt)./cellVol4(1:end-1,1:end-1,:,:)))));
advTz=squeeze(nansum(nansum(nansum(...
    (AdvWt(:,:,2:46,:)-AdvWt(:,:,1:45,:)).*inWag3(:,:,1:45,1:nt)./cellVol4(:,:,1:45,:)))));
difTz=squeeze(nansum(nansum(nansum((DifzEt(:,:,2:46,:)+DifIt(:,:,2:46,:)+KpT(:,:,2:46,:)-DifzEt(:,:,1:45,:)-DifIt(:,:,1:45,:)-KpT(:,:,1:45,:)).*inWag3(:,:,1:45,1:nt)./cellVol4(:,:,1:45,:)))));
difTh=squeeze(nansum(nansum(nansum((DifxEt(1:end-1,1:end-1,:,:)-DifxEt(2:end,1:end-1,:,:)+DifyEt(1:end-1,1:end-1,:,:)...
    -DifyEt(1:end-1,2:end,:,:)).*inWag3(1:end-1,1:end-1,:,1:nt)./cellVol4(1:end-1,1:end-1,:,:)))));
for i=1:nt
surfTtend3(i)=nansum(nansum(squeeze(inWag3(:,:,1,i)).*(-squeeze(WT(:,:,1,i)))./(dZ(1,1,1).*hFacC(:,:,1))));
end

load('varying148ts16levelsRho.mat', 'Rho')
Rho=squeeze(Rho(:,:,1,1:nt));
load('cp16levels148.mat', 'cp')
cp=squeeze(cp(:,:,1,1:nt));
surfQ=squeeze(nansum(nansum(Tflux.*inWag3(:,:,1,1:nt)./(dZ(1,1,1).*repmat(hFacC(:,:,1),[1 1 nt]).*Rho.*cp))));
clear Rho cp Tflux

load('varying148ts16levelsRho.mat', 'T');
T=T(:,:,1:15,1:nt);
load('dwadz148levels15.mat','dWdZ')
dWdZ=dWdZ(:,:,:,1:nt);
tdwdz=squeeze(nansum(nansum(nansum(T.*dWdZ.*inWag3(:,:,1:15,1:nt)))));
clear T dWdZ cellVol DXC DYC hFacC RAC

load('lagrangeWAGcontentsTS58.mat','dTdt*')
save('lagrangeWagTempBudget58.mat')
%% S fluxes, Eulerian way, load
clear
load('distancesAreas.mat','RAC','hFacC','dxc','dyc')
DXC=reshape(dxc,[700 200]); DYC=reshape(dyc,[700 200]);
load('geometrySpinupSteady.mat','dInterface')
nt=58;

dZ(1,1,1:46)=diff(dInterface);
cellVol=repmat(RAC,[1 1 46]).*hFacC.*repmat(dZ,[700 200 1]);
cellVol4=repmat(cellVol,[1 1 1 nt]);

load('fluxesS58day.mat')
load('lagrangeWAGboundary.mat','inWag3')
hold1=zeros([700 200 46 nt]);
[a1,a2,a3,a4]=size(inWag3);
hold1(1:a1,1:a2,1:a3,1:nt)=inWag3(:,:,:,1:nt);
inWag3=hold1; clear hold1

%% S fluxes from diagnostics
AdvSh=squeeze(nansum(nansum(nansum(...
    (AdvUs(1:end-1,1:end-1,:,:)-AdvUs(2:end,1:end-1,:,:)+AdvVs(1:end-1,1:end-1,:,:)...
    -AdvVs(1:end-1,2:end,:,:)).*inWag3(1:end-1,1:end-1,:,1:nt)./cellVol4(1:end-1,1:end-1,:,:)...
        ))));
AdvSz=squeeze(nansum(nansum(nansum(...
    (AdvWs(:,:,2:46,:)-AdvWs(:,:,1:45,:)).*inWag3(:,:,1:45,1:nt)./cellVol4(:,:,1:45,:)...
    ))));
difSz=squeeze(nansum(nansum(nansum((DifEs(:,:,2:46,:)+DifIs(:,:,2:46,:)+KpS(:,:,2:46,:)-DifEs(:,:,1:45,:)-DifIs(:,:,1:45,:)-KpS(:,:,1:45,:)).*inWag3(:,:,1:45,1:nt)./cellVol4(:,:,1:45,:)))));
difSh=squeeze(nansum(nansum(nansum((DifxEs(1:end-1,1:end-1,:,:)-DifxEs(2:end,1:end-1,:,:)+DifyEs(1:end-1,1:end-1,:,:)...
    -DifyEs(1:end-1,2:end,:,:)).*inWag3(1:end-1,1:end-1,:,1:nt)./cellVol4(1:end-1,1:end-1,:,:)))));
for i=1:nt
surfStend3(i)=nansum(nansum(inWag3(:,:,1,i).*(-squeeze(WS(:,:,1,i)))./(dZ(1,1,1).*hFacC(:,:,1))));
end
load('varying148ts16levelsRho.mat', 'Rho')
Rho=squeeze(Rho(:,:,1,1:nt));
surfF=squeeze(nansum(nansum(Sflux.*squeeze(inWag3(:,:,1,1:nt))./(dZ(1,1,1).*repmat(hFacC(:,:,1),[1 1 nt]).*Rho))));

load('varying148ts16levelsRho.mat', 'S');
S=S(:,:,1:15,1:nt);
load('dwadz148levels15.mat','dWdZ')
dWdZ=dWdZ(:,:,:,1:nt);
sdwdz=squeeze(nansum(nansum(nansum(S.*dWdZ.*inWag3(:,:,1:15,1:nt)))));
clear S dWdZ Sflux cellVol DXC DYC hFacC RAC

load('lagrangeWAGcontentsTS58.mat','dSdt*')
save('lagrangeWagSalBudget58.mat')
%% plots budget
load('lagrangeWAGcontentsTS58.mat','dSdt*')
load('lagrangeWagSalBudget58.mat','AdvS*','sdwdz','difS*','surfStend*','surfF')
figure; plot(-dSdt); hold all
plot(AdvSh+sdwdz); plot(AdvSz-sdwdz);
plot(difSh); plot(difSz); plot(surfStend3'+surfF)
plot(AdvSh(1:45)+AdvSz(1:45)+difSh(1:45)+difSz(1:45)+surfStend3(1:45)'+surfF(1:45))
legend('-dS/dt','AdvH','AdvV','difH','difV','surf','tot fluxes')

figure; plot(-dSdt2); hold all
plot(AdvSh+sdwdz); plot(AdvSz-sdwdz);
plot(difSh); plot(difSz); plot(surfStend3'+surfF)
plot(AdvSh(1:45)+AdvSz(1:45)+difSh(1:45)+difSz(1:45)+surfStend3(1:45)'+surfF(1:45))
legend('-dS/dt','AdvH','AdvV','difH','difV','surf','tot fluxes')

%% gate advective fluxes from lagrange calculations before, compared to the advective ones from diagnostics




