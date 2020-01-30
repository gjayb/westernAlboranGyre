%gathers snapshots of T, S, SSH @day-ends
%saves, gets daily-mean fluxes, saves
%calculates density, specific heat
%calculates Euler WAG contents

% times=8640:8640:(14*8640);
% 
% S1=rdmds('S',times);
% T1=rdmds('T',times);
% SSH1=rdmds('SSH',times);
% save('stssh14daysnap.mat','-v7.3')
% 
% clear S1 T1 SSH1
%%
addpath('/nobackup1/gbrett/mStuff/')
%addpath('/nobackup1/gbrett/addforce/')
times=8640:8640:(162*8640);
days=1:59;
times2=[1 8640:8640:(162*8640)];
% 
% AdvUs=rdmds('AdvUs',times);
% AdvVs=rdmds('AdvVs',times);
% AdvWs=rdmds('AdvWs',times);
% save('fluxesAdvS59day.mat','-v7.3');
% clear Adv*
% 
% DifEs=rdmds('DifZeS',times);
% DifIs=rdmds('DifZiS',times);
% DifxEs=rdmds('DifXeS',times);
% DifyEs=rdmds('DifYeS',times);
% disp('half of S read')
% save('fluxesDifS59day.mat','-v7.3')
% clear Dif*
% WS=rdmds('WSM',times);
% US=rdmds('USM',times);
% VS=rdmds('VSM',times);
% save('fluxesVelS59day.mat','-v7.3')
% clear *S
% Stend=rdmds('Stend',times);
% Sflux=rdmds('Sflux',times);
% KpS=rdmds('KppGS',times);
% Fsw=rdmds('Foce',times);
% save('fluxesOtherS59day.mat','-v7.3')
% % % 
%  clear S* *S Fsw
%%
% clear
% times=8640:8640:(148*8640)
% days=1:148;
% nt=length(days);
% 
% AdvUt=rdmds('AdvUt',times);
% AdvVt=rdmds('AdvVt',times);
% AdvWt=rdmds('AdvWt',times);
% save('fluxesAdvT59day.mat','-v7.3');
% clear Adv*
% DifzEt=rdmds('DifZeT',times);
% DifIt=rdmds('DifZiT',times);
% DifxEt=rdmds('DifXeT',times);
% DifyEt=rdmds('DifYeT',times);
% disp('half of T read')
% save('fluxesDifT59day.mat','-v7.3')
% clear Dif*
% WT=rdmds('WTM',times);
% UT=rdmds('UTM',times);%missing b/c misnumbered data.diagnostics
% VT=rdmds('VTM',times);
% save('fluxesVelT59day.mat','-v7.3')
% clear *T
% Ttend=rdmds('Ttend',times);
% Tflux=rdmds('Tflux',times);
% KpT=rdmds('KppGT',times);
% %Qsw=rdmds('Qoce',times);
% Qnet=rdmds('QNETave',times);
% save('fluxesOtherT59day.mat','-v7.3')
% clear Q* K* Ttend Tflux
% 
% SSH=rdmds('SSHave',times);
% %Ta=rdmds('Tave',times);
% RhoA=rdmds('RhoAAve',times);
% Rho0=rdmds('RhoRef');
% Rho=RhoA+repmat(Rho0,[700 200 1 nt]); clear RhoA
% save('rhoSSH148day.mat','-v7.3')
% 
% clear A* D*
% Ssnap=rdmds('S',times2);
% Tsnap=rdmds('T',times2);
% Stendsnap=diff(Ssnap,1,4); clear Ssnap
% Ttendsnap=diff(Tsnap,1,4); clear Tsnap
% load('fluxesOtherS59day.mat','Stend')
% restartS=Stendsnap-Stend; clear Stend
% load('fluxesOtherT59day.mat','Ttend')
% restartT=Ttendsnap-Ttend; clear Ttend
% save('tsTendsnap.mat','-v7.3'); %clear Stendsnap Ttendsnap
%%
% load('edgesWAGeuler2017NF.mat','inWag')
% load('tsTendsnap.mat')
% load('fluxesOtherT59day.mat','Ttend')
% load('fluxesOtherS59day.mat','Stend')
% tendT=squeeze(nansum(nansum(nansum(Ttend.*repmat(inWag,[1 1 1 59])))));
% tendS=squeeze(nansum(nansum(nansum(Stend.*repmat(inWag,[1 1 1 59])))));
% resT=squeeze(nansum(nansum(nansum(restartT.*repmat(inWag,[1 1 1 59])))));
% resS=squeeze(nansum(nansum(nansum(restartS.*repmat(inWag,[1 1 1 59])))));
% clear *end restart*
% save('tendRestart.mat')
%%
% clear
% load('stssh14daysnap.mat','S1','T1')
% load('geometrySpinupSteady.mat')
% 
% disp('load done')
% 
% dBin=0.5*dInterface(1:end-1)+0.5*dInterface(2:end);
% p(1,1,1:46)=gsw_p_from_z(-dBin,36);
% p=repmat(p,[700 200 1 14]);
% disp('p made')
% CT=gsw_CT_from_pt(S1,T1);
% disp('1/4')
% Rho=gsw_rho(S1,CT,p);
% disp('2/4')
% Sigma=gsw_sigma0(S1,CT);
% disp('3/4')
% CP=gsw_cp_t_exact(S1,T1,p);
% 
% disp('saving')
% save('snapshots14days.mat','-v7.3')


%% quantities in eulerian wag over time

% load('edgesWAGeuler2017NR.mat','dInterface','inWag')
% inWag=inWag(:,:,1:16);
% load('distancesAreas.mat','RAC','hFacC')
% load('ssh148.mat','ssh148');
% RAC3=repmat(RAC,[1 1 16]);
% 
% clear dZ3
% dZ3(1,1,1:16)=dInterface(2:17)-dInterface(1:16);
% dZ3=repmat(dZ3,[700 200 1]);
% 
% load('contentsWAGeuler2017.mat','volWagEuler1')
% load('stssh14daysnap.mat','SSH1')
% volWagEuler=volWagEuler1*ones(14,1)+squeeze(sum(sum(SSH1.*repmat(inWag(:,:,1).*RAC,[1 1 14]))));
% 
% load('snapshots14days.mat','Rho');
% 
% massWagEuler=squeeze(sum(sum(sum(repmat(RAC3.*inWag.*dZ3,[1 1 1 14]).*Rho(:,:,1:16,:)))))+squeeze(sum(sum(SSH1.*squeeze(Rho(:,:,1,:)).*repmat(inWag(:,:,1).*RAC,[1 1 14]))));
% 
% load('snapshots14days.mat','S1');
% S1(S1<10)=NaN;
% saltWagEuler=squeeze(nansum(nansum(nansum(repmat(RAC3.*inWag.*dZ3,[1 1 1 14]).*Rho(:,:,1:16,:).*S1(:,:,1:16,:)))))+squeeze(nansum(nansum(SSH1.*squeeze(Rho(:,:,1,:).*S1(:,:,1,:)).*repmat(inWag(:,:,1).*RAC,[1 1 14]))));
% 
% 
% clear S S2
% load('snapshots14days.mat', 'CP')
% load('snapshots14days.mat','T1');
% T1(T1<1)=NaN;
% T=T1+273.15;%kelvin!
% heatWagEuler=squeeze(nansum(nansum(nansum(repmat(RAC3.*inWag.*dZ3,[1 1 1 14]).*Rho(:,:,1:16,:).*CP(:,:,1:16,:).*T(:,:,1:16,:)))))+squeeze(nansum(nansum(SSH1.*squeeze(Rho(:,:,1,:).*CP(:,:,1,:).*T(:,:,1,:)).*repmat(inWag(:,:,1).*RAC,[1 1 14]))));
% 
% clear T Rho cp *3 
% 
% save('contentsWAGeuler14days.mat','*Euler*')
%% trial heat budget, Euler WAG
% 
% load('fluxesT11day.mat')
% load('edgesWAGeuler2017.mat','inWag','open*','dInterface','XC','YC','*Coast')
% load('distancesAreas.mat','RAC','hFacC')
% 
% dZ(1,1,1:46)=diff(dInterface);
% cellVol=repmat(RAC,[1 1 46]).*hFacC.*repmat(dZ,[700 200 1]);
% cellVol4=repmat(cellVol,[1 1 1 11]);
% 
% advTh=squeeze(nansum(nansum(nansum((AdvUt(1:end-1,1:end-1,:,:).*repmat(openW(1:end-1,1:end-1,:,:),[1 1 1 11])...
%     -AdvUt(2:end,1:end-1,:,:).*repmat(openE(1:end-1,1:end-1,:,:),[1 1 1 11])+AdvVt(1:end-1,1:end-1,:,:).*repmat(openS(1:end-1,1:end-1,:,:),[1 1 1 11])...
%     -AdvVt(1:end-1,2:end,:,:).*repmat(openN(1:end-1,1:end-1,:,:),[1 1 1 11]))./cellVol4(1:end-1,1:end-1,:,:)))));
% advTz=squeeze(nansum(nansum(nansum(AdvWt(:,:,2:46,:).*repmat(openZd(:,:,1:45,:),[1 1 1 11])./cellVol4(:,:,1:45,:)))));
% 
% difTz=squeeze(nansum(nansum(nansum((DifzEt(:,:,2:46,:)+DifIt(:,:,2:46,:)+KpT(:,:,2:46,:)).*repmat(openZd(:,:,1:45,:),[1 1 1 11])./cellVol4(:,:,1:45,:)))));
% difTh=squeeze(nansum(nansum(nansum((DifxEt(1:end-1,1:end-1,:,:).*repmat(openW(1:end-1,1:end-1,:,:),[1 1 1 11])...
%     -DifxEt(2:end,1:end-1,:,:).*repmat(openE(1:end-1,1:end-1,:,:),[1 1 1 11])+DifyEt(1:end-1,1:end-1,:,:).*repmat(openS(1:end-1,1:end-1,:,:),[1 1 1 11])...
%     -DifyEt(1:end-1,2:end,:,:).*repmat(openN(1:end-1,1:end-1,:,:),[1 1 1 11]))./cellVol4(1:end-1,1:end-1,:,:)))));
% 
% tsurfcor=squeeze(nansum(nansum(squeeze(WT(:,:,1,:)).*repmat(RAC,[1 1 11]))))./sum(sum(RAC));
% %short1(1,1,1:46)=0.62*exp(dInterface(1:46)/0.6)+0.38*exp(dInterface(1:46)/20);%no shortwave!
% %short2(1,1,1:46)=0.62*exp(dInterface(2:47)/0.6)+0.38*exp(dInterface(2:47)/20);
% %shortwaveApprox=Qsw/(1026*4000)/(repmat(dZ,[700 200 1 11]).*repmat(hFacC,[1 1 1 11]))*repmat(short1-short2,[700 200 1 11]);
% for i=1:11
% surfTtend(i)=nansum(nansum(inWag(:,:,1).*(tsurfcor(i).*ones([700 200])-squeeze(WT(:,:,1,i)))./(dZ(1,1,1).*hFacC(:,:,1))));
% end
% surfQapprox=squeeze(nansum(nansum(Tflux.*repmat(inWag(:,:,1),[1 1 11])./(1026*4000*dZ(1,1,1).*repmat(hFacC(:,:,1),[1 1 11])))));%approx bc not sure what rhoConst or Cp should be, using typical values
% 
% dTdt=squeeze(nansum(nansum(nansum(Ttend.*repmat(inWag,[1 1 1 11])/86400))));
% 
% 
% figure;
% plot(advTh); hold all
% plot(advTz); plot(difTh); plot(difTz); plot(surfTtend); plot(surfQapprox)
% plot(dTdt,'linewidth',2)
% plot(advTh+advTz+difTh+difTz+surfTtend.'+surfQapprox,'--')
% legend('horizontal advection','vertical advection','horizontal diffusion','vertical diffusion','surfTtend','surfQ','dT/dt','sum(fluxes)')
% %clear cellVol4 Adv* Dif* KpT Ttend Qnet
% %save('heatBudgetTrial11day.mat','-v7.3')

%% trial single cell budget
% load('fluxesT11day.mat')
% load('edgesWAGeuler2017.mat','dInterface','XC','YC','*Coast')
% load('distancesAreas.mat','RAC','hFacC')
% 
% dZ(1,1,1:46)=diff(dInterface);
% cellVol=repmat(RAC,[1 1 46]).*hFacC.*repmat(dZ,[700 200 1]);
% cellVol4=repmat(cellVol,[1 1 1 11]);
% 
% advTh=squeeze((AdvUt(300,100,5,:)...
%     -AdvUt(301,100,5,:)+AdvVt(300,100,5,:)...
%     -AdvVt(300,101,5,:))./cellVol4(300,100,5,:));
% advTz=squeeze((AdvWt(300,100,6,:)-AdvWt(300,100,5,:))./cellVol4(300,100,5,:));
% difTz=squeeze((DifzEt(300,100,6,:)+DifIt(300,100,6,:)+KpT(300,100,6,:)-DifzEt(300,100,5,:)-DifIt(300,100,5,:)-KpT(300,100,5,:))./cellVol4(300,100,5,:));
% difTh=squeeze((DifxEt(300,100,5,:)    -DifxEt(301,100,5,:)...
%     +DifyEt(300,100,5,:)    -DifyEt(300,101,5,:))./cellVol4(300,100,5,:));
% 
% dTdt=squeeze(Ttend(300,100,5,:)./86400);
% 
% figure;
% plot(advTh); hold all
% plot(advTz); plot(difTh); plot(difTz); 
% plot(dTdt,'linewidth',2)
% plot(advTh+advTz+difTh+difTz,'--')
% legend('horizontal advection','vertical advection','horizontal diffusion','vertical diffusion','dT/dt','sum(fluxes)')

%% trial single column budget


% advTh=nansum(squeeze((AdvUt(300,100,:,:)...
%     -AdvUt(301,100,:,:)+AdvVt(300,100,:,:)...
%     -AdvVt(300,101,:,:))./cellVol4(300,100,:,:)));
% advTz=nansum(squeeze((AdvWt(300,100,2:46,:)-AdvWt(300,100,1:45,:))./cellVol4(300,100,1:45,:)));
% difTz=nansum(squeeze((DifzEt(300,100,2:46,:)+DifIt(300,100,2:46,:)+KpT(300,100,2:46,:)-DifzEt(300,100,1:45,:)-DifIt(300,100,1:45,:)-KpT(300,100,1:45,:))./cellVol4(300,100,1:45,:)));
% difTh=nansum(squeeze((DifxEt(300,100,:,:)    -DifxEt(301,100,:,:)...
%     +DifyEt(300,100,:,:)    -DifyEt(300,101,:,:))./cellVol4(300,100,:,:)));
% 
% tsurfcor=squeeze(nansum(nansum(squeeze(WT(:,:,1,:)).*repmat(RAC,[1 1 11]))))./sum(sum(RAC));
% for i=1:11
% surfTtend(i)=(tsurfcor(i)-squeeze(WT(300,100,1,i)))./(dZ(1,1,1).*hFacC(300,100,1));
% end
% surfQapprox=squeeze(Tflux(300,100,:)./(1026*4000*dZ(1,1,1).*repmat(hFacC(300,100,1),[1 1 11])));%approx bc not sure what rhoConst or Cp should be, using typical values
% 
% 
% dTdt=nansum(squeeze(Ttend(300,100,:,:)./86400));
% 
% figure;
% plot(advTh); hold all
% plot(advTz); plot(difTh); plot(difTz); plot(surfTtend); plot(surfQapprox)
% plot(dTdt,'linewidth',2)
% plot(advTh+advTz+difTh+difTz+surfTtend+surfQapprox.','--')
% legend('horizontal advection','vertical advection','horizontal diffusion','vertical diffusion','surfTtend','surfQ','dT/dt','sum(fluxes)')
% 
% sumflux=advTh+advTz+difTh+difTz+surfTtend+surfQapprox.';
% perrorT=100*(dTdt-sumflux)./dTdt;
% figure; plot(perrorT)

%% trial 2 Euler WAG: include the whole blame thing
%clear
addpath('/nobackup1/gbrett/mStuff/')

load('edgesWAGlagrange2017.mat','inWag','open*','*edge')
load('geometrySpinupSteady.mat','dInterface','Angle*','*C*')
load('distancesAreas.mat','RAC','hFacC','dxc','dyc')
DXC=reshape(dxc,[700 200]); DYC=reshape(dyc,[700 200]);
nt=146;%148 euler
days=1:146;%148 euler
Rho=999.8*ones([700 200 nt]);
load('saTCtCpSigma162NR.mat', 'CP')
cp=squeeze(CP(:,:,1,days)); clear CP


dInterface
dZ(1,1,1:46)=diff(dInterface);
cellVol=repmat(RAC,[1 1 46]).*hFacC.*repmat(dZ,[700 200 1]);
cellVol4=repmat(cellVol,[1 1 1 nt]);

euler=false; lagrange=true;

%old
%addpath('/nobackup1/gbrett/addforce/')
%load('fluxesT58and39day.mat')
%load('rhoSSH148day.mat', 'Rho')
%squeeze(Rho(:,:,1,days));
%% euler
if euler
load('fluxesAdvT148day.mat')
advTz=squeeze(nansum(nansum(nansum((AdvWt(:,:,2:46,:)-AdvWt(:,:,1:45,:)).*repmat(inWag(:,:,1:45,:),[1 1 1 nt])))));%./cellVol4(:,:,1:45,:)))));
advTh=squeeze(nansum(nansum(nansum((AdvUt(1:end-1,1:end-1,:,:)-AdvUt(2:end,1:end-1,:,:)+AdvVt(1:end-1,1:end-1,:,:)...
    -AdvVt(1:end-1,2:end,:,:)).*repmat(inWag(1:end-1,1:end-1,:,:),[1 1 1 nt])))));%./cellVol4(1:end-1,1:end-1,:,:)))));
% runningTotal=AdvWt(1:end-1,1:end-1,2:46,:)-AdvWt(1:end-1,1:end-1,1:45,:)...
%     +AdvUt(1:end-1,1:end-1,1:45,:)-AdvUt(2:end,1:end-1,1:45,:)+AdvVt(1:end-1,1:end-1,1:45,:)...
%     -AdvVt(1:end-1,2:end,1:45,:);
advTHopen=squeeze(nansum(nansum(nansum(...
    (AdvUt(1:end-1,1:end-1,:,:).*repmat(openW(1:end-1,1:end-1,:),[1 1 1 nt])-AdvUt(2:end,1:end-1,:,:).*repmat(openE(1:end-1,1:end-1,:),[1 1 1 nt])...
    +AdvVt(1:end-1,1:end-1,:,:).*repmat(openS(1:end-1,1:end-1,:),[1 1 1 nt])-AdvVt(1:end-1,2:end,:,:).*repmat(openN(1:end-1,1:end-1,:),[1 1 1 nt]))...
    )))); 
advTZopen=squeeze(nansum(nansum(nansum(...
    (AdvWt(:,:,2:46,:).*repmat(openZd(:,:,1:45),[1 1 1 nt])-AdvWt(:,:,1:45,:).*repmat(openZu(:,:,1:45),[1 1 1 nt]))...%./CellVol4(:,:,1:45,:)...
    ))));
advTBopen=squeeze(nansum(nansum(nansum(...
    (AdvUt(1:end-1,1:end-1,1:45,:).*repmat(openW(1:end-1,1:end-1,1:45),[1 1 1 nt])-AdvUt(2:end,1:end-1,1:45,:).*repmat(openE(1:end-1,1:end-1,1:45),[1 1 1 nt])...
    +AdvVt(1:end-1,1:end-1,1:45,:).*repmat(openS(1:end-1,1:end-1,1:45),[1 1 1 nt])-AdvVt(1:end-1,2:end,1:45,:).*repmat(openN(1:end-1,1:end-1,1:45),[1 1 1 nt]))...
    .*repmat(bottomedge(1:end-1,1:end-1,1:45),[1 1 1 nt])...
    +(AdvWt(1:end-1,1:end-1,2:46,:).*repmat(openZd(1:end-1,1:end-1,1:45),[1 1 1 nt])-AdvWt(1:end-1,1:end-1,1:45,:).*repmat(openZu(1:end-1,1:end-1,1:45),[1 1 1 nt]))...%./CellVol4(:,:,1:45,:)...
    )))); 
advTBopenMag=squeeze(nansum(nansum(nansum(...
    (abs(AdvUt(1:end-1,1:end-1,1:45,:)).*repmat(openW(1:end-1,1:end-1,1:45),[1 1 1 nt])+abs(AdvUt(2:end,1:end-1,1:45,:)).*repmat(openE(1:end-1,1:end-1,1:45),[1 1 1 nt])...
    +abs(AdvVt(1:end-1,1:end-1,1:45,:)).*repmat(openS(1:end-1,1:end-1,1:45),[1 1 1 nt])+abs(AdvVt(1:end-1,2:end,1:45,:)).*repmat(openN(1:end-1,1:end-1,1:45),[1 1 1 nt]))...
    .*repmat(bottomedge(1:end-1,1:end-1,1:45),[1 1 1 nt])...
    +(abs(AdvWt(1:end-1,1:end-1,2:46,:)).*repmat(openZd(1:end-1,1:end-1,1:45),[1 1 1 nt])+abs(AdvWt(1:end-1,1:end-1,1:45,:)).*repmat(openZu(1:end-1,1:end-1,1:45),[1 1 1 nt]))...%./CellVol4(:,:,1:45,:)...
    )))); 
advTHopenS=squeeze(nansum(nansum(nansum(...
    (AdvUt(1:end-1,1:end-1,:,:).*repmat(openW(1:end-1,1:end-1,:),[1 1 1 nt])-AdvUt(2:end,1:end-1,:,:).*repmat(openE(1:end-1,1:end-1,:),[1 1 1 nt])...
    +AdvVt(1:end-1,1:end-1,:,:).*repmat(openS(1:end-1,1:end-1,:),[1 1 1 nt])-AdvVt(1:end-1,2:end,:,:).*repmat(openN(1:end-1,1:end-1,:),[1 1 1 nt]))...
    .*repmat(sideedge(1:end-1,1:end-1,:),[1 1 1 nt]))))); 
advTHopenSMag=squeeze(nansum(nansum(nansum(...
    (abs(AdvUt(1:end-1,1:end-1,:,:)).*repmat(openW(1:end-1,1:end-1,:),[1 1 1 nt])+abs(AdvUt(2:end,1:end-1,:,:)).*repmat(openE(1:end-1,1:end-1,:),[1 1 1 nt])...
    +abs(AdvVt(1:end-1,1:end-1,:,:)).*repmat(openS(1:end-1,1:end-1,:),[1 1 1 nt])+abs(AdvVt(1:end-1,2:end,:,:)).*repmat(openN(1:end-1,1:end-1,:),[1 1 1 nt]))...
    .*repmat(sideedge(1:end-1,1:end-1,:),[1 1 1 nt]))))); 
clear Adv*
load('fluxesDifT148day.mat','DifzEt','DifIt')
load('fluxesOtherT148day.mat','KpT')
difTz=squeeze(nansum(nansum(nansum((DifzEt(:,:,2:46,:)+DifIt(:,:,2:46,:)+KpT(:,:,2:46,:)).*repmat(openZd(:,:,1:45,:),[1 1 1 nt])))));%./cellVol4(:,:,1:45,:)))));
difTz2=squeeze(nansum(nansum(nansum((DifzEt(:,:,2:46,:)+DifIt(:,:,2:46,:)+KpT(:,:,2:46,:)-DifzEt(:,:,1:45,:)-DifIt(:,:,1:45,:)-KpT(:,:,1:45,:)).*repmat(inWag(:,:,1:45,:),[1 1 1 nt])))));%./cellVol4(:,:,1:45,:)))));
difTZopen1=squeeze(nansum(nansum(nansum(...
      ((DifzEt(1:end-1,1:end-1,2:46,:)+DifIt(1:end-1,1:end-1,2:46,:)+KpT(1:end-1,1:end-1,2:46,:)).*repmat(openZd(1:end-1,1:end-1,1:45),[1 1 1 nt])...
      -(DifzEt(1:end-1,1:end-1,1:45,:)+DifIt(1:end-1,1:end-1,1:45,:)+KpT(1:end-1,1:end-1,1:45,:)).*repmat(openZu(1:end-1,1:end-1,1:45),[1 1 1 nt]))))));
difTZopen1Mag=squeeze(nansum(nansum(nansum(...
      ((abs(DifzEt(1:end-1,1:end-1,2:46,:))+abs(DifIt(1:end-1,1:end-1,2:46,:))+abs(KpT(1:end-1,1:end-1,2:46,:))).*repmat(openZd(1:end-1,1:end-1,1:45),[1 1 1 nt])...
      -(abs(DifzEt(1:end-1,1:end-1,1:45,:))+abs(DifIt(1:end-1,1:end-1,1:45,:))+abs(KpT(1:end-1,1:end-1,1:45,:))).*repmat(openZu(1:end-1,1:end-1,1:45),[1 1 1 nt]))))));
maxDifTZopen1=max(difTZopen1)
% runningTotal=runningTotal+DifzEt(1:end-1,1:end-1,2:46,:)+DifIt(1:end-1,1:end-1,2:46,:)...
%     +KpT(1:end-1,1:end-1,2:46,:)-DifzEt(1:end-1,1:end-1,1:45,:)-DifIt(1:end-1,1:end-1,1:45,:)-KpT(1:end-1,1:end-1,1:45,:);
clear KpT Dif*
load('fluxesDifT148day.mat','DifxEt','DifyEt')
difTHopenS=squeeze(nansum(nansum(nansum(...
    (DifxEt(1:end-1,1:end-1,1:45,:).*repmat(openW(1:end-1,1:end-1,1:45),[1 1 1 nt])-DifxEt(2:end,1:end-1,1:45,:).*repmat(openE(1:end-1,1:end-1,1:45),[1 1 1 nt])...
    +DifyEt(1:end-1,1:end-1,1:45,:).*repmat(openS(1:end-1,1:end-1,1:45),[1 1 1 nt])-DifyEt(1:end-1,2:end,1:45,:).*repmat(openN(1:end-1,1:end-1,1:45),[1 1 1 nt]))...
    .*repmat(sideedge(1:end-1,1:end-1,1:45),[1 1 1 nt])...
    ))));
difTHopenSMag=squeeze(nansum(nansum(nansum(...
    (abs(DifxEt(1:end-1,1:end-1,1:45,:)).*repmat(openW(1:end-1,1:end-1,1:45),[1 1 1 nt])+abs(DifxEt(2:end,1:end-1,1:45,:)).*repmat(openE(1:end-1,1:end-1,1:45),[1 1 1 nt])...
    +abs(DifyEt(1:end-1,1:end-1,1:45,:)).*repmat(openS(1:end-1,1:end-1,1:45),[1 1 1 nt])+abs(DifyEt(1:end-1,2:end,1:45,:)).*repmat(openN(1:end-1,1:end-1,1:45),[1 1 1 nt]))...
    .*repmat(sideedge(1:end-1,1:end-1,1:45),[1 1 1 nt])...
    ))));
difTZopen2=squeeze(nansum(nansum(nansum(...
    (DifxEt(1:end-1,1:end-1,1:45,:).*repmat(openW(1:end-1,1:end-1,1:45),[1 1 1 nt])-DifxEt(2:end,1:end-1,1:45,:).*repmat(openE(1:end-1,1:end-1,1:45),[1 1 1 nt])...
    +DifyEt(1:end-1,1:end-1,1:45,:).*repmat(openS(1:end-1,1:end-1,1:45),[1 1 1 nt])-DifyEt(1:end-1,2:end,1:45,:).*repmat(openN(1:end-1,1:end-1,1:45),[1 1 1 nt]))...
    .*repmat(bottomedge(1:end-1,1:end-1,1:45),[1 1 1 nt])...
    ))));
difTZopen2Mag=squeeze(nansum(nansum(nansum(...
    (abs(DifxEt(1:end-1,1:end-1,1:45,:)).*repmat(openW(1:end-1,1:end-1,1:45),[1 1 1 nt])+abs(DifxEt(2:end,1:end-1,1:45,:)).*repmat(openE(1:end-1,1:end-1,1:45),[1 1 1 nt])...
    +abs(DifyEt(1:end-1,1:end-1,1:45,:)).*repmat(openS(1:end-1,1:end-1,1:45),[1 1 1 nt])+abs(DifyEt(1:end-1,2:end,1:45,:)).*repmat(openN(1:end-1,1:end-1,1:45),[1 1 1 nt]))...
    .*repmat(bottomedge(1:end-1,1:end-1,1:45),[1 1 1 nt])...
    ))));
difTZopen=difTZopen1+difTZopen2;
difTh=squeeze(nansum(nansum(nansum((DifxEt(1:end-1,1:end-1,:,:)-DifxEt(2:end,1:end-1,:,:)+DifyEt(1:end-1,1:end-1,:,:)...
    -DifyEt(1:end-1,2:end,:,:)).*repmat(inWag(1:end-1,1:end-1,:,:),[1 1 1 nt])))));%./cellVol4(1:end-1,1:end-1,:,:)))));
% runningTotal=runningTotal+DifxEt(1:end-1,1:end-1,1:45,:)-DifxEt(2:end,1:end-1,1:45,:)+DifyEt(1:end-1,1:end-1,1:45,:)...
%     -DifyEt(1:end-1,2:end,1:45,:);
clear Dif*
globalArea=2.130281953410312e45;

% %this should be put back in!!
load('divTS162NR.mat','ThDiv');
% load('varying148ts16levelsRho.mat', 'T');
% T=T(:,:,1:15,1:nt);
% load('dwadz148levels15.mat','dWdZ')
% dWdZ=dWdZ(:,:,:,1:nt);
tdwdz=-squeeze(nansum(nansum(nansum(ThDiv(:,:,1:15,1:nt).*repmat(inWag(:,:,1:15,:),[1 1 1 nt])))));
advTh2=advTh+tdwdz;
advTz2=advTz-tdwdz; 
clear ThDiv

load('fluxesVelT148day.mat','WT')
tsurfcor=squeeze(nansum(nansum(squeeze(WT(:,:,1,:)).*repmat(RAC,[1 1 nt]))))./globalArea;
tsurfcor2=squeeze(nansum(nansum(squeeze(WT(:,:,1,:)).*repmat(RAC.*inWag(:,:,1),[1 1 nt]))))./globalArea;
for i=1:nt
surfTtend3(i)=nansum(nansum(inWag(:,:,1).*(-squeeze(WT(:,:,1,i))).*cellVol(:,:,1)./(dZ(1,1,1).*hFacC(:,:,1))));
surfTtend3Mag(i)=nansum(nansum(inWag(:,:,1).*(squeeze(abs(WT(:,:,1,i)))).*cellVol(:,:,1)./(dZ(1,1,1).*hFacC(:,:,1))));
end
% runningTotal(1:end-1,1:end-1,1,:)=runningTotal(1:end-1,1:end-1,1,:)-WT(1:end-1,1:end-1,1,i).*cellVol4(1:end-1,1:end-1,1,:)./(dZ(1,1,1).*hFacC(1:end-1,1:end-1,1));
clear WT
load('fluxesOtherT148day.mat')
surfQ=squeeze(nansum(nansum(Tflux(:,:,1:nt).*repmat(inWag(:,:,1).*cellVol(:,:,1),[1 1 nt])./(dZ(1,1,1).*repmat(hFacC(:,:,1),[1 1 nt]).*Rho.*cp))));
dTdt=squeeze(nansum(nansum(nansum(Ttend(:,:,:,1:nt).*cellVol4.*repmat(inWag,[1 1 1 nt])/86400))));
surfQMag=squeeze(nansum(nansum(abs(Tflux(:,:,1:nt)).*repmat(inWag(:,:,1).*cellVol(:,:,1),[1 1 nt])./(dZ(1,1,1).*repmat(hFacC(:,:,1),[1 1 nt]).*Rho.*cp))));
dTdtMag=squeeze(nansum(nansum(nansum(abs(Ttend(:,:,:,1:nt)).*cellVol4.*repmat(inWag,[1 1 1 nt])/86400))));
surfQnet=squeeze(nansum(nansum(Qnet(:,:,1:nt).*repmat(inWag(:,:,1).*cellVol(:,:,1),[1 1 nt])./(dZ(1,1,1).*repmat(hFacC(:,:,1),[1 1 nt]).*Rho.*cp))));%tflux=qnet+qprecipevap
% runningTotal=runningTotal+Tflux.*cellVol4./(dZ(1,1,1).*repmat(hFacC(:,:,1),[1 1 nt]).*Rho.*cp)...
%     -Ttend.*cellVol4./86400;
clear Ttend Qnet Tflux
%load('tsTendsnap.mat','restartT')
%reT=squeeze(nansum(nansum(nansum(restartT.*cellVol4.*repmat(inWag,[1 1 1 nt])/86400))));


sumflux=advTh+advTz+difTh+difTz2+surfTtend3.'+surfQ;%+reT;
perrorT=100*(dTdt-sumflux)./dTdt;


surfTot=surfTtend3.'+surfQ;
errorT=dTdt-(advTh+advTz+difTz2+difTh+surfTot);
days=1:148;

disp('saving ttends')
%save('runningTbudget.mat','runningTotal','-v7.3');
clear open* *edge
save('temperatureBudgetEulerWAG148dayNF.mat')%,'-v7.3')%only want vectors
end
%% lagrange
if lagrange
%     load('fluxesAdvT148day.mat')
% advTz=squeeze(nansum(nansum(nansum((AdvWt(:,:,2:46,:)-AdvWt(:,:,1:45,:)).*inWag(:,:,1:45,:)))));
% advTh=squeeze(nansum(nansum(nansum((AdvUt(1:end-1,1:end-1,:,:)-AdvUt(2:end,1:end-1,:,:)+AdvVt(1:end-1,1:end-1,:,:)...
%     -AdvVt(1:end-1,2:end,:,:)).*inWag(1:end-1,1:end-1,:,:)))));
% advTHopen=squeeze(nansum(nansum(nansum(...
%     (AdvUt(1:end-1,1:end-1,:,:).*openW(1:end-1,1:end-1,:,:)-AdvUt(2:end,1:end-1,:,:).*openE(1:end-1,1:end-1,:,:)...
%     +AdvVt(1:end-1,1:end-1,:,:).*openS(1:end-1,1:end-1,:,:)-AdvVt(1:end-1,2:end,:,:).*openN(1:end-1,1:end-1,:,:))...
%     )))); 
% advTZopen=squeeze(nansum(nansum(nansum(...
%     (AdvWt(:,:,2:46,:).*openZd(:,:,1:45,:)-AdvWt(:,:,1:45,:).*openZu(:,:,1:45,:))...
%     ))));
% advTBopen=squeeze(nansum(nansum(nansum(...
%     (AdvUt(1:end-1,1:end-1,1:45,:).*openW(1:end-1,1:end-1,1:45,:)-AdvUt(2:end,1:end-1,1:45,:).*openE(1:end-1,1:end-1,1:45,:)...
%     +AdvVt(1:end-1,1:end-1,1:45,:).*openS(1:end-1,1:end-1,1:45,:)-AdvVt(1:end-1,2:end,1:45,:).*openN(1:end-1,1:end-1,1:45,:))...
%     .*bottomedge(1:end-1,1:end-1,1:45,:)...
%     +(AdvWt(1:end-1,1:end-1,2:46,:).*openZd(1:end-1,1:end-1,1:45,:)-AdvWt(1:end-1,1:end-1,1:45,:).*openZu(1:end-1,1:end-1,1:45,:))...
%     )))); 
% advTBopenMag=squeeze(nansum(nansum(nansum(...
%     (abs(AdvUt(1:end-1,1:end-1,1:45,:)).*openW(1:end-1,1:end-1,1:45,:)+abs(AdvUt(2:end,1:end-1,1:45,:)).*openE(1:end-1,1:end-1,1:45,:)...
%     +abs(AdvVt(1:end-1,1:end-1,1:45,:)).*openS(1:end-1,1:end-1,1:45,:)+abs(AdvVt(1:end-1,2:end,1:45,:)).*openN(1:end-1,1:end-1,1:45,:))...
%     .*repmat(bottomedge(1:end-1,1:end-1,1:45),[1 1 1 nt])...
%     +(abs(AdvWt(1:end-1,1:end-1,2:46,:)).*openZd(1:end-1,1:end-1,1:45,:)+abs(AdvWt(1:end-1,1:end-1,1:45,:)).*openZu(1:end-1,1:end-1,1:45,:))...
%     )))); 
% advTHopenS=squeeze(nansum(nansum(nansum(...
%     (AdvUt(1:end-1,1:end-1,:,:).*openW(1:end-1,1:end-1,:,:)-AdvUt(2:end,1:end-1,:,:).*openE(1:end-1,1:end-1,:,:)...
%     +AdvVt(1:end-1,1:end-1,:,:).*openS(1:end-1,1:end-1,:,:)-AdvVt(1:end-1,2:end,:,:).*openN(1:end-1,1:end-1,:,:))...
%     .*sideedge(1:end-1,1:end-1,:,:))))); 
% advTHopenSMag=squeeze(nansum(nansum(nansum(...
%     (abs(AdvUt(1:end-1,1:end-1,:,:)).*openW(1:end-1,1:end-1,:,:)+abs(AdvUt(2:end,1:end-1,:,:)).*openE(1:end-1,1:end-1,:,:)...
%     +abs(AdvVt(1:end-1,1:end-1,:,:)).*openS(1:end-1,1:end-1,:,:)+abs(AdvVt(1:end-1,2:end,:,:)).*openN(1:end-1,1:end-1,:,:))...
%     .*sideedge(1:end-1,1:end-1,:,:))))); 
% clear Adv*
% load('fluxesDifT148day.mat','DifzEt','DifIt')
% Difsum=DifzEt+DifIt; clear DifzEt DifIt
% load('fluxesOtherT148day.mat','KpT')
% Difsum=Difsum+KpT; clear KpT
% disp('dif z')
% %difTz=squeeze(nansum(nansum(nansum((DifzEt(:,:,2:46,:)+DifIt(:,:,2:46,:)+KpT(:,:,2:46,:)).*openZd(:,:,1:45,:)))));%./cellVol4(:,:,1:45,:)))));
% difTz2=squeeze(nansum(nansum(nansum((Difsum(:,:,2:20,1:nt)-Difsum(:,:,1:19,1:nt)).*inWag(:,:,1:19,1:nt)))));%squeeze(nansum(nansum(nansum((DifzEt(:,:,2:46,:)+DifIt(:,:,2:46,:)+KpT(:,:,2:46,:)-DifzEt(:,:,1:45,:)-DifIt(:,:,1:45,:)-KpT(:,:,1:45,:)).*inWag(:,:,1:45,:)))));
% difTZopen1=squeeze(nansum(nansum(nansum(...
%       ((Difsum(1:end-1,1:end-1,2:20,1:nt)).*openZd(1:end-1,1:end-1,1:19,1:nt)...
%       -(Difsum(1:end-1,1:end-1,1:19,1:nt)).*openZu(1:end-1,1:end-1,1:19,1:nt))))));
% difTZopen1Mag=squeeze(nansum(nansum(nansum(...
%       ((abs(Difsum(1:end-1,1:end-1,2:20,1:nt))).*openZd(1:end-1,1:end-1,1:19,1:nt)...
%       -(abs(Difsum(1:end-1,1:end-1,1:19,1:nt))).*openZu(1:end-1,1:end-1,1:19,1:nt))))));
% maxDifTZopen1=max(difTZopen1)
% % runningTotal=runningTotal+DifzEt(1:end-1,1:end-1,2:46,:)+DifIt(1:end-1,1:end-1,2:46,:)...
% %     +KpT(1:end-1,1:end-1,2:46,:)-DifzEt(1:end-1,1:end-1,1:45,:)-DifIt(1:end-1,1:end-1,1:45,:)-KpT(1:end-1,1:end-1,1:45,:);
% clear KpT Dif*
% save('temperatureLagrangeWAGdifZ.mat','difT*')
% disp('dif z done')
% load('fluxesDifT148day.mat','DifxEt','DifyEt')
% disp('dif h')
% difTHopenS=squeeze(nansum(nansum(nansum(...
%     (DifxEt(1:end-1,1:end-1,1:19,1:nt).*openW(1:end-1,1:end-1,1:19,1:nt)-DifxEt(2:end,1:end-1,1:19,1:nt).*openE(1:end-1,1:end-1,1:19,1:nt)...
%     +DifyEt(1:end-1,1:end-1,1:19,1:nt).*openS(1:end-1,1:end-1,1:19,1:nt)-DifyEt(1:end-1,2:end,1:19,1:nt).*openN(1:end-1,1:end-1,1:19,1:nt))...
%     .*sideedge(1:end-1,1:end-1,1:19,1:nt)...
%     ))));
% difTHopenSMag=squeeze(nansum(nansum(nansum(...
%     (abs(DifxEt(1:end-1,1:end-1,1:19,1:nt)).*openW(1:end-1,1:end-1,1:19,1:nt)+abs(DifxEt(2:end,1:end-1,1:19,1:nt)).*openE(1:end-1,1:end-1,1:19,1:nt)...
%     +abs(DifyEt(1:end-1,1:end-1,1:19,1:nt)).*openS(1:end-1,1:end-1,1:19,1:nt)+abs(DifyEt(1:end-1,2:end,1:19,1:nt)).*openN(1:end-1,1:end-1,1:19,1:nt))...
%     .*sideedge(1:end-1,1:end-1,1:19,1:nt)...
%     ))));
% difTZopen2=squeeze(nansum(nansum(nansum(...
%     (DifxEt(1:end-1,1:end-1,1:19,1:nt).*openW(1:end-1,1:end-1,1:19,1:nt)-DifxEt(2:end,1:end-1,1:19,1:nt).*openE(1:end-1,1:end-1,1:19,1:nt)...
%     +DifyEt(1:end-1,1:end-1,1:19,1:nt).*openS(1:end-1,1:end-1,1:19,1:nt)-DifyEt(1:end-1,2:end,1:19,1:nt).*openN(1:end-1,1:end-1,1:19,1:nt))...
%     .*bottomedge(1:end-1,1:end-1,1:19,1:nt)...
%     ))));
% difTZopen2Mag=squeeze(nansum(nansum(nansum(...
%     (abs(DifxEt(1:end-1,1:end-1,1:19,1:nt)).*openW(1:end-1,1:end-1,1:19,1:nt)+abs(DifxEt(2:end,1:end-1,1:19,1:nt)).*openE(1:end-1,1:end-1,1:19,1:nt)...
%     +abs(DifyEt(1:end-1,1:end-1,1:19,1:nt)).*openS(1:end-1,1:end-1,1:19,1:nt)+abs(DifyEt(1:end-1,2:end,1:19,1:nt)).*openN(1:end-1,1:end-1,1:19,1:nt))...
%     .*bottomedge(1:end-1,1:end-1,1:19,1:nt)...
%     ))));
% difTZopen=difTZopen1+difTZopen2;
% difTh=squeeze(nansum(nansum(nansum((DifxEt(1:end-1,1:end-1,1:20,1:nt)-DifxEt(2:end,1:end-1,1:20,1:nt)+DifyEt(1:end-1,1:end-1,1:20,1:nt)...
%     -DifyEt(1:end-1,2:end,1:20,1:nt)).*inWag(1:end-1,1:end-1,1:20,1:nt)))));
% clear Dif*
% save('temperatureLagrangeWAGdif.mat','difT*')
disp('done dif h')
nt=146
globalArea=2.130281953410312e45;

% %this should be put back in!!
% load('divTS162NR.mat','ThDiv');
% % load('varying148ts16levelsRho.mat', 'T');
% % T=T(:,:,1:15,1:nt);
% % load('dwadz148levels15.mat','dWdZ')
% % dWdZ=dWdZ(:,:,:,1:nt);
% tdwdz=-squeeze(nansum(nansum(nansum(ThDiv(:,:,1:15,1:nt).*inWag(:,:,1:15,:)))));
% advTh2=advTh+tdwdz;
% advTz2=advTz-tdwdz; 
% clear ThDiv

load('fluxesVelT148day.mat','WT')
disp('surf')
%tsurfcor=squeeze(nansum(nansum(squeeze(WT(:,:,1,:)).*repmat(RAC,[1 1 nt]))))./globalArea;
%tsurfcor2=squeeze(nansum(nansum(squeeze(WT(:,:,1,:).*inWag(:,:,1,:)).*repmat(RAC,[1 1 nt]))))./globalArea;
for i=1:nt
surfTtend3(i)=nansum(nansum(squeeze((-inWag(:,:,1,i).*WT(:,:,1,i))).*cellVol(:,:,1)./(dZ(1,1,1).*hFacC(:,:,1))));
surfTtend3Mag(i)=nansum(nansum((squeeze(inWag(:,:,1,i).*abs(WT(:,:,1,i)))).*cellVol(:,:,1)./(dZ(1,1,1).*hFacC(:,:,1))));
end
disp('surfTtend3 done')
% runningTotal(1:end-1,1:end-1,1,:)=runningTotal(1:end-1,1:end-1,1,:)-WT(1:end-1,1:end-1,1,i).*cellVol4(1:end-1,1:end-1,1,:)./(dZ(1,1,1).*hFacC(1:end-1,1:end-1,1));
clear WT
load('fluxesOtherT148day.mat')
nt=146
sizeTflux=size(Tflux)
sizeinWag=size(inWag)
sizeCV4=size(cellVol4)
size(Rho)
size(cp)
size(hFacC)
nt
denom=dZ(1,1,1).*repmat(hFacC(1:700,1:200,1),[1 1 nt]).*Rho.*cp;
disp('denom done')
surfQ=squeeze(nansum(nansum(Tflux(1:700,1:200,1:nt).*squeeze(inWag(1:700,1:200,1,1:nt).*cellVol4(1:700,1:200,1,1:nt))./denom)));
disp('surfQ done')
surfQMag=squeeze(nansum(nansum(abs(Tflux(:,:,1:nt)).*squeeze(inWag(:,:,1,:).*cellVol4(:,:,1,:))./(denom))));
disp('surfQmag done')

sizeTtend=size(Ttend)
dTdt=squeeze(nansum(nansum(nansum(Ttend(:,:,1:20,1:nt).*cellVol4(:,:,1:20,:).*inWag/86400))));
disp('dTdt done')
dTdtMag=squeeze(nansum(nansum(nansum(abs(Ttend(:,:,1:20,1:nt)).*cellVol4(:,:,1:20,:).*inWag/86400))));
disp('dTdtmag done')

surfQnet=squeeze(nansum(nansum(Qnet(:,:,1:nt).*squeeze(inWag(:,:,1,:).*cellVol4(:,:,1,:))./(dZ(1,1,1).*repmat(hFacC(:,:,1),[1 1 nt]).*Rho.*cp))));%tflux=qnet+qprecipevap
disp('surfqnet done')

clear Ttend Qnet Tflux
save('temperatureLagrangeWAGsurface.mat','surfTtend3*','surfQ*','dTdt*','surfQnet')
disp('temperature surface done')
%sumflux=advTh+advTz+difTh+difTz2+surfTtend3.'+surfQ;
%perrorT=100*(dTdt-sumflux)./dTdt;


surfTot=surfTtend3.'+surfQ;
%errorT=dTdt-(advTh+advTz+difTz2+difTh+surfTot);
%days=1:148;

disp('saving ttends')
%save('runningTbudget.mat','runningTotal','-v7.3');
clear open* *edge
%save('temperatureBudgetLagrangeWAG146dayNF.mat')%,'-v7.3')%only want vectors
      
end
%% old junk
% advTh=squeeze(nansum(nansum(nansum((AdvUt(1:end-1,1:end-1,:,:)-AdvUt(2:end,1:end-1,:,:)+AdvVt(1:end-1,1:end-1,:,:)...
%     -AdvVt(1:end-1,2:end,:,:)).*repmat(inWag(1:end-1,1:end-1,:,:),[1 1 1 11])./cellVol4(1:end-1,1:end-1,:,:)))));

%advTh=squeeze(nansum(nansum(nansum( ((UT(1:end-1,1:end-1,:,:)-UT(2:end,1:end-1,:,:))./repmat(DXC(1:end-1,1:end-1),[1 1 46 nt])...
%    +(VT(1:end-1,1:end-1,:,:)-VT(1:end-1,2:end,:,:))./repmat(DYC(1:end-1,1:end-1),[1 1 46 nt])).*repmat(inWag(1:end-1,1:end-1,:,:),[1 1 1 nt])./repmat(hFacC(1:end-1,1:end-1,:),[1 1 1 nt])))));
%advTz=squeeze(nansum(nansum(nansum((WT(:,:,2:46,:)-WT(:,:,1:45,:)).*repmat(inWag(:,:,1:45,:),[1 1 1 nt])./repmat(dZ(:,:,1:45),[700 200 1 nt])))));

% numB=((UT(1:end-1,1:end-1,:,:)-UT(2:end,1:end-1,:,:))./repmat(DXC(1:end-1,1:end-1),[1 1 46 nt])...
%     +(VT(1:end-1,1:end-1,:,:)-VT(1:end-1,2:end,:,:))./repmat(DYC(1:end-1,1:end-1),[1 1 46 nt])).*repmat(inWag(1:end-1,1:end-1,:,:),[1 1 1 nt])./repmat(hFacC(1:end-1,1:end-1,:),[1 1 1 nt]);


%difTz3=difTz2+numZ;
%difTh3=difTh+numH;
% numH=squeeze(nansum(nansum(nansum(...
%     numA-numB...
%         ))));
%    clear numA numB
% numZ=squeeze(nansum(nansum(nansum(...
%     (AdvWt(:,:,2:46,:)-AdvWt(:,:,1:45,:)).*repmat(inWag(:,:,1:45,:),[1 1 1 nt])./cellVol4(:,:,1:45,:)...
%     -(WT(:,:,2:46,:)-WT(:,:,1:45,:)).*repmat(inWag(:,:,1:45,:),[1 1 1 nt])./repmat(dZ(:,:,1:45),[700 200 1 nt])...
%     ))));

%surfQapprox=squeeze(nansum(nansum(Tflux.*repmat(inWag(:,:,1),[1 1 45])./(3994*dZ(1,1,1).*repmat(hFacC(:,:,1),[1 1 45]).*Rho))));%approx bc using const Cp
%surfTtend(i)=nansum(nansum(inWag(:,:,1).*(tsurfcor(i).*ones([700 200])-squeeze(WT(:,:,1,i)))./(dZ(1,1,1).*hFacC(:,:,1))));
%surfTtend2(i)=nansum(nansum(inWag(:,:,1).*(tsurfcor2(i).*ones([700 200])-squeeze(WT(:,:,1,i)))./(dZ(1,1,1).*hFacC(:,:,1))));

%short1(1,1,1:46)=0.62*exp(dInterface(1:46)/0.6)+0.38*exp(dInterface(1:46)/20);%no shortwave!
%short2(1,1,1:46)=0.62*exp(dInterface(2:47)/0.6)+0.38*exp(dInterface(2:47)/20);
%shortwaveApprox=Qsw/(1026*4000)/(repmat(dZ,[700 200 1 45]).*repmat(hFacC,[1 1 1 45]))*repmat(short1-short2,[700 200 1 45]);
% advTz2=advTz+numZ-tdwdz;
% advTh2=advTh+tdwdz;
%% plots
% load('temperatureBudgetEulerWAG148dayNF.mat')
% 
% % %this one!
% % figure;
% % plot(advTh+advTz); hold all
% % %plot(advTz);
% % plot(difTh); plot(difTz2); plot(surfTtend3+surfQ.')
% % plot(-dTdt,'linewidth',2)
% % plot(advTh+advTz+difTh+difTz2+surfTtend3.'+surfQ-dTdt,'--','linewidth',2)
% % legend('advection','horizontal diffusion','vertical diffusion','surface','-dT/dt','sum(fluxes)')
% % 
% %surfaceT=-(advTHopenS+advTBopen+difTHopenS+difTZopen)+dTdt;
surfaceT=surfQ+surfTtend3.';
figure; plot(advTHopenS,'linewidth',2); hold all
plot(advTBopen,'linewidth',2); 
plot(difTZopen,'linewidth',2); 
plot(surfaceT,'linewidth',2)
plot(-dTdt,'linewidth',2)
plot(difTHopenS,'linewidth',2); %plot(difZopen)
plot(advTHopenS+advTBopen+difTHopenS+difTZopen+surfaceT-dTdt,'k--','linewidth',2)
xlabel('days','fontsize',20); ylabel('potential temperature flux, \circ C m^3/s','fontsize',20)
legend('advection through sides','advection through bottom','diffusion through bottom','surface incl SSH','dT/dt','diffusion through sides','sum(fluxes)')
title('Temperature budget, Euler WAG','fontsize',22)
set(gca,'fontsize',20)

rhocp=1000*4000; T0=273.15;
figure; plot(rhocp*(advTHopenS+surfTtend3.'),'linewidth',2); hold all
plot(rhocp*advTBopen,'linewidth',2); 
plot(rhocp*surfQ,'linewidth',2)
plot(-rhocp*dTdt,'linewidth',2)
plot(rhocp*difTZopen,'linewidth',2); 
plot(rhocp*difTHopenS,'linewidth',2); %plot(difZopen)
plot(rhocp*(advTHopenS+advTBopen+difTHopenS+difTZopen+surfaceT-dTdt),'k--','linewidth',2)
xlabel('days','fontsize',20); ylabel('Heat transports, J/s','fontsize',20)
legend('advection through sides','advection through bottom','surface','-dH/dt','diffusion through bottom','diffusion through sides','sum')
title('Heat budget, Euler WAG','fontsize',22)
set(gca,'fontsize',18)
%%
%rho0=1000; 
tbar=-16.99;
%dvdtPerfect=(vfluxhs+vfluxb+precip);
dvdt=[volumeSSHflux;0];
vcorrection=-(vfluxhs+vfluxb+precip-[volumeSSHflux;0]);
figure;
plot(rhocp*(advTHopenS+surfTtend3.')+rhocp*tbar*(vfluxhs+vcorrection),'linewidth',2); hold all 
plot(advTBopen*rhocp+rhocp*tbar*vfluxb,'linewidth',2)%-advSZu
%plot(advSZu,'linewidth',2)%,'advection through top'
plot(rhocp*difTZopen,'linewidth',2); plot(rhocp*(surfQ+precip*tbar),'linewidth',2)
plot(-rhocp*(dTdt+tbar*dvdt),'linewidth',2)
plot(rhocp*difTHopenS,'linewidth',2);
plot(rhocp*(advTHopenS+advTBopen+difTHopenS+surfTtend3.'+surfQ-dTdt)+rhocp*tbar*(vfluxhs+vcorrection+vfluxb+precip-dvdt),'k--','linewidth',4)%+difSZopen
legend('advection through sides','advection through bottom','diffusion through bottom','surface','-dH/dt','diffusion through sides','sum')
xlabel('days','fontsize',20)
axis tight
ylabel('integrated changes in heat, J/s','fontsize',20)
title('Euler WAG Heat Budget with T=16.99 removed','fontsize',22)
set(gca,'fontsize',18)
%%
figure;
plot(rhocp*(advTHopenS+surfTtend3.')+rhocp*tbar*(vfluxhs+vcorrection)+advTBopen*rhocp+rhocp*tbar*vfluxb,'linewidth',2); hold all 
plot(rhocp*difTZopen,'linewidth',2); plot(rhocp*(surfQ+precip*tbar),'linewidth',2)
plot(-rhocp*(dTdt+tbar*dvdt),'linewidth',2)
plot(rhocp*difTHopenS,'linewidth',2);
plot(rhocp*(advTHopenS+advTBopen+difTHopenS+surfTtend3.'+surfQ-dTdt)+rhocp*tbar*(vfluxhs+vcorrection+vfluxb+precip-dvdt),'k--','linewidth',4)%+difSZopen
legend('advection','diffusion through bottom','surface','-dH/dt','diffusion through sides','sum')
xlabel('days','fontsize',20)
axis tight
ylabel('integrated changes in heat, J/s','fontsize',20)
title('Euler WAG Heat Budget with T=16.99 removed','fontsize',22)
set(gca,'fontsize',18)
%%
figure;
plot(rhocp*(advTHopenS+surfTtend3.')+rhocp*tbar*(vfluxhs+vcorrection)+advTBopen*rhocp+rhocp*tbar*vfluxb-rhocp*(dTdt+tbar*dvdt),'linewidth',2); hold all 
plot(rhocp*difTZopen,'linewidth',2); plot(rhocp*(surfQ+precip*tbar),'linewidth',2)
plot(rhocp*difTHopenS,'linewidth',2);
plot(rhocp*(advTHopenS+advTBopen+difTHopenS+surfTtend3.'+surfQ-dTdt)+rhocp*tbar*(vfluxhs+vcorrection+vfluxb+precip-dvdt)...
    -mean(rhocp*(advTHopenS+advTBopen+difTHopenS+surfTtend3.'+surfQ-dTdt)+rhocp*tbar*(vfluxhs+vcorrection+vfluxb+precip-dvdt)),'k--','linewidth',4)%+difSZopen
legend('advection-dH/dt','diffusion through bottom','surface','diffusion through sides','sum')
xlabel('days','fontsize',20)
axis tight
ylabel('integrated changes in heat, J/s','fontsize',20)
title('Euler WAG Heat Budget with T=16.99 removed','fontsize',22)
set(gca,'fontsize',18)
%%
figure; plot(rhocp*(advTHopenS+surfTtend3.')+rhocp*T0*vfluxhs,'linewidth',2); hold all
plot(rhocp*advTBopen+rhocp*T0*vfluxb,'linewidth',2); 
plot(rhocp*difTZopen,'linewidth',2); 
plot(rhocp*surfQ,'linewidth',2)
plot(-rhocp*dTdt,'linewidth',2)
plot(rhocp*difTHopenS,'linewidth',2); %plot(difZopen)
plot(rhocp*(advTHopenS+advTBopen+difTHopenS+difTZopen+surfaceT-dTdt+rhocp*T0*vfluxhs+rhocp*T0*vfluxb),'k--','linewidth',2)
xlabel('days','fontsize',20); ylabel('Heat transports, J/s','fontsize',20)
legend('advection through sides','advection through bottom','diffusion through bottom','surface','dH/dt','diffusion through sides','sum(fluxes)')
title('Heat budget, Euler WAG','fontsize',22)
set(gca,'fontsize',20)
% 
% figure; plot(advTHopenSMag,'linewidth',2); hold all
% plot(advTBopenMag,'linewidth',2); 
% plot(difTZopen1Mag+difTZopen2Mag,'linewidth',2); 
% plot(surfQMag,'linewidth',2); 
% plot(surfTtend3Mag,'linewidth',2); 
% plot(dTdtMag,'linewidth',2)
% plot(difTHopenSMag,'linewidth',2); %plot(difZopen)
% title('Temperature flux magnitudes, Euler WAG','fontsize',22)
% xlabel('days','fontsize',20); ylabel('potential temperature flux, \circ C m^3/s','fontsize',20)
% legend('advection through sides','advection through bottom','diffusion through bottom','surface','SSH','dT/dt','diffusion through sides')
% set(gca,'YScale','log')
% set(gca,'fontsize',20)
% 
% figure; plot(advTHopenS+advTBopen+surfTtend3.','linewidth',2); hold all %plot(advTBopen);
% plot(nan,nan)
%  %plot(difTHopenS,'linewidth',2); %plot(difZopen)
% plot(difTZopen,'linewidth',2); plot(surfQ,'linewidth',2); plot(-dTdt,'linewidth',2)
% plot(advTHopenS+advTBopen+difTHopenS+difTZopen+surfaceT-dTdt,'k--','linewidth',2)
% xlabel('days','fontsize',20); ylabel('potential temperature flux, \circ C m^3/s','fontsize',20)
% legend('advection incl SSH','diffusion through bottom','surface','dT/dt','sum(fluxes)')
% title('temperature budget, Euler WAG','fontsize',22)
% set(gca,'fontsize',20)
% 
% figure; plot(advTHopenS+advTBopen-dTdt+surfTtend3.','linewidth',2); hold all %plot(advTBopen);
%  %plot(difTHopenS,'linewidth',2); %plot(difZopen)
%  plot(nan,nan)
% plot(difTZopen,'linewidth',2); 
% plot(surfQ,'linewidth',2)%plot(surfaceT,'linewidth',2); %plot(-dTdt,'linewidth',2)
% plot(advTHopenS+advTBopen+difTHopenS+difTZopen+surfaceT-dTdt,'k--','linewidth',2)
% xlabel('days','fontsize',20); ylabel('potential temperature flux, \circ C m^3/s','fontsize',20)
% legend('advection incl SSH-dTdt','diffusion through bottom','surface','sum(fluxes)')
% title('temperature budget, Euler WAG','fontsize',22)
% set(gca,'fontsize',20)
% 
meanTermsC={'adv','h diff','z diff','surface','-dTdt','total'};
meanTerms=[mean(advTHopenS+advTBopen+surfTtend3.'),mean(difTHopenS),mean(difTZopen),mean(surfQ),mean(-dTdt),mean(advTHopenS+advTBopen+difTHopenS+difTZopen+surfaceT-dTdt)];
figure; bar(meanTerms)
set(gca,'XTickLabel',meanTermsC)
title('Mean Euler WAG Temperature Fluxes','fontsize',14)
ylabel('\circ C m^3/s','fontsize',14)
set(gca,'fontsize',14)

meanTermsC={'adv','h diff','z diff','surface','-dHdt','total'};
meanTerms=[mean(rhocp*(advTHopenS+advTBopen+surfTtend3.')),mean(rhocp*difTHopenS),mean(rhocp*difTZopen),mean(rhocp*surfQ),mean(-rhocp*dTdt),mean(rhocp*(advTHopenS+advTBopen+difTHopenS+difTZopen+surfaceT-dTdt))];
figure; bar(meanTerms)
set(gca,'XTickLabel',meanTermsC)
title('Mean Euler WAG Heat Transports','fontsize',14)
ylabel('J/s','fontsize',14)
set(gca,'fontsize',14)

meanTermsC={'adv','h diff','z diff','surface','-dHdt','total'};
meanTerms=[mean(rhocp*(advTHopenS+surfTtend3.')+rhocp*tbar*(vfluxhs+vcorrection)+advTBopen*rhocp+rhocp*tbar*vfluxb),...
    mean(rhocp*difTHopenS),mean(rhocp*difTZopen),mean(rhocp*(surfQ+precip*tbar)),mean(-rhocp*(dTdt+tbar*dvdt)),...
    mean(rhocp*(advTHopenS+advTBopen+difTHopenS+difTZopen+surfQ+surfTtend3.'-dTdt))];
figure; bar(meanTerms)
set(gca,'XTickLabel',meanTermsC)
title('Mean Euler WAG Heat Transports','fontsize',14)
ylabel('J/s','fontsize',14)
set(gca,'fontsize',14)

% 
% % meanTermsC={'h adv','z adv','h diff','z diff','surface','-dTdt','total'};
% % meanTerms=[mean(advTHopenS),mean(advTBopen),mean(difTHopenS),mean(difTZopen),mean(surfaceT),mean(-dTdt),mean(advTHopenS+advTBopen+difTHopenS+difTZopen+surfaceT-dTdt)];
% % figure; bar(meanTerms)
% % set(gca,'XTickLabel',meanTermsC)
% % title('Mean Euler WAG Temperature Fluxes')
% 
% meanTermsC={'adv side','adv bottom','h diff','z diff','surface','SSH','dTdt'};
% meanTerms=[mean(advTHopenSMag),mean(advTBopenMag),mean(difTHopenSMag),mean(difTZopen1Mag+difTZopen2Mag),mean(surfQMag),mean(surfTtend3Mag),mean(dTdtMag)];
% figure; bar(meanTerms)
% set(gca,'XTickLabel',meanTermsC)
% title('Mean Euler WAG Temperature Flux Magnitudes','fontsize',14)
% set(gca,'yscale','log')
% ylabel('\circ C m^3/s','fontsize',14)
% set(gca,'fontsize',14)
% 
% % figure; plot(advTh2); hold all
% % plot(advTz2); plot(difTh); plot(difTz2);
% % plot(surfTot); plot(dTdt,'linewidth',2)
% % plot(sumflux,'--'); 
% % legend('horizontal advection','vertical advection','horizontal diffusion','vertical diffusion','surface','dTdt','total fluxes')
% % 
% % figure; plot(errorT); title('error T')
% % 
% % figure; plot(-100*(dTdt-(advTh2+advTz2+difTz2+difTh+surfTtend3.'+surfQ))./dTdt); title('percent error T')
% % 
% % figure; plot(advTh2); hold all
% % plot(advTz2); plot(difTh); plot(difTz2);
% % plot(surfTot); plot(-dTdt,'linewidth',2)
% % plot(sumflux-dTdt,'r--','linewidth',2); 
% % legend('horizontal advection','vertical advection','horizontal diffusion','vertical diffusion','surface','-dTdt','sum')
% % % 
% % sumflux=advTh+advTz+difTh+difTz2+surfTtend2.'+surfQ;
% % perrorT=100*(dTdt-sumflux)./dTdt;
% % %figure; plot(perrorT) %less than 1% error for times 5-11
% % % sumflux2=advTh+advTz+difTh+difTz2+surfTtend3.'+surfQ;
% % % perrorT2=100*(dTdt-sumflux2)./dTdt;
% % % figure; plot(perrorT); hold all; plot(perrorT2)
% % sumflux3=advTh+advTz+difTh+difTz3+surfTtend3.'+surfQ;
% % perrorT2=100*(dTdt-sumflux2)./dTdt;
% % figure; plot(perrorT); hold all; plot(perrorT2)
% % 
% 
% % figure;
% % plot(advTh); hold all
% % plot(advTz); 
% % plot(difTh); %plot(difTz3); 
% % plot(dTdt-(advTh+advTz+difTh+surfTtend3.'+surfQ),'--')
% % plot(surfTtend3.'+surfQ)
% % plot(dTdt,'linewidth',2)
% % plot(advTh+advTz+difTh+surfTtend3.'+surfQ,'--','linewidth',2)
% % 
% % legend('horiz advection','vertical adv','horizontal diffusion','assumed vert diff','surf','dT/dt','sum(fluxes)')
% 
% 
% % 
% % figure;
% % plot(advTh+advTz); hold all
% % %plot(advTz); 
% % plot(difTh); plot(difTz2); plot(surfTtend3+surfQ.')
% % plot(dTdt,'linewidth',2)
% % plot(sumflux,'--','linewidth',2)
% % legend('advection','horizontal diffusion','vertical diffusion','surface','dT/dt','sum(fluxes)')
% 
% % 
% % figure;
% % plot(advTh+advTz+numH+numZ); hold all 
% % plot(dTdt-(advTh+advTz+numH+numZ+difTz2+surfTtend3.'+surfQ));
% % plot(difTz2); plot(surfTtend3.'+surfQ)
% % plot(-dTdt,'linewidth',2)
% % plot(sumflux-dTdt,'--','linewidth',2)
% % legend('advection','horizontal diffusion','vertical diffusion','surface','-dT/dt','sum')

%% new euler boundary
% load('steady148tswRho.mat','Sal','Sigma');
% load('wagBoundarySteadySurfaceInt14.mat')
% inWagH=Sal(:,:,1)<36.475 &Sal(:,:,1)>20 &XC<-3.05 &XC>-5.4;
% inWagZ=Sigma<27.5;
% inWag=inWagZ.*repmat(inWagH,[1 1 46]);
% 
% figure; [c1,h1]=contourf(XC,YC,Sal(:,:,1),36.2:0.05:36.9);
% %clabel(c1,h1,'labelspacing',130)
% hold all; plot(lontrF,lattrF,'r','linewidth',2)
%  plot(lontrB,lattrB,'b','linewidth',2)
% contour(XC,YC,inWagH,[1 1],'m--','linewidth',2)
% caxis([36.2 36.8])
% axis([-6.5 -1 34.9 37.2])
% legend('mean salinity','mean unstable manifold','mean stable manifold','Euler WAG boundary')
% title('Surface WAG and Salinity','fontsize',18)
% xlabel('longitude','fontsize',18)
% ylabel('latitude','fontsize',18)
% set(gca,'fontsize',16)

%% trial 2 Euler WAG: salinity
clear
addpath('/nobackup1/gbrett/mStuff/')
%addpath('/nobackup1/gbrett/addforce/')
%load('fluxesS58and39day.mat')
load('edgesWAGeuler2017NF.mat','inWag','open*','*edge','XC','YC')
% inWag=zeros(size(inWag));
% inWag(350:351,100:101,5:6)=1;
% openW=inWag; openE=inWag; openZd=inWag; openZu=inWag;
% openS=inWag; openN=inWag;
% openZd(350:351,100:101,6)=1; 
% openZu(350:351,100:101,5)=1;
% openW(350,100:101,5:6)=1;
% openE(351,100:101,5:6)=1;
% openS(350:351,100,5:6)=1;
% openN(350:351,101,5:6)=1;
%%
load('geometrySpinupSteady.mat','dInterface','*Coast')
load('distancesAreas.mat','RAC','hFac*','dxg','dyg')
DXG=repmat(reshape(dxg,[700 200]),[1 1 46]); DYG=repmat(reshape(dyg,[700 200]),[1 1 46]);
days=1:148;
nt=148;

% load('varying148ts16levelsRho.mat', 'Rho')
% Rho=squeeze(Rho(:,:,1,days));

%load('rhoSSH148day.mat','Rho','SSH');
%Rho1=squeeze(Rho(:,:,1,days));

%load('tsPress162NR.mat','Rho0')
Rho1=999.8*ones([700 200 nt]);%rhoConst
%Rho=repmat(Rho0,[700 200 1 nt]);

dZ(1,1,1:46)=diff(dInterface);
dZ4=repmat(dZ,[700 200 1 nt]);
cellVol=repmat(RAC,[1 1 46]).*hFacC.*repmat(dZ,[700 200 1]);
cellVol4=repmat(cellVol,[1 1 1 nt]);
%cellVol4(:,:,1,:)=repmat(RAC,[1 1 nt]).*(SSH+dZ(1));
%% volume integral of salinity, units psu m^3/s
load('fluxesAdvS148day.mat') %matches for 8-cell block
advSh=squeeze(nansum(nansum(nansum(...
    (AdvUs(1:end-1,1:end-1,:,:)-AdvUs(2:end,1:end-1,:,:)+AdvVs(1:end-1,1:end-1,:,:)...
    -AdvVs(1:end-1,2:end,:,:)).*repmat(inWag(1:end-1,1:end-1,:),[1 1 1 nt])...%./cellVol4(1:end-1,1:end-1,:,:)...
    ))));
advSz=squeeze(nansum(nansum(nansum(...
    (AdvWs(:,:,2:46,:)-AdvWs(:,:,1:45,:)).*repmat(inWag(:,:,1:45),[1 1 1 nt])...%./cellVol4(:,:,1:45,:)...
    ))));
advSHopen=squeeze(nansum(nansum(nansum(...
    (AdvUs(1:end-1,1:end-1,:,:).*repmat(openW(1:end-1,1:end-1,:),[1 1 1 nt])-AdvUs(2:end,1:end-1,:,:).*repmat(openE(1:end-1,1:end-1,:),[1 1 1 nt])...
    +AdvVs(1:end-1,1:end-1,:,:).*repmat(openS(1:end-1,1:end-1,:),[1 1 1 nt])-AdvVs(1:end-1,2:end,:,:).*repmat(openN(1:end-1,1:end-1,:),[1 1 1 nt]))...
    )))); %./nansum(nansum(nansum(cellVol4(1:end-1,1:end-1,:,:).*repmat(inWag(1:end-1,1:end-1,:),[1 1 1 nt])...
   % ))));
advSZopen=squeeze(nansum(nansum(nansum(...
    (AdvWs(:,:,2:46,:).*repmat(openZd(:,:,1:45),[1 1 1 nt])-AdvWs(:,:,1:45,:).*repmat(openZu(:,:,1:45),[1 1 1 nt]))...%./CellVol4(:,:,1:45,:)...
    ))));%./nansum(nansum(nansum(cellVol4(:,:,1:45,:).*repmat(inWag(:,:,1:45),[1 1 1 nt])...    
    %))));
advSBopen=squeeze(nansum(nansum(nansum(...
    (AdvUs(1:end-1,1:end-1,1:45,:).*repmat(openW(1:end-1,1:end-1,1:45),[1 1 1 nt])-AdvUs(2:end,1:end-1,1:45,:).*repmat(openE(1:end-1,1:end-1,1:45),[1 1 1 nt])...
    +AdvVs(1:end-1,1:end-1,1:45,:).*repmat(openS(1:end-1,1:end-1,1:45),[1 1 1 nt])-AdvVs(1:end-1,2:end,1:45,:).*repmat(openN(1:end-1,1:end-1,1:45),[1 1 1 nt]))...
    .*repmat(bottomedge(1:end-1,1:end-1,1:45),[1 1 1 nt])...
    +(AdvWs(1:end-1,1:end-1,2:46,:).*repmat(openZd(1:end-1,1:end-1,1:45),[1 1 1 nt])-AdvWs(1:end-1,1:end-1,1:45,:).*repmat(openZu(1:end-1,1:end-1,1:45),[1 1 1 nt]))...%./CellVol4(:,:,1:45,:)...
    )))); 
advSBopenMag=squeeze(nansum(nansum(nansum(...
    (abs(AdvUs(1:end-1,1:end-1,1:45,:)).*repmat(openW(1:end-1,1:end-1,1:45),[1 1 1 nt])+abs(AdvUs(2:end,1:end-1,1:45,:)).*repmat(openE(1:end-1,1:end-1,1:45),[1 1 1 nt])...
    +abs(AdvVs(1:end-1,1:end-1,1:45,:)).*repmat(openS(1:end-1,1:end-1,1:45),[1 1 1 nt])+abs(AdvVs(1:end-1,2:end,1:45,:)).*repmat(openN(1:end-1,1:end-1,1:45),[1 1 1 nt]))...
    .*repmat(bottomedge(1:end-1,1:end-1,1:45),[1 1 1 nt])...
    +(abs(AdvWs(1:end-1,1:end-1,2:46,:)).*repmat(openZd(1:end-1,1:end-1,1:45),[1 1 1 nt])+abs(AdvWs(1:end-1,1:end-1,1:45,:)).*repmat(openZu(1:end-1,1:end-1,1:45),[1 1 1 nt]))...%./CellVol4(:,:,1:45,:)...
    )))); 
advSHopenS=squeeze(nansum(nansum(nansum(...
    (AdvUs(1:end-1,1:end-1,:,:).*repmat(openW(1:end-1,1:end-1,:),[1 1 1 nt])-AdvUs(2:end,1:end-1,:,:).*repmat(openE(1:end-1,1:end-1,:),[1 1 1 nt])...
    +AdvVs(1:end-1,1:end-1,:,:).*repmat(openS(1:end-1,1:end-1,:),[1 1 1 nt])-AdvVs(1:end-1,2:end,:,:).*repmat(openN(1:end-1,1:end-1,:),[1 1 1 nt]))...
    .*repmat(sideedge(1:end-1,1:end-1,:,:),[1 1 1 nt])...
    )))); 
advSHopenSMag=squeeze(nansum(nansum(nansum(...
    (abs(AdvUs(1:end-1,1:end-1,:,:)).*repmat(openW(1:end-1,1:end-1,:),[1 1 1 nt])+abs(AdvUs(2:end,1:end-1,:,:)).*repmat(openE(1:end-1,1:end-1,:),[1 1 1 nt])...
    +abs(AdvVs(1:end-1,1:end-1,:,:)).*repmat(openS(1:end-1,1:end-1,:),[1 1 1 nt])+abs(AdvVs(1:end-1,2:end,:,:)).*repmat(openN(1:end-1,1:end-1,:),[1 1 1 nt]))...
    .*repmat(sideedge(1:end-1,1:end-1,:,:),[1 1 1 nt])...
    )))); 

% advSaltH=squeeze(nansum(nansum(nansum(...
%     (AdvUs(1:end-1,1:end-1,:,:)-AdvUs(2:end,1:end-1,:,:)+AdvVs(1:end-1,1:end-1,:,:)...
%     -AdvVs(1:end-1,2:end,:,:)).*Rho(1:end-1,1:end-1,:,:).*repmat(inWag(1:end-1,1:end-1,:),[1 1 1 nt])...
%     ))));
% advSaltZ=squeeze(nansum(nansum(nansum(...
%     (AdvWs(:,:,2:46,:)-AdvWs(:,:,1:45,:)).*Rho(:,:,1:45,:).*repmat(inWag(:,:,1:45),[1 1 1 nt])...
%     ))));
% advSaltHopen=squeeze(nansum(nansum(nansum(...
%     (AdvUs(1:end-1,1:end-1,:,:).*repmat(openW(1:end-1,1:end-1,:),[1 1 1 nt])-AdvUs(2:end,1:end-1,:,:).*repmat(openE(1:end-1,1:end-1,:),[1 1 1 nt])...
%     +AdvVs(1:end-1,1:end-1,:,:).*repmat(openS(1:end-1,1:end-1,:),[1 1 1 nt])-AdvVs(1:end-1,2:end,:,:).*repmat(openN(1:end-1,1:end-1,:),[1 1 1 nt]))...
%     .*Rho(1:end-1,1:end-1,:,:)...
%     ))));
% advSaltZopen=squeeze(nansum(nansum(nansum(...
%     (AdvWs(:,:,2:46,:).*repmat(openZd(:,:,1:45),[1 1 1 nt])-AdvWs(:,:,1:45,:).*repmat(openZu(:,:,1:45),[1 1 1 nt])).*Rho(:,:,1:45,:)...
%     ))));
% advSaltHopenS=squeeze(nansum(nansum(nansum(...
%     (AdvUs(1:end-1,1:end-1,:,:).*repmat(openW(1:end-1,1:end-1,:).*sideedge(1:end-1,1:end-1,:),[1 1 1 nt])-AdvUs(2:end,1:end-1,:,:).*repmat(openE(1:end-1,1:end-1,:).*sideedge(1:end-1,1:end-1,:),[1 1 1 nt])...
%     +AdvVs(1:end-1,1:end-1,:,:).*repmat(openS(1:end-1,1:end-1,:).*sideedge(1:end-1,1:end-1,:),[1 1 1 nt])-AdvVs(1:end-1,2:end,:,:).*repmat(openN(1:end-1,1:end-1,:).*sideedge(1:end-1,1:end-1,:),[1 1 1 nt]))...
%     .*Rho(1:end-1,1:end-1,:,:)...
%     ))));
% advSaltBopen=squeeze(nansum(nansum(nansum(...
%     (AdvWs(1:end-1,1:end-1,2:46,:).*repmat(openZd(1:end-1,1:end-1,1:45),[1 1 1 nt])-AdvWs(1:end-1,1:end-1,1:45,:).*repmat(openZu(1:end-1,1:end-1,1:45),[1 1 1 nt])).*Rho(1:end-1,1:end-1,1:45,:)...
%     +(AdvUs(1:end-1,1:end-1,1:45,:).*repmat(openW(1:end-1,1:end-1,1:45).*bottomedge(1:end-1,1:end-1,1:45),[1 1 1 nt])-AdvUs(2:end,1:end-1,1:45,:).*repmat(openE(1:end-1,1:end-1,1:45).*bottomedge(1:end-1,1:end-1,1:45),[1 1 1 nt])...
%     +AdvVs(1:end-1,1:end-1,1:45,:).*repmat(openS(1:end-1,1:end-1,1:45).*bottomedge(1:end-1,1:end-1,1:45),[1 1 1 nt])-AdvVs(1:end-1,2:end,1:45,:).*repmat(openN(1:end-1,1:end-1,1:45).*bottomedge(1:end-1,1:end-1,1:45),[1 1 1 nt]))...
%     .*Rho(1:end-1,1:end-1,1:45,:)...
%     ))));
clear Adv*
disp('s adv done')

load('fluxesDifS148day.mat','DifEs','DifIs')
load('fluxesOtherS148day.mat','KpS')
%difSz=squeeze(nansum(nansum(nansum((DifEs(:,:,2:46,:)+DifIs(:,:,2:46,:)+KpS(:,:,2:46,:)).*repmat(openZd(:,:,1:45,:),[1 1 1 nt])./cellVol4(:,:,1:45,:)))));
difSz2=squeeze(nansum(nansum(nansum((DifEs(:,:,2:46,:)+DifIs(:,:,2:46,:)+KpS(:,:,2:46,:)-DifEs(:,:,1:45,:)-DifIs(:,:,1:45,:)-KpS(:,:,1:45,:)).*repmat(inWag(:,:,1:45),[1 1 1 nt])))));%./cellVol4(:,:,1:45,:)))));
difSZopen1=squeeze(nansum(nansum(nansum(...
      ((DifEs(1:end-1,1:end-1,2:46,:)+DifIs(1:end-1,1:end-1,2:46,:)+KpS(1:end-1,1:end-1,2:46,:)).*repmat(openZd(1:end-1,1:end-1,1:45),[1 1 1 nt])...
      -(DifEs(1:end-1,1:end-1,1:45,:)+DifIs(1:end-1,1:end-1,1:45,:)+KpS(1:end-1,1:end-1,1:45,:)).*repmat(openZu(1:end-1,1:end-1,1:45),[1 1 1 nt]))))));
difSZopen1Mag=squeeze(nansum(nansum(nansum(...
      ((abs(DifEs(1:end-1,1:end-1,2:46,:))+abs(DifIs(1:end-1,1:end-1,2:46,:))+abs(KpS(1:end-1,1:end-1,2:46,:))).*repmat(openZd(1:end-1,1:end-1,1:45),[1 1 1 nt])...
      +(abs(DifEs(1:end-1,1:end-1,1:45,:))+abs(DifIs(1:end-1,1:end-1,1:45,:))+abs(KpS(1:end-1,1:end-1,1:45,:))).*repmat(openZu(1:end-1,1:end-1,1:45),[1 1 1 nt]))))));
% difSaltZ2=squeeze(nansum(nansum(nansum((DifEs(:,:,2:46,:)+DifIs(:,:,2:46,:)+KpS(:,:,2:46,:)-DifEs(:,:,1:45,:)-DifIs(:,:,1:45,:)-KpS(:,:,1:45,:)).*Rho(:,:,1:45,:).*repmat(inWag(:,:,1:45),[1 1 1 nt])))));
% difSaltBopen1=squeeze(nansum(nansum(nansum(...
%       ((DifEs(1:end-1,1:end-1,2:46,:)+DifIs(1:end-1,1:end-1,2:46,:)+KpS(1:end-1,1:end-1,2:46,:)).*repmat(openZd(1:end-1,1:end-1,1:45),[1 1 1 nt])...
%       -(DifEs(1:end-1,1:end-1,1:45,:)+DifIs(1:end-1,1:end-1,1:45,:)+KpS(1:end-1,1:end-1,1:45,:)).*repmat(openZu(1:end-1,1:end-1,1:45),[1 1 1 nt])).*Rho(1:end-1,1:end-1,1:45,:)))));
% 
%     (DifEs(1:end-1,1:end-1,2:46,:).*repmat(openZd(1:end-1,1:end-1,1:45),[1 1 1 nt])+DifIs(1:end-1,1:end-1,2:46,:).*repmat(openZd(1:end-1,1:end-1,1:45),[1 1 1 nt])+KpS(1:end-1,1:end-1,2:46,:).*repmat(openZd(1:end-1,1:end-1,1:45),[1 1 1 nt])...
%     -DifEs(1:end-1,1:end-1,1:45,:).*repmat(openZu(1:end-1,1:end-1,1:45),[1 1 1 nt])-DifIs(1:end-1,1:end-1,1:45,:).*repmat(openZu(1:end-1,1:end-1,1:45),[1 1 1 nt])-KpS(1:end-1,1:end-1,1:45,:).*repmat(openZu(1:end-1,1:end-1,1:45),[1 1 1 nt])).*Rho(1:end-1,1:end-1,1:45,:)))));
clear Dif* Kp*
load('fluxesDifS148day.mat','DifxEs','DifyEs')
difSh=squeeze(nansum(nansum(nansum((DifxEs(1:end-1,1:end-1,:,:)-DifxEs(2:end,1:end-1,:,:)+DifyEs(1:end-1,1:end-1,:,:)...
    -DifyEs(1:end-1,2:end,:,:)).*repmat(inWag(1:end-1,1:end-1,:),[1 1 1 nt])))));%./cellVol4(1:end-1,1:end-1,:,:)))));
difSHopen=squeeze(nansum(nansum(nansum(...
    (DifxEs(1:end-1,1:end-1,1:45,:).*repmat(openW(1:end-1,1:end-1,1:45),[1 1 1 nt])-DifxEs(2:end,1:end-1,1:45,:).*repmat(openE(1:end-1,1:end-1,1:45),[1 1 1 nt])...
    +DifyEs(1:end-1,1:end-1,1:45,:).*repmat(openS(1:end-1,1:end-1,1:45),[1 1 1 nt])-DifyEs(1:end-1,2:end,1:45,:).*repmat(openN(1:end-1,1:end-1,1:45),[1 1 1 nt]))...
    .*repmat(sideedge(1:end-1,1:end-1,1:45),[1 1 1 nt])...
    ))));
difSHopenMag=squeeze(nansum(nansum(nansum(...
    (abs(DifxEs(1:end-1,1:end-1,1:45,:)).*repmat(openW(1:end-1,1:end-1,1:45),[1 1 1 nt])+abs(DifxEs(2:end,1:end-1,1:45,:)).*repmat(openE(1:end-1,1:end-1,1:45),[1 1 1 nt])...
    +abs(DifyEs(1:end-1,1:end-1,1:45,:)).*repmat(openS(1:end-1,1:end-1,1:45),[1 1 1 nt])+abs(DifyEs(1:end-1,2:end,1:45,:)).*repmat(openN(1:end-1,1:end-1,1:45),[1 1 1 nt]))...
    .*repmat(sideedge(1:end-1,1:end-1,1:45),[1 1 1 nt])...
    ))));
difSZopen2=squeeze(nansum(nansum(nansum(...
    (DifxEs(1:end-1,1:end-1,1:45,:).*repmat(openW(1:end-1,1:end-1,1:45),[1 1 1 nt])-DifxEs(2:end,1:end-1,1:45,:).*repmat(openE(1:end-1,1:end-1,1:45),[1 1 1 nt])...
    +DifyEs(1:end-1,1:end-1,1:45,:).*repmat(openS(1:end-1,1:end-1,1:45),[1 1 1 nt])-DifyEs(1:end-1,2:end,1:45,:).*repmat(openN(1:end-1,1:end-1,1:45),[1 1 1 nt]))...
    .*repmat(bottomedge(1:end-1,1:end-1,1:45),[1 1 1 nt])...
    ))));
difSZopen2Mag=squeeze(nansum(nansum(nansum(...
    (abs(DifxEs(1:end-1,1:end-1,1:45,:)).*repmat(openW(1:end-1,1:end-1,1:45),[1 1 1 nt])+abs(DifxEs(2:end,1:end-1,1:45,:)).*repmat(openE(1:end-1,1:end-1,1:45),[1 1 1 nt])...
    +abs(DifyEs(1:end-1,1:end-1,1:45,:)).*repmat(openS(1:end-1,1:end-1,1:45),[1 1 1 nt])+abs(DifyEs(1:end-1,2:end,1:45,:)).*repmat(openN(1:end-1,1:end-1,1:45),[1 1 1 nt]))...
    .*repmat(bottomedge(1:end-1,1:end-1,1:45),[1 1 1 nt])...
    ))));
difSZopen=difSZopen1+difSZopen2;
% difSaltH=squeeze(nansum(nansum(nansum((DifxEs(1:end-1,1:end-1,:,:)-DifxEs(2:end,1:end-1,:,:)+DifyEs(1:end-1,1:end-1,:,:)...
%     -DifyEs(1:end-1,2:end,:,:)).*Rho(1:end-1,1:end-1,:,:).*repmat(inWag(1:end-1,1:end-1,:),[1 1 1 nt])))));
% difSaltHopenS=squeeze(nansum(nansum(nansum(...
%     (DifxEs(1:end-1,1:end-1,:,:).*repmat(openW(1:end-1,1:end-1,:),[1 1 1 nt])-DifxEs(2:end,1:end-1,:,:).*repmat(openE(1:end-1,1:end-1,:),[1 1 1 nt])...
%     +DifyEs(1:end-1,1:end-1,:,:).*repmat(openS(1:end-1,1:end-1,:),[1 1 1 nt])-DifyEs(1:end-1,2:end,:,:).*repmat(openN(1:end-1,1:end-1,:),[1 1 1 nt]))...
%     .*Rho(1:end-1,1:end-1,:,:).*repmat(sideedge(1:end-1,1:end-1,:),[1 1 1 nt])...
%     ))));
% difSaltBopen2=squeeze(nansum(nansum(nansum(...
%     (DifxEs(1:end-1,1:end-1,1:45,:).*repmat(openW(1:end-1,1:end-1,1:45),[1 1 1 nt])-DifxEs(2:end,1:end-1,1:45,:).*repmat(openE(1:end-1,1:end-1,1:45),[1 1 1 nt])...
%     +DifyEs(1:end-1,1:end-1,1:45,:).*repmat(openS(1:end-1,1:end-1,1:45),[1 1 1 nt])-DifyEs(1:end-1,2:end,1:45,:).*repmat(openN(1:end-1,1:end-1,1:45),[1 1 1 nt]))...
%     .*Rho(1:end-1,1:end-1,1:45,:).*repmat(bottomedge(1:end-1,1:end-1,1:45),[1 1 1 nt])...
%     ))));
clear Dif*
disp('s dif done')

globalArea=2.130281953410312e45;
load('fluxesVelS148day.mat','WS')
ssurfcor=squeeze(nansum(nansum(squeeze(WS(:,:,1,:)).*repmat(RAC,[1 1 nt]))))./globalArea;
ssurfcor2=squeeze(nansum(nansum(squeeze(WS(:,:,1,:)).*repmat(RAC.*inWag(:,:,1),[1 1 nt]))))./globalArea;
for i=1:nt
surfStend3(i)=nansum(nansum(inWag(:,:,1).*(-squeeze(WS(:,:,1,i))).*cellVol4(:,:,1,i)./(dZ(1,1,1).*hFacC(:,:,1)))); %do I need hfac here
surfStend3Mag(i)=nansum(nansum(inWag(:,:,1).*(squeeze(abs(WS(:,:,1,i)))).*cellVol4(:,:,1,i)./(dZ(1,1,1).*hFacC(:,:,1)))); %do I need hfac here
%surfSaltTend3(i)=nansum(nansum(inWag(:,:,1).*(-squeeze(WS(:,:,1,i))).*Rho(:,:,1,i).*cellVol4(:,:,1,i)./(dZ(1,1,1).*hFacC(:,:,1))));
end
load('fluxesVelS148day.mat','US','VS')
% advSHopenS=squeeze(nansum(nansum(nansum(...
%     (US(1:end-1,1:end-1,:,:).*repmat(openW(1:end-1,1:end-1,:)./DXG(1:end-1,1:end-1,:),[1 1 1 nt])-US(2:end,1:end-1,:,:).*repmat(openE(1:end-1,1:end-1,:)./DXG(1:end-1,1:end-1,:),[1 1 1 nt])...
%     +VS(1:end-1,1:end-1,:,:).*repmat(openS(1:end-1,1:end-1,:)./DYG(1:end-1,1:end-1,:),[1 1 1 nt])-VS(1:end-1,2:end,:,:).*repmat(openN(1:end-1,1:end-1,:)./DYG(1:end-1,1:end-1,:),[1 1 1 nt]))...
%     .*repmat(hFacC(1:end-1,1:end-1,:).*sideedge(1:end-1,1:end-1,:),[1 1 1 nt])...
%     ))));
% advSBopen=squeeze(nansum(nansum(nansum(...
%     (WS(1:end-1,1:end-1,2:46,:).*repmat(openZd(1:end-1,1:end-1,1:45),[1 1 1 nt])-WS(1:end-1,1:end-1,1:45,:).*repmat(openZu(1:end-1,1:end-1,1:45),[1 1 1 nt])).*repmat(hFacC(1:end-1,1:end-1,1:45),[1 1 1 nt])./dZ4(1:end-1,1:end-1,1:45,:)...
%     +(US(1:end-1,1:end-1,1:45,:).*repmat(openW(1:end-1,1:end-1,1:45)./DXG(1:end-1,1:end-1,1:45),[1 1 1 nt])-US(2:end,1:end-1,1:45,:).*repmat(openE(1:end-1,1:end-1,1:45)./DXG(1:end-1,1:end-1,1:45),[1 1 1 nt])...
%     +VS(1:end-1,1:end-1,1:45,:).*repmat(openS(1:end-1,1:end-1,1:45)./DYG(1:end-1,1:end-1,1:45),[1 1 1 nt])-VS(1:end-1,2:end,1:45,:).*repmat(openN(1:end-1,1:end-1,1:45)./DYG(1:end-1,1:end-1,1:45),[1 1 1 nt]))...
%     .*repmat(hFacC(1:end-1,1:end-1,1:45).*bottomedge(1:end-1,1:end-1,1:45),[1 1 1 nt])...
%     ))));
% advSHopenWeighted=squeeze(nansum(nansum(nansum(...
%     (US(1:end-1,1:end-1,:,:).*repmat(openW(1:end-1,1:end-1,:)./DXG(1:end-1,1:end-1,:),[1 1 1 nt])-US(2:end,1:end-1,:,:).*repmat(openE(1:end-1,1:end-1,:)./DXG(1:end-1,1:end-1,:),[1 1 1 nt])...
%     +VS(1:end-1,1:end-1,:,:).*repmat(openS(1:end-1,1:end-1,:)./DYG(1:end-1,1:end-1,:),[1 1 1 nt])-VS(1:end-1,2:end,:,:).*repmat(openN(1:end-1,1:end-1,:)./DYG(1:end-1,1:end-1,:),[1 1 1 nt]))...
%     .*repmat(hFacC(1:end-1,1:end-1,:),[1 1 1 nt])...
%     ))));
% advSZopenWeighted=squeeze(nansum(nansum(nansum(...
%     (WS(1:end-1,1:end-1,2:46,:).*repmat(openZd(1:end-1,1:end-1,1:45),[1 1 1 nt])-WS(1:end-1,1:end-1,1:45,:).*repmat(openZu(1:end-1,1:end-1,1:45),[1 1 1 nt])).*repmat(hFacC(1:end-1,1:end-1,1:45),[1 1 1 nt])./dZ4(1:end-1,1:end-1,1:45,:)...
%     ))));
% usH=squeeze(nansum(nansum(nansum(...
%     ((US(1:end-1,1:end-1,:,:)-US(2:end,1:end-1,:,:))./repmat(DXG(1:end-1,1:end-1,:),[1 1 1 nt])...
%     +(VS(1:end-1,1:end-1,:,:)-VS(1:end-1,2:end,:,:))./repmat(DYG(1:end-1,1:end-1,:),[1 1 1 nt])).*repmat(inWag(1:end-1,1:end-1,:),[1 1 1 nt])...
%     ))));
% wsZ=squeeze(nansum(nansum(nansum(...
%     (WS(:,:,2:46,:)-WS(:,:,1:45,:)).*repmat(inWag(:,:,1:45),[1 1 1 nt])./dZ4(:,:,1:45,:)...
%     ))));
clear US VS WS
load('fluxesOtherS148day.mat')
surfF=squeeze(nansum(nansum(Sflux.*squeeze(cellVol4(:,:,1,:)).*repmat(inWag(:,:,1),[1 1 nt])./(dZ(1,1,1).*repmat(hFacC(:,:,1),[1 1 nt]).*Rho1))));%Sflux is g/(s*m^2), salinity flux g/(s*kg)
surfFMag=squeeze(nansum(nansum(abs(Sflux).*squeeze(cellVol4(:,:,1,:)).*repmat(inWag(:,:,1),[1 1 nt])./(dZ(1,1,1).*repmat(hFacC(:,:,1),[1 1 nt]).*Rho1))));%Sflux is g/(s*m^2), salinity flux g/(s*kg)
%dSdtNew=squeeze(nansum(nansum(nansum(Stend.*cellVol4.*repmat(inWag,[1 1 1 nt])/86400))));%./nansum(nansum(nansum(cellVol4.*repmat(inWag,[1 1 1 nt])/86400))));
dSdt=squeeze(nansum(nansum(nansum(Stend.*cellVol4.*repmat(inWag,[1 1 1 nt])/86400))));
dSdtMag=squeeze(nansum(nansum(nansum(abs(Stend).*cellVol4.*repmat(inWag,[1 1 1 nt])/86400))));
%surfFsalt=squeeze(nansum(nansum(Sflux.*squeeze(cellVol4(:,:,1,:)).*repmat(inWag(:,:,1),[1 1 nt])./(dZ(1,1,1).*repmat(hFacC(:,:,1),[1 1 nt])))));
%dSaltdt=squeeze(nansum(nansum(nansum(Stend.*cellVol4.*Rho.*repmat(inWag,[1 1 1 nt])/86400))));

surfFnet=squeeze(nansum(nansum(Fsw.*squeeze(cellVol4(:,:,1,:)).*repmat(inWag(:,:,1),[1 1 nt])./(dZ(1,1,1).*repmat(hFacC(:,:,1),[1 1 nt]).*Rho1))));
%load('tsTendsnap.mat','restartS')
%reS=squeeze(nansum(nansum(nansum(restartS.*cellVol4.*repmat(inWag,[1 1 1 nt])/86400))));

sumflux=advSh+advSz+difSh+difSz2+surfStend3.'+surfF;%+reS;
perrorS=100*(dSdt-sumflux)./dSdt;
%sumSaltFlux=advSaltH+advSaltZ+difSaltH+difSaltZ2+surfSaltTend3.'+surfFsalt;
%perrorSalt=100*(dSaltdt-sumSaltFlux)./dSaltdt;
disp('s surf done')
clear Sflux Stend Fsw

% %add back in when run finishes %doesnt have cellVol4 in yet
% load('divTS162NR.mat','ShDiv')
% sdwdz=-squeeze(nansum(nansum(nansum(ShDiv(:,:,:,1:nt).*repmat(inWag,[1 1 1 nt])))));
% %saltdwdz=-squeeze(nansum(nansum(nansum(ShDiv(:,:,:,1:nt).*repmat(inWag,[1 1 1 nt]).*Rho.*cellVol4))));
% clear ShDiv
%  advSz2=advSz-sdwdz;
%  advSh2=advSh+sdwdz;
%  % advSaltZ2=advSaltZ-saltdwdz;
%  %advSaltH2=advSaltH+saltdwdz;

surfStot=surfStend3.'+surfF;
%errorS=dSdt-(advSh2+advSz2+difSz2+difSh+surfStend3.'+surfF);


disp('not saving stends')
clear dZ4 dZ hFac* DX* dx* cellVol* Rho* open* *edge
%save('salinityBudgetEulerWag148daysIntegralNF.mat')%,'-v7.3')
%% plots
% %difSaltBopen=difSaltBopen1+difSaltBopen2;
% 
figure;
plot(advSHopenS,'linewidth',2); hold all 
plot(advSBopen,'linewidth',2)%-advSZu
%plot(advSZu,'linewidth',2)%,'advection through top'
plot(difSZopen,'linewidth',2); plot(surfStend3.'+surfF,'linewidth',2)
plot(-dSdt,'linewidth',2)
plot(difSHopen,'linewidth',2);
plot(advSHopenS+advSBopen+difSHopen+surfStend3.'+surfF-dSdt,'k--','linewidth',4)%+difSZopen
legend('advection through sides','advection through bottom','diffusion through bottom','surface incl SSH','-dS/dt','diffusion through sides','sum')
xlabel('days','fontsize',20)
ylabel('integrated changes in salinity, psu m^3/s','fontsize',20)
title('Euler WAG Salinity Budget','fontsize',22)
set(gca,'fontsize',20)
%%
rho0=1000; sbar=-36.61;
dvdtPerfect=(vfluxhs+vfluxb+precip);
dvdt=[volumeSSHflux;0];
vcorrection=-(vfluxhs+vfluxb+precip-[volumeSSHflux;0]);
figure;
plot(rho0*(advSHopenS+surfStend3.')+rho0*sbar*(vfluxhs+vcorrection),'linewidth',2); hold all 
plot(advSBopen*rho0+rho0*sbar*vfluxb,'linewidth',2)%-advSZu
%plot(advSZu,'linewidth',2)%,'advection through top'
plot(rho0*difSZopen,'linewidth',2); plot(rho0*(surfF+precip*sbar),'linewidth',2)
plot(-rho0*(dSdt+sbar*dvdt),'linewidth',2)
plot(rho0*difSHopen,'linewidth',2);
plot(rho0*(advSHopenS+advSBopen+difSHopen+surfStend3.'+surfF-dSdt)+rho0*sbar*(vfluxhs+vcorrection+vfluxb+precip-dvdt),'k--','linewidth',4)%+difSZopen
legend('advection through sides','advection through bottom','diffusion through bottom','surface incl SSH','-dS/dt','diffusion through sides','sum')
xlabel('days','fontsize',20)
ylabel('integrated changes in salt, g/s','fontsize',20)
title('Euler WAG Salt Budget with salinity 36.61 removed','fontsize',22)
set(gca,'fontsize',18)

%%
meanTermsC={'h adv','z adv','z diff','surface','-dSdt','total'};
meanTerms=[mean(rho0*(advSHopenS+surfStend3.')+rho0*sbar*(vfluxhs+vcorrection)),mean(advSBopen*rho0+rho0*sbar*vfluxb)...
    mean(rho0*difSZopen),mean(rho0*(surfF+precip*sbar)),mean(-rho0*(dSdt+sbar*dvdt)),...
    mean(rho0*(advSHopenS+advSBopen+difSHopen+surfStend3.'+surfF-dSdt)+rho0*sbar*(vfluxhs+vcorrection+vfluxb+precip-dvdt))/10];
figure; bar(meanTerms)
set(gca,'XTickLabel',meanTermsC)
title('Mean Euler WAG Salt Fluxes, reference salinity 36.61','fontsize',14)
set(gca,'fontsize',14)
ylabel('g/s','fontsize',14)

% figure;
% plot(advSHopenSMag,'linewidth',2); hold all 
% plot(advSBopenMag,'linewidth',2)
% plot(difSZopen1Mag+difSZopen2Mag,'linewidth',2); plot(surfFMag,'linewidth',2)
% plot(dSdtMag,'linewidth',2)
% plot(surfStend3Mag.','linewidth',2)
% plot(difSHopenMag,'linewidth',2);
% legend('advection through sides','advection through bottom','diffusion through bottom','surface','-dS/dt','SSH','diffusion through sides')
% xlabel('days','fontsize',20)
% ylabel('magnitude of salinity fluxes, psu m^3/s','fontsize',20)
% set(gca,'YScale','log')
% title('Euler WAG Salinity Budget Flux Magnitudes','fontsize',22)
% set(gca,'fontsize',20)
% 
% % figure;
% % plot(advSHopenS+advSBopen,'linewidth',2); hold all 
% % plot(difSHopen,'linewidth',2);
% % plot(difSZopen,'linewidth',2); plot(surfStend3.'+surfF,'linewidth',2)
% % plot(-dSdt,'linewidth',2)
% % plot(advSHopenS+advSBopen+difSHopen+difSZopen+surfStend3.'+surfF-dSdt,'--','linewidth',4)
% % legend('advection','diffusion through sides','diffusion through bottom','surface','-dS/dt','sum')
% % xlabel('days','fontsize',20)
% % ylabel('integrated changes in salinity, psu m^3/s','fontsize',20)
% % title('Euler WAG Salinity Budget','fontsize',22)
% % set(gca,'fontsize',20)
%% 
figure;
plot(advSHopenS+advSBopen,'linewidth',2); hold all 
%plot(advSZu+surfStend3.','linewidth',2)
plot(nan,nan)
plot(difSHopen,'linewidth',2);
plot(difSZopen,'linewidth',2); %plot(surfF,'linewidth',2)
plot(-dSdt,'linewidth',2)
plot(advSHopenS+advSBopen+difSHopen+difSZopen+surfStend3.'+surfF-dSdt,'k--','linewidth',4)
legend('advection','h diffusion','z diffusion','-dS/dt','sum')
xlabel('days','fontsize',20)
ylabel('integrated changes in salinity, psu m^3/s','fontsize',20)
title('Euler WAG Salinity Budget','fontsize',22)
set(gca,'fontsize',20)
%% 
% figure;
% plot(advSHopenS+advSBopen+surfStend3.'-dSdt,'linewidth',2); hold all 
% plot(nan,nan)
% %plot(difSHopen,'linewidth',2);
% plot(difSZopen,'linewidth',2); plot(surfF,'linewidth',2)
% plot(advSHopenS+advSBopen+difSHopen+difSZopen+surfStend3.'+surfF-dSdt,'k--','linewidth',4)
% legend('advection including SSH-dSdt','diffusion through bottom','surface','sum')
% xlabel('days','fontsize',20)
% ylabel('integrated changes in salinity, psu m^3/s','fontsize',20)
% title('Euler WAG Salinity Budget','fontsize',22)
% set(gca,'fontsize',20)
% 
% 
% % figure;
% % plot(advSHopenS+advSBopen+surfStend3.'-dSdt+surfF,'linewidth',2); hold all 
% % plot(difSHopen,'linewidth',2);
% % plot(difSZopen,'linewidth',2); %plot(surfF,'linewidth',2)
% % plot(advSHopenS+advSBopen+difSHopen+difSZopen+surfStend3.'+surfF-dSdt,'--','linewidth',4)
% % legend('advection including SSH+surface-dSdt','diffusion through sides','diffusion through bottom','sum')
% % xlabel('days','fontsize',20)
% % ylabel('integrated changes in salinity, psu m^3/s','fontsize',20)
% % title('Euler WAG Salinity Budget','fontsize',22)
% % set(gca,'fontsize',20)
% 
% % figure;
% % plot(advSh+advSz,'linewidth',2); hold all 
% % plot(difSh,'linewidth',2);
% % plot(difSz2,'linewidth',2); plot(surfStend3.'+surfF,'linewidth',2)
% % plot(-dSdt,'linewidth',2)
% % plot(sumflux-dSdt,'--','linewidth',4)
% % legend('advection','horizontal diffusion','vertical diffusion','surface','-dS/dt','sum')
% % xlabel('days')
% % ylabel('changes in salinity, psu/s')
% % title('Euler salinity budget')
% % 
% % figure;
% % plot(advSaltH+advSaltZ,'linewidth',2); hold all 
% % plot(difSaltH,'linewidth',2);
% % plot(difSaltZ2,'linewidth',2); plot(surfSaltTend3.'+surfFsalt,'linewidth',2)
% % plot(-dSaltdt,'linewidth',2)
% % plot(sumSaltFlux-dSaltdt,'--','linewidth',4)
% % legend('advection','horizontal diffusion','vertical diffusion','surface','-dS/dt','sum')
% % xlabel('days')
% % ylabel('changes in salt content, g/s')
% % title('Euler salt budget')
% % 
% % difSaltBopen=difSaltBopen1+difSaltBopen2;
% % figure;
% % plot(advSaltHopenS,'linewidth',2); hold all 
% % plot(advSaltBopen,'linewidth',2)
% % plot(difSaltHopenS,'linewidth',2);
% % plot(difSaltBopen,'linewidth',2); plot(surfSaltTend3.'+surfFsalt,'linewidth',2)
% % plot(-dSaltdt,'linewidth',2)
% % plot(advSaltHopenS+advSaltBopen+difSaltHopenS+difSaltBopen+surfSaltTend3.'+surfFsalt-dSaltdt,'--','linewidth',4)
% % legend('advection through sides','advection through bottom','diffusion through sides','diffusion through bottom','surface','-dS/dt','sum')
% % xlabel('days')
% % ylabel('changes in salt content, g/s')
% % title('Euler salt budget')
% % 
% % figure;
% % plot(advSaltHopenS+advSaltBopen,'linewidth',2); hold all 
% % %plot(advSaltBopen,'linewidth',2)
% % plot(difSaltHopenS,'linewidth',2);
% % plot(difSaltBopen,'linewidth',2); plot(surfSaltTend3.'+surfFsalt,'linewidth',2)
% % plot(-dSaltdt,'linewidth',2)
% % plot(advSaltHopenS+advSaltBopen+difSaltHopenS+difSaltBopen+surfSaltTend3.'+surfFsalt-dSaltdt,'--','linewidth',4)
% % legend('advection through sides and bottom','diffusion through sides','diffusion through bottom','surface','-dS/dt','sum')
% % xlabel('days')
% % ylabel('changes in salt content, g/s')
% % title('Euler salt budget')
% % 
% % meanTermsC={'side adv/10','bottom adv/10','side diff','bottom diff','surface','-dSdt','total'};
% % meanTerms=[mean(advSaltHopenS)/10,mean(advSaltBopen)/10,mean(difSaltHopenS),mean(difSaltBopen),mean(surfSaltTend3.'+surfFsalt),mean(-dSaltdt),mean(advSaltHopenS+advSaltBopen+difSaltHopenS+difSaltBopen+surfSaltTend3.'+surfFsalt-dSaltdt)];
% % figure; bar(meanTerms)
% % set(gca,'XTickLabel',meanTermsC)
% % title('Mean Euler WAG Salt Fluxes')
% % 
% % meanTermsC={'side adv/10','bottom adv/10','side diff','bottom diff','-dSdt','total'};
% % meanTerms=[mean(advSaltHopenS)/10,mean(advSaltBopen)/10,mean(difSaltHopenS),mean(difSaltBopen),mean(-dSaltdt),mean(advSaltHopenS+advSaltBopen+difSaltHopenS+difSaltBopen-dSaltdt)];
% % figure; bar(meanTerms)
% % set(gca,'XTickLabel',meanTermsC)
% % title('Mean Euler WAG Salt Fluxes')
% % 
% % meanTermsC={'adv','side diff','bottom diff','surface','-dSdt','total'};
% % meanTerms=[mean(advSaltHopenS+advSaltBopen),mean(difSaltHopenS),mean(difSaltBopen),mean(surfSaltTend3.'+surfFsalt),mean(-dSaltdt),mean(advSaltHopenS+advSaltBopen+difSaltHopenS+difSaltBopen+surfSaltTend3.'+surfFsalt-dSaltdt)];
% % figure; bar(meanTerms)
% % set(gca,'XTickLabel',meanTermsC)
% % title('Mean Euler WAG Salt Fluxes')
% % 
% % meanTermsC={'adv','side diff','bottom diff','surface','-dSdt','total'};
% % meanTerms=[mean(advSaltH+advSaltZ),mean(difSaltHopenS),mean(difSaltBopen),mean(surfSaltTend3.'+surfFsalt),mean(-dSaltdt),mean(advSaltH+advSaltZ+difSaltHopenS+difSaltBopen+surfSaltTend3.'+surfFsalt-dSaltdt)];
% % figure; bar(meanTerms)
% % set(gca,'XTickLabel',meanTermsC)
% % title('Mean Euler WAG Salt Fluxes')
% % 
meanTermsC={'adv','z diff','surface','-dSdt','total'};
meanTerms=[mean(advSHopenS+advSBopen+surfStend3.'),mean(difSZopen),mean(surfF),mean(-dSdt),mean(advSHopenS+advSBopen+difSHopen+difSZopen+surfStend3.'+surfF-dSdt)];
figure; bar(rho0*meanTerms)
set(gca,'XTickLabel',meanTermsC)
title('Mean Euler WAG Salt Transports','fontsize',14)
set(gca,'fontsize',14)
ylabel('g/s','fontsize',14)

%% 
meanTermsC={'adv top','adv bot','adv side','z diff','-dSdt','total'};
meanTerms=[mean(advSZu),mean(advSBopen),mean(advSHopen),mean(difSZopen),mean(-dSdt),mean(advSHopenS+advSBopen+difSHopen+difSZopen+surfStend3.'+surfF-dSdt)];
figure; bar(rho0*meanTerms)
set(gca,'XTickLabel',meanTermsC)
title('Mean Euler WAG Salt Transports','fontsize',14)
set(gca,'fontsize',14)
ylabel('g/s','fontsize',14)
meanTermsC={'adv','z diff','-dSdt','total'};
meanTerms=[mean(advSHopenS+advSBopen+surfStend3.'),mean(difSZopen),mean(-dSdt),mean(advSHopenS+advSBopen+difSHopen+difSZopen+surfStend3.'+surfF-dSdt)];
figure; bar(rho0*meanTerms)
set(gca,'XTickLabel',meanTermsC)
title('Mean Euler WAG Salt Transports','fontsize',14)
set(gca,'fontsize',14)
ylabel('g/s','fontsize',14)
% 
% % meanTermsC={'h adv','z adv','h diff','z diff','surface','-dSdt','total'};
% % meanTerms=[mean(advSHopenS),mean(advSBopen),mean(difSHopen),mean(difSZopen),mean(surfStend3.'+surfF),mean(-dSdt),mean(advSHopenS+advSBopen+difSHopen+difSZopen+surfStend3.'+surfF-dSdt)];
% % figure; bar(meanTerms)
% % set(gca,'XTickLabel',meanTermsC)
% % title('Mean Euler WAG Salinity Fluxes','fontsize',14)
% % ylabel('psu m^3/s','fontsize',14)
% % set(gca,'fontsize',14)
% 
% meanTermsC={'adv','SSH','h diff','z diff','surface','-dSdt'};
% meanTerms=[mean(advSHopenSMag+advSBopenMag),mean(surfStend3Mag.'),mean(difSHopenMag),mean(difSZopen1Mag+difSZopen2Mag),mean(surfFMag),mean(dSdtMag)];
% figure; bar(meanTerms)
% set(gca,'XTickLabel',meanTermsC)
% title('Mean Euler WAG Salinity Flux Magnitudes','fontsize',14)
% ylabel('psu m^3/s','fontsize',14)
% set(gca,'fontsize',14)
% set(gca,'yscale','log')
% 
% % meanTermsC={'adv','h diff','z diff','surface','dSdt','total'};
% % meanTerms=[mean(abs(advSHopenS+advSBopen)),mean(abs(difSHopen)),mean(abs(difSZopen)),mean(abs(surfStend3).'+abs(surfF)),mean(abs(dSdt)),mean(abs(advSHopenS+advSBopen+difSHopen+difSZopen+surfStend3.'+surfF-dSdt))];
% % figure; bar(meanTerms)
% % set(gca,'XTickLabel',meanTermsC)
% % title('Mean Euler WAG Salinity Flux Magnitudes')
%% old junk
% advSh=squeeze(nansum(nansum(nansum( ((US(1:end-1,1:end-1,:,:)-US(2:end,1:end-1,:,:))./repmat(DXC(1:end-1,1:end-1),[1 1 46 nt])...
%     +(VS(1:end-1,1:end-1,:,:)-VS(1:end-1,2:end,:,:))./repmat(DYC(1:end-1,1:end-1),[1 1 46 nt])).*repmat(inWag(1:end-1,1:end-1,:),[1 1 1 nt])./repmat(hFacC(1:end-1,1:end-1,:),[1 1 1 nt])))));
% advSz=squeeze(nansum(nansum(nansum((WS(:,:,2:46,:)-WS(:,:,1:45,:)).*repmat(inWag(:,:,1:45),[1 1 1 nt])./repmat(dZ(:,:,1:45),[700 200 1 nt])))));
% 
% difSz3=difSz2+numZ;
% difSh3=difSh+numH;
%short1(1,1,1:46)=0.62*exp(dInterface(1:46)/0.6)+0.38*exp(dInterface(1:46)/20);%no shortwave!
%short2(1,1,1:46)=0.62*exp(dInterface(2:47)/0.6)+0.38*exp(dInterface(2:47)/20);
%shortwaveApprox=Qsw/(1026*4000)/(repmat(dZ,[700 200 1 45]).*repmat(hFacC,[1 1 1 45]))*repmat(short1-short2,[700 200 1 45]);
%surfStend(i)=nansum(nansum(inWag(:,:,1).*(ssurfcor(i).*ones([700 200])-squeeze(WS(:,:,1,i)))./(dZ(1,1,1).*hFacC(:,:,1))));
%surfStend2(i)=nansum(nansum(inWag(:,:,1).*(ssurfcor2(i).*ones([700 200])-squeeze(WS(:,:,1,i)))./(dZ(1,1,1).*hFacC(:,:,1))));
%%surfQapprox=squeeze(nansum(nansum(Tflux.*repmat(inWag(:,:,1),[1 1 45])./(3994*dZ(1,1,1).*repmat(hFacC(:,:,1),[1 1 45]).*Rho))));%approx bc using const Cp
