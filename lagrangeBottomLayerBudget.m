%sigma 27 to 27.5 budget for volume
load('sigma27WAGareas.mat', 'inwagSareas')
inwag27areas=inwagSareas;
load('sigma2675WAGareas.mat', 'inwagSareas')
inwag2675areas=inwagSareas;
fillneeded=find(squeeze(nansum(nansum(inwag27areas)))==0);
fillneeded=fillneeded(14:end);
fillneeded27=fillneeded;
fillneeded2=find(squeeze(nansum(nansum(inwag2675areas)))==0);
fillneeded2675=fillneeded2(14:end);
load('sigma275WAGareas.mat', 'inwagSareas')
inwag275areas=inwagSareas;
fillneeded2=find(squeeze(nansum(nansum(inwag275areas)))==0);
fillneeded275=fillneeded2(14:end);
load('sigma265WAGareas.mat', 'inwagSareas')
inwag265areas=inwagSareas;
fillneeded2=find(squeeze(nansum(nansum(inwag265areas)))==0);
fillneeded265=fillneeded2(14:end);
load('sigma263WAGareas.mat', 'inwagSareas')
inwag263areas=inwagSareas;
fillneeded2=find(squeeze(nansum(nansum(inwag263areas)))==0);
fillneeded263=fillneeded2(14:end);
load('surfaceWAGareas.mat', 'inwagSareas')
fillneeded2=find(squeeze(nansum(nansum(inwagSareas)))==0);
fillneededS=fillneeded2(14:end);

allfill=unique([fillneededS;fillneeded263;fillneeded265;fillneeded2675;fillneeded27]);

load('gateAdvectionSigma27indexedMar.mat', 'gateFlux263265')
gateFlux27275=gateFlux263265;
load('sigma27EdgeCheck3NEWnocoast.mat', 'volNoGate2')
load('diffusiveMovedIsopycnal27.mat','isoDepth')
iso27moved=isoDepth;
load('diffusiveMovedIsopycnal275.mat','isoDepth')
iso275moved=isoDepth;
load('iso275depthNFrev','isoDepth')
iso275=isoDepth;
load('iso27depthNFrev','isoDepth')
iso27=isoDepth;

load('diffusiveMovedIsopycnal2675.mat','isoDepth')
iso2675moved=isoDepth;
load('diffusiveMovedIsopycnal265.mat','isoDepth')
iso265moved=isoDepth;
load('diffusiveMovedIsopycnal263.mat','isoDepth')
iso263moved=isoDepth;

load('iso2675depthNFrev','isoDepth')
iso2675=isoDepth;
load('iso265depthNFrev','isoDepth')
iso265=isoDepth;
load('iso263depthNFrev','isoDepth')
iso263=isoDepth;
load('uvwSSHDailyDepth1rotated148F.mat', 'SSHa')

nt=146;
load('geometrySpinupSteady','dInterface')
%% filling in layers in order!
inwag263areas(:,:,fillneeded263)=inwagSareas(:,:,fillneeded263);
inwag265areas(:,:,fillneeded265)=inwag263areas(:,:,fillneeded265);
inwag2675areas(:,:,fillneeded2675)=inwag265areas(:,:,fillneeded2675);
inwag27areas(:,:,fillneeded27)=inwag2675areas(:,:,fillneeded27);
inwag275areas(:,:,fillneeded275)=inwag27areas(:,:,fillneeded275);
%% isopycnal motion volume flux diffusion
 volDiff275in27=inwag27areas.*(iso275moved(:,:,1:nt)-iso275(:,:,1:nt));%positive is deepening, so volume change is into gyre
 dvdtDiff275in27=squeeze(nansum(nansum(volDiff275in27)))/86400;
 volDiff27in27=inwag27areas.*(iso27moved(:,:,1:nt)-iso27(:,:,1:nt));%positive is deepening, so volume change is up into gyre
 dvdtDiff27in27=squeeze(nansum(nansum(volDiff27in27)))/86400;
 %so dvdtDiff275in27 positive is into the layer, want that -dvdtDiff27in27 (b/c positive there is out the top)
volDiff275in275=inwag275areas.*(iso275moved(:,:,1:nt)-iso275(:,:,1:nt));%positive is deepening, so volume change is into gyre
 dvdtDiff275in275=squeeze(nansum(nansum(volDiff275in275)))/86400;
 
  volDiff27in2675=inwag2675areas.*(iso27moved(:,:,1:nt)-iso27(:,:,1:nt));%positive is deepening, so volume change is into gyre
 dvdtDiff27in2675=squeeze(nansum(nansum(volDiff27in2675)))/86400;
   volDiff2675in2675=inwag2675areas.*(iso2675moved(:,:,1:nt)-iso2675(:,:,1:nt));%positive is deepening, so volume change is into gyre
 dvdtDiff2675in2675=squeeze(nansum(nansum(volDiff2675in2675)))/86400;
  volDiff2675in265=inwag265areas.*(iso2675moved(:,:,1:nt)-iso2675(:,:,1:nt));%positive is deepening, so volume change is into gyre
 dvdtDiff2675in265=squeeze(nansum(nansum(volDiff2675in265)))/86400;
   volDiff265in265=inwag265areas.*(iso265moved(:,:,1:nt)-iso265(:,:,1:nt));%positive is deepening, so volume change is into gyre
 dvdtDiff265in265=squeeze(nansum(nansum(volDiff265in265)))/86400;
 
 %% volume in the layer
 load('saTCtCpSigma162NF.mat', 'Sigma')
Sigma=Sigma(:,:,1:30,1:146);
dZ=diff(dInterface);
for k=1:30%through 16 is enough!
    k
inWag3area(:,:,k,1:146)=double(squeeze(Sigma(:,:,k,1:146))<27.5).*double(squeeze(Sigma(:,:,k,1:146))>=27).*inwag27areas(:,:,1:146).*dZ(k);
end
vol27275=squeeze(nansum(nansum(nansum(inWag3area))));
vol27275(vol27275==0)=nan;
dvdt1=diff(vol27275)/86400;
%
vol27275b=squeeze(nansum(nansum((inwag27areas.*(iso275(:,:,1:146)-iso27(:,:,1:146))))));
vol27275b(vol27275b==0)=nan;
dvdt2=diff(vol27275b)/86400;


vol27b=squeeze(nansum(nansum((inwag27areas.*(0.5*(iso275(:,:,1:146)+iso27(:,:,1:146))-0.5*(iso27(:,:,1:146)+iso2675(:,:,1:146)))))));
vol275b=squeeze(nansum(nansum((inwag275areas.*(iso275(:,:,1:146)-0.5*(iso27(:,:,1:146)+iso275(:,:,1:146)))))));
vol27b275b=vol27b+vol275b;
dvdt2centered=diff(vol27b275b)/86400;
dvdt27=diff(vol27)/86400; dvdt27(13)=0;
%% volume from snapshot sigma
%  load('tsSigmaSnapshots162NF.mat', 'Sigma')
%  Sigma=Sigma(:,:,1:30,1:148);
% load('iso275depthNFsnapRev','isoDepth')
% iso275s=isoDepth;
% load('iso27depthNFsnapRev','isoDepth')
% iso27s=isoDepth;

dZ=diff(dInterface);
for k=1:30%through 16 is enough!
    k
inWag3area(:,:,k,1:146)=double(squeeze(Sigma(:,:,k,1:146))<27.5).*double(squeeze(Sigma(:,:,k,1:146))>=27).*inwag27areas(:,:,1:146).*dZ(k);
end
vol27275c=squeeze(nansum(nansum(nansum(inWag3area))));
%vol27275c(vol27275c==0)=nan;
dvdt1s=diff(vol27275c)/86400;
%
vol27275d=squeeze(nansum(nansum((inwag27areas.*(iso275s(:,:,1:146)-iso27s(:,:,1:146))))));
%vol27275d(vol27275d==0)=nan;
dvdt2s=diff(vol27275d)/86400;

figure; plot(vol27275); hold all; plot(vol27275b); plot(vol27275c); plot(vol27275d); legend('Sigma ave','Isodepth ave','Sigma snap','isodepth snap')
figure; plot(dvdt1); hold all; plot(dvdt2); plot(dvdt1s); plot(dvdt2s); legend('Sigma ave','Isodepth ave','Sigma snap','isodepth snap')

%% errors from cutting manifolds
inwag27areasNew=inwag27areas;
 load('lagrangeWagAreasFilledOld.mat', 'inwag27areas')
 vol27275uncut=squeeze(nansum(nansum((inwag27areas.*(iso275(:,:,1:146)-iso27(:,:,1:146))))));
 dvdt2uncut=diff(vol27275uncut)/86400;
 largecuts=(abs(-dvdt2+dvdt2uncut)>2e5);
 %% more layers
 %old, uncut
vol267527o=squeeze(nansum(nansum((inwag2675areas.*(iso27(:,:,1:146)-iso2675(:,:,1:146))))));
vol2652675o=squeeze(nansum(nansum((inwag265areas.*(iso2675(:,:,1:146)-iso265(:,:,1:146))))));
vol263265o=squeeze(nansum(nansum((inwag263areas.*(iso265(:,:,1:146)-iso263(:,:,1:146))))));
volS263o=squeeze(nansum(nansum((inwagSareas.*(SSHa(:,:,1:146)+iso263(:,:,1:146))))));
 %new, cut
vol267527n=squeeze(nansum(nansum((inwag2675areasNew.*(iso27(:,:,1:146)-iso2675(:,:,1:146))))));
vol265267n=squeeze(nansum(nansum((inwag265areasNew.*(iso2675(:,:,1:146)-iso265(:,:,1:146))))));
vol263265n=squeeze(nansum(nansum((inwag263areasNew.*(iso265(:,:,1:146)-iso263(:,:,1:146))))));
volS263n=squeeze(nansum(nansum((inwagSareasNew.*(SSHa(:,:,1:146)+iso263(:,:,1:146))))));
%snap
vol2652675d=squeeze(nansum(nansum((inwag265areasNew.*(iso2675s(:,:,1:146)-iso265s(:,:,1:146))))));
vol267527d=squeeze(nansum(nansum((inwag265areasNew.*(iso27s(:,:,1:146)-iso2675s(:,:,1:146))))));
 %% total volume
vol27275b=squeeze(nansum(nansum((inwag27areas.*(iso275(:,:,1:146)-iso27(:,:,1:146))))));
vol267527b=squeeze(nansum(nansum((inwag2675areas.*(iso27(:,:,1:146)-iso2675(:,:,1:146))))));
vol2652675b=squeeze(nansum(nansum((inwag265areas.*(iso2675(:,:,1:146)-iso265(:,:,1:146))))));
vol263265b=squeeze(nansum(nansum((inwag263areas.*(iso265(:,:,1:146)-iso263(:,:,1:146))))));
volS263b=squeeze(nansum(nansum((inwagSareas.*(SSHa(:,:,1:146)+iso263(:,:,1:146))))));

volTotMar=volS263b+vol263265b+vol2652675b+vol267527b+vol27275b;
 dvdtMar=diff(volTotMar)./86400; dvdtMar(13)=0;
 %% 
 %11 days large cuts, 9 days filled
 notclosed=-dvdt2(14:145)+gateFlux27275(15:146).'+dvdtDiff275in27(14:145)-dvdtDiff27in27(14:145)+volNoGate2(14:145);
 gooddays(14:145)=abs(notclosed)<1e5;
 sum(gooddays(14:145))
 %% plot
 figure; plot(-dvdt2,'linewidth',2); hold all; plot(gateFlux27275(2:end),'linewidth',2); plot(dvdtDiff275in27,'linewidth',2)
 plot(-dvdtDiff27in27,'linewidth',2); plot(volNoGate2,'linewidth',2); 
 plot(14:145,-dvdt2(14:145)+gateFlux27275(15:146).'+dvdtDiff275in27(14:145)-dvdtDiff27in27(14:145)+volNoGate2(14:145),'--','linewidth',2)
 legend('-dV/dt via isodepth','gate','diffusion bottom','diffusion top','cross-manifold','total','known problem days')
 
%   figure; plot(-dvdt1,'linewidth',2); hold all; plot(gateFlux27275(2:end),'linewidth',2); plot(dvdtDiff275in27,'linewidth',2)
%  plot(-dvdtDiff27in27,'linewidth',2); plot(volNoGate2,'linewidth',2); 
%  plot(14:145,-dvdt1(14:145)+gateFlux27275(15:146).'+dvdtDiff275in27(14:145)-dvdtDiff27in27(14:145)+volNoGate2(14:145),'--','linewidth',2)
%  legend('-dV/dt via whole cell depth','gate','diffusion bottom','diffusion top','cross-manidold','total')
%%
figure; plot(-dvdt27,'linewidth',2); hold all; plot(gateFlux27(2:end),'linewidth',2); plot(dvdtDiff275in27,'linewidth',2)
 plot(-dvdtDiff27in27,'linewidth',2); plot(volNoGate2,'linewidth',2); 
 plot(14:145,-dvdt27(14:145)+gateFlux27(15:146).'+dvdtDiff275in27(14:145)-dvdtDiff27in27(14:145)+volNoGate2(14:145),'--','linewidth',2)
 plot(baddays27,3e5*ones(size(baddays27)),'r*')
 ylim([-4e5 4e5])
 legend('-dV/dt via isodepth','gate','diffusion bottom','diffusion top','cross-manifold','total','known problem days')
 %%
 figure;  plot(14:145,-dvdt2(14:145)+gateFlux27275(15:146).'+dvdtDiff275in27(14:145)-dvdtDiff27in27(14:145)+volNoGate2(14:145),'--','linewidth',2)
hold all;  plot(14:145,-dvdt27(14:145)+gateFlux27(15:146).'+gateFlux275(15:146).'+dvdtDiff275in275(14:145)-dvdtDiff27in27(14:145)+volNoGate2(14:145),'--','linewidth',2)

 title('Lower layer budget, centered in vertical')
 %% good sections only
 section1=17:23; section2=30:36; section3=40:50; section4=88:98; section5=134:142; section3issues=[43 48];
 gf=gateFlux267527(2:end);
 total1=-dvdt2675(1:145)+gf(2:146).'+dvdtDiff27in2675(1:145)-dvdtDiff2675in2675(1:145);%+volNoGate2(1:145);
 dvdt2=dvdt2675;
figure; plot(section1,-dvdt2(section1),'b','linewidth',2); hold all; plot(section1,gf(section1),'r','linewidth',2); plot(section1,dvdtDiff275in27(section1),'m','linewidth',2)
 plot(section1,-dvdtDiff27in27(section1),'c','linewidth',2)
 plot(section1,total1(section1),'k--','linewidth',2); plot(section1,volNoGate2(section1),'linewidth',2,'Color',[0.75 0 0.75]); 
 legend('-dV/dt','gate','diffusion bottom','diffusion top','total','cross-manifold error')
plot(section2,-dvdt2(section2),'b','linewidth',2); hold all; plot(section2,gf(section2),'r','linewidth',2); plot(section2,dvdtDiff275in27(section2),'m','linewidth',2)
 plot(section2,-dvdtDiff27in27(section2),'c','linewidth',2); plot(section2,volNoGate2(section2),'linewidth',2,'Color',[0.75 0 0.75]); 
 plot(section2,total1(section2),'k--','linewidth',2)
 plot(section3,-dvdt2(section3),'b','linewidth',2); hold all; plot(section3,gf(section3),'r','linewidth',2); plot(section3,dvdtDiff275in27(section3),'m','linewidth',2)
 plot(section3,-dvdtDiff27in27(section3),'c','linewidth',2); plot(section3,volNoGate2(section3),'linewidth',2,'Color',[0.75 0 0.75]); 
 plot(section3,total1(section3),'k--','linewidth',2)
 plot(section4,-dvdt2(section4),'b','linewidth',2); hold all; plot(section4,gf(section4),'r','linewidth',2); plot(section4,dvdtDiff275in27(section4),'m','linewidth',2)
 plot(section4,-dvdtDiff27in27(section4),'c','linewidth',2); plot(section4,volNoGate2(section4),'linewidth',2,'Color',[0.75 0 0.75]); 
 plot(section4,total1(section4),'k--','linewidth',2)
 plot(section5,-dvdt2(section5),'b','linewidth',2); hold all; plot(section5,gf(section5),'r','linewidth',2); plot(section5,dvdtDiff275in27(section5),'m','linewidth',2)
 plot(section5,-dvdtDiff27in27(section5),'c','linewidth',2); plot(section5,volNoGate2(section5),'linewidth',2,'Color',[0.75 0 0.75]); 
 plot(section5,total1(section5),'k--','linewidth',2)