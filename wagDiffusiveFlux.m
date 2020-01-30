%% explanation
%horizontal diffusion set to be run for each isopycnal layer
%vertical diffusion for sigma 27.5 to 28

%% location
%load('wagAreaAndFluxGate.mat', 'inWAG')%adjust by adding 26, 265, or 27 to filename
load('wagAreaAndFluxGate275.mat', '*Closed')
load('uvwAve148.mat', 'XC')
load('uvwAve148.mat', 'YC')
load('distancesAreas.mat','RAC')
figure;
for i=1:131
    c=[];
    if latWagClosed(i,1)>0
        np=find(latWagClosed(i,:)>0,1,'last');
        inWAG=inpolygon(XC,YC,lonWagClosed(i,1:np),latWagClosed(i,1:np));
        [c,h]=contour(XC,YC,double(inWAG),[1 1]);
    end
    if ~isempty(c)
    npointsi=find(c(2,:)==floor(c(2,:)));
    varhold=max(c(2,npointsi));
    i1=1+find(c(2,:)==varhold);
    iend=varhold+i1-1;
    inWAG2(:,:,i)=inpolygon(XC,YC,c(1,i1:iend),c(2,i1:iend));
    end
end

areaWAG2=squeeze(nansum(nansum(inWAG2.*repmat(RAC,[1 1 131]))));

clear inWAG
%inWAG2=sum(inWAG,3)>60;

%% horizontal diffusion, loading
kh=1e-9;%m^2/s numerical applied
kh2=100;%reasonable estimate maybe for ocean?? okubo suggests 10m^2/s for lengthscales of 10km

load('distancesAreas.mat','DXG')
load('distancesAreas.mat','DYG')
load('isopycnalDepths2.mat')
%upper surface
load('tsrhocpNativeGridIsoDepth27.mat')
salS=sIso;
tS=tIso;
cpS=cpIso;
rhoS=rhoIso;

% load('varying148ts16levelsRho.mat','S')
% salS=squeeze(S(:,:,1,:));
% clear S
% load('varying148ts16levelsRho.mat','Rho')
% rhoS=squeeze(Rho(:,:,1,:));
% clear Rho
% load('varying148ts16levelsRho.mat','T')
% tS=squeeze(T(:,:,1,:));
% clear T
% load('cp16levels148.mat');
% cpS=squeeze(cp(:,:,1,:)); 
% clear cp

%lower surface
load('tsrhocpNativeGridIsoDepth275.mat')

layerdepth=iso275depth-iso27depth;
%% horizontal diffusion: salt, heat
saltDiffHFlux=zeros([148,1]);
heatDiffHFlux=zeros([148,1]);

nh=1;
for k=1:148
    k
    csal=salS(:,:,k).*rhoS(:,:,k);
    ch=tS(:,:,k).*cpS(:,:,k).*rhoS(:,:,k);
    dZ=layerdepth(:,:,k);
    inWAG=inWAG2(:,:,k);
for i=200:550
    for j=30:199
        if inWAG(i,j)
            if ~inWAG(i-1,j) &&inWAG(i+1,j) &&csal(i-nh,j)>0 &&csal(i+nh-1,j)>0%left open
                %holdvar1(i,j,k)=-kh*(csal(i+10,j)-csal(i-11,j))*dZ(i,j)*DYG(i,j)/sum(DXG(i-11:i+10,j));
                saltDiffHFlux(k)=saltDiffHFlux(k)-kh*(csal(i+nh-1,j)-csal(i-nh,j))*dZ(i,j)*DYG(i,j)/sum(DXG(i-nh:i+nh-1,j));
                heatDiffHFlux(k)=heatDiffHFlux(k)-kh*(ch(i+nh-1,j)-ch(i-nh,j))*dZ(i,j)*DYG(i,j)/sum(DXG(i-nh:i+nh-1,j));
            end
            
            if ~inWAG(i+1,j) &&inWAG(i-1,j) &&csal(i+nh,j)>0 &&csal(i-nh+1,j)>0%right open
                saltDiffHFlux(k)=saltDiffHFlux(k)-kh*(csal(i-nh+1,j)-csal(i+nh,j))*dZ(i,j)*DYG(i,j)/sum(DXG(i-nh+1:i+nh,j));
                heatDiffHFlux(k)=heatDiffHFlux(k)-kh*(ch(i-nh+1,j)-ch(i+nh,j))*dZ(i,j)*DYG(i,j)/sum(DXG(i-nh+1:i+nh,j));

            end
            
            if ~inWAG(i,j+1) &&inWAG(i,j-1) &&csal(i,j-nh+1)>0 &&csal(i,j+nh)>0 %top open
                saltDiffHFlux(k)=saltDiffHFlux(k)-kh*(csal(i,j-nh+1)-csal(i,j+nh))*dZ(i,j)*DXG(i,j)/sum(DYG(i-nh+1:i+nh,j));
                heatDiffHFlux(k)=heatDiffHFlux(k)-kh*(ch(i,j-nh+1)-ch(i,j+nh))*dZ(i,j)*DXG(i,j)/sum(DYG(i-nh+1:i+nh,j));
            end
            
        end
    end
    %saltDiffHFlux(k)
end
%saltDiffHFlux(k)
end

%% horizontal diffusion: salt, surface to sigma=26, plot
figure; plot(saltDiffHFlux*kh2/kh); title('k=100m^2/s, salt flux (g/s)')
figure; plot(heatDiffHFlux*kh2/kh); title('k=100m^2/s, heat flux (W)')
%% save horizontal
fn='diffFluxHoriz27275.mat';
save(fn,'inWAG2','saltDiffHFlux','heatDiffHFlux','inWAG2','k*')
%% plot horizontal
load('diffFluxHoriz27275.mat')
saltDH27275=saltDiffHFlux;
heatDH27275=heatDiffHFlux;
load('diffFluxHoriz26527.mat')
saltDH26527=saltDiffHFlux;
heatDH26527=heatDiffHFlux;
load('diffFluxHoriz26265.mat')
saltDH26265=saltDiffHFlux;
heatDH26265=heatDiffHFlux;
load('diffFluxHorizS26.mat')
saltDHS26=saltDiffHFlux;
heatDHS26=heatDiffHFlux;

saltDHtot=saltDHS26+saltDH26265+saltDH26527+saltDH27275;
heatDHtot=heatDHS26+heatDH26265+heatDH26527+heatDH27275;

kh2
kh3=1e4

figure; plot(saltDHS26*kh3/kh); hold all; plot(saltDH26265*kh3/kh);
plot(saltDH26527*kh3/kh); plot(saltDH27275*kh3/kh); 
plot((saltDHS26(1:131)+saltDH26265(1:131)+saltDH26527(1:131)+saltDH27275(1:131))*kh3/kh,'k','LineWidth',2)
title('Diffusive Salt Flux, horizontal, k=10^4 m^2/s')
legend('surface','\sigma=26','\sigma=26.5','\sigma=27','total')
xlabel('time (days)')
ylabel('salt flux, g/s')

figure; plot(heatDHS26*kh3/kh); hold all; plot(heatDH26265*kh3/kh);
plot(heatDH26527*kh3/kh); plot(heatDH27275*kh3/kh); 
plot((heatDHS26(1:131)+heatDH26265(1:131)+heatDH26527(1:131)+heatDH27275(1:131))*kh3/kh,'k','LineWidth',2)
title('Diffusive Heat Flux, horizontal, k=10m^2/s')
legend('surface','\sigma=26','\sigma=26.5','\sigma=27','total')
xlabel('time (days)')
ylabel('heat flux, W')
%%

%% vertical salt&heat flux, load
load('tsrhocpNativeGridIsoDepth28.mat')
 s28=sIso; r28=rhoIso; t28=tIso; cp28=cpIso;
load('tsrhocpNativeGridIsoDepth275.mat')
 s275=sIso; r275=rhoIso; t275=tIso; cp275=cpIso;
 load('isopycnalDepths2.mat')
 layerDepth=iso28depth-iso275depth;
 %using location of sigma=27 WAG
 load('diffFluxHoriz27275.mat','inWAG2')
 
%% vertical salt&heat flux, calculate
kvI=1e-5;%imposed
kv=0.02;%based on typical values from kpp vertical diffusivity
for i=1:131
   csalgrad=(s275(:,:,i).*r275(:,:,i)-s28(:,:,i).*r28(:,:,i))./(layerDepth(:,:,i));
   saltDiffVFlux(i)=-kv*nansum(csalgrad(inWAG2(:,:,i)).*RAC(inWAG2(:,:,i))); %g salt/s
   cHgrad=(t275(:,:,i).*cp275(:,:,i).*r275(:,:,i)-t28(:,:,i).*cp28(:,:,i).*r28(:,:,i))./(layerDepth(:,:,i));
   heatDiffVFlux(i)=-kv*nansum(cHgrad(inWAG2(:,:,i)).*RAC(inWAG2(:,:,i))); %g salt/s
end
figure; plot(saltDiffVFlux)
%% plot
kh4=1e4;
gateSfluxT=(gateSflux26265+gateSflux26527+gateSflux27275+gateSfluxS26(1:138));
figure; plot(saltDiffVFlux); hold all; plot(saltDHtot*kh4/kh); plot(gateSfluxT) 
plot((gateSfluxT(1:131)+saltDiffVFlux(1:131)+(kh4/kh)*saltDHtot(1:131).'),'k','LineWidth',2) 
legend('vertical diffusion','horizontal diffusion','advection')

kh3=1e4;
ds1=saltDHS26*kh3/kh; ds1n=-ds1; ds1(ds1<0)=nan; ds1n(ds1n<0)=nan;
ds2=saltDH26265*kh3/kh; ds2n=-ds2; ds2(ds2<0)=nan; ds2n(ds2n<0)=nan;
ds3=saltDH26527*kh3/kh; ds3n=-ds3; ds3(ds3<0)=nan; ds3n(ds3n<0)=nan;
ds4=saltDH27275*kh3/kh; ds4n=-ds4; ds4(ds4<0)=nan; ds4n(ds4n<0)=nan;
figure; plot(ds1,'b'); hold on; 
plot(ds2,'Color',[0 0.5 0]);
plot(ds3,'r'); plot(ds4,'c'); 
plot(saltDHtot*kh3/kh,'LineWidth',2,'Color',[0.75 0 0.75])
plot(saltDiffVFlux,'LineWidth',2,'Color',[0.16 0.38 0.27]);
plot(gateSfluxT,'LineWidth',2,'Color',[0.75 0.75 0]) 
plot(ds1n,'bo'); plot(ds2n,'o','Color',[0 0.5 0]); plot(ds1n,'ro'); plot(ds1n,'co');
plot(-saltDHtot*kh3/kh,'o','Color',[0.75 0 0.75])
plot(-saltDiffVFlux,'o','Color',[0.16 0.38 0.27]);
plot(-gateSfluxT,'o','Color',[0.75 0.75 0]) 
title('Salt Fluxes','Fontsize',24)
legend('diffusion surface to \sigma=26','\sigma=26-\sigma=26.5','\sigma=26.5-\sigma=27','\sigma=27-\sigma=27.5','total horizontal diffusion','vertical diffusion','advection')
xlabel('time (days)','Fontsize',24')
ylabel('salt flux, g/s','Fontsize',24)
set(gca,'YScale','Log')
set(gca,'FontSize',24)

%% plot heat budget
load('gateAdvection26527.mat')
load('gateAdvection26265.mat')
load('gateAdvection27275.mat')
load('gateAdvectionSurface.mat')
%lagrangeWAGheatInt(isnan(lagrangeWAGheatInt))=0;
gateHfluxT=gateHflux26265+gateHflux26527+gateHflux27275+gateHfluxS26(1:138);
advectiveOffset=1.12e14*sign(gateHfluxT);
heatDiffVFlux(heatDiffVFlux==0)=nan;
figure; plot(gateHfluxT+advectiveOffset,'LineWidth',2,'Color',[0.75 0.75 0])
hold all
plot(heatDiffVFlux,'LineWidth',2,'Color',[0.16 0.38 0.27])
plot(heatDHtot*kh4/kh,'LineWidth',2,'Color',[0.75 0 0.75])
plot(lagrangeWAGheatInt,'LineWidth',2)%made in surfaceForcingForBudgets.m
plot(lagrangeWAGheatInt(1:131)+gateHfluxT(1:131)+advectiveOffset(1:131)+heatDiffVFlux+heatDHtot(1:131).'*kh4/kh,'k','LineWidth',2)
legend('advection','vertical diffusion','horizontal diffusion','surface forcing','total')
title('Heat Fluxes','FontSize',24); xlabel('time (days)','FontSize',24); ylabel('heat flux (W)','FontSize',24)
set(gca,'FontSize',24)

%% horizontal diffusion: heat, surface to sigma=26, loading
% kh=10;% 1e-9;%m^2/s numerical applied
% %kh2=100;%reasonable estimate maybe for ocean?? okubo suggests 10m^2/s for lengthscales of 10km
% 
% load('varying148ts16levelsRho.mat','T')
% temp=squeeze(T(:,:,10,:));
% clear T
% load('varying148ts16levelsRho.mat','Rho')
% rho=squeeze(Rho(:,:,10,:));
% clear Rho
% load('cp16levels148.mat');
% cp=squeeze(cp(:,:,10,:));
% load('distancesAreas.mat','DXG')
% load('distancesAreas.mat','DYG')
% %% horizontal diffusion: salt, surface to sigma=26
% heatDiffHFluxS275=zeros([148,1]);
% %holdvar1=zeros([700,200]);
% for k=1:148
%     k
%     csal=temp(:,:,k).*rho(:,:,k).*cp(:,:,k);
%     dZ=iso275depth(:,:,k);
% for i=200:500
%     for j=30:199
%         if inWAG2(i,j)
%             if ~inWAG2(i-1,j) &&inWAG2(i+1,j) &&csal(i-11,j)>0 &&csal(i+10,j)>0%left open
%                 %holdvar1(i,j,k)=-kh*(csal(i+10,j)-csal(i-11,j))*dZ(i,j)*DYG(i,j)/sum(DXG(i-11:i+10,j));
%                 heatDiffHFluxS275(k)=saltDiffHFlux(k)-kh*(csal(i+10,j)-csal(i-11,j))*dZ(i,j)*DYG(i,j)/sum(DXG(i-11:i+10,j));
%             end
%             
%             if ~inWAG2(i+1,j) &&inWAG2(i-1,j) &&csal(i+11,j)>0 &&csal(i-10,j)>0%right open
%                 heatDiffHFluxS275(k)=saltDiffHFlux(k)-kh*(csal(i-10,j)-csal(i+11,j))*dZ(i,j)*DYG(i,j)/sum(DXG(i-10:i+11,j));
%             end
%             
%             if ~inWAG2(i,j+1) &&inWAG2(i,j-1) &&csal(i,j-10)>0 &&csal(i,j+11)>0 %top open
%                 heatDiffHFluxS275(k)=saltDiffHFlux(k)-kh*(csal(i,j-10)-csal(i,j+11))*dZ(i,j)*DXG(i,j)/DYG(i+1,j);
%             end
%             
%         end
%     end
% end
% heatDiffHFluxS275(k)
% end

%% vertical salt flux, load
load('tsrhocpNativeGridIsoDepth28.mat')
 t28=tIso; r28=rhoIso; cp28=cpIso;
load('tsrhocpNativeGridIsoDepth275.mat')
 t275=tIso; r275=rhoIso; cp275=cpIso; clear sIso;
%% vertical salt flux, calculate
kv=1e-5;
for i=1:148
   csalgrad=(t275(:,:,i).*r275(:,:,i).*cp275(:,:,i)-t28(:,:,i).*r28(:,:,i).*cp28(:,:,i))./(iso28depth(:,:,i)-iso275depth(:,:,i));
   saltDiffVFlux(i)=-kv*nansum(csalgrad(inWAG2).*RAC(inWAG2)); %g salt/s
    
end
figure; plot(saltDiffVFlux)

%% total salt and heat, eulerian, very rough estimate
for k=1:148
    csal=sal(:,:,k).*rho(:,:,k);
    cheat=temp(:,:,k).*rho(:,:,k).*cp(:,:,k);
    dZ=iso275depth(:,:,k);
    saltEst(k)=nansum(nansum(csal.*dZ));
    heatEst(k)=nansum(nansum(cheat.*dZ));
end
figure; plot(saltEst); hold all; plot(diff(saltEst));
figure; plot(heatEst); hold all; plot(diff(heatEst));