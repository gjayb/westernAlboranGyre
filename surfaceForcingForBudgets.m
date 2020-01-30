%% load surface wind stress and heat flux fields

load('C:\Users\JayB\Documents\MATLAB\MITgcm\addForceFull\windForEkman.mat')

tauN148=tauN(:,:,iNov12007:i148days);
tauE148=tauE(:,:,iNov12007:i148days);

dEdY=zeros(size(tauN148));
dEdY(:,2:32,:,:)=0.5*(tauN148(:,3:33,:,:)-tauN148(:,1:31,:,:));
dEdY(:,[1 33],:,:)=0.5*tauN148(:,[1 33],:,:);
dEdY=dEdY/(0.25*110000*cosd(mean(double(latC))));

dNdX=zeros(size(tauE148));
dNdX(2:48,:,:,:)=0.5*(tauE148(3:49,:,:,:)-tauE148(1:47,:,:,:));
dNdX([1 49],:,:,:)=tauE148([2 49],:,:,:)-tauE148([1 48],:,:,:);
dNdX=dNdX/(0.25*110000);

curlNE=dNdX-dEdY;
%% load inWAG, grab single curve
load('wagAreaAndFluxGate.mat', 'inWAG','XC','YC')
load('distancesAreas.mat','RAC')
xmin=min(min(XC)); ymin=min(min(YC));
xcm=111000*cosd(YC).*(XC-xmin*ones(size(XC)));
ycm=111000*(YC-ymin*ones(size(YC)));
in1=sum(inWAG,3)>60;
%% steady area wind stress curl
for i=1:148
   curlNow=griddata(double(lonCg).',double(latCg).',mean(curlNE(:,:,(4*i-3):4*i),3),XC(in1),YC(in1));
   eulerWAGwindInt(i)=nansum(curlNow.*RAC(in1));
end

figure; plot(eulerWAGwindInt)
grid on
xlabel('simulation day'); ylabel('integrated wind stress curl, m^2/s')
title('Constant WAG boundary wind stress forcing')
%% changing area v2

load('wagAreaAndFluxGate.mat', 'inWAG')
load('uvwAve148.mat', 'XC')
load('uvwAve148.mat', 'YC')
load('distancesAreas.mat','RAC')
figure;
for i=9:148
    [c,h]=contour(XC,YC,double(inWAG(:,:,i)),[1 1]);
    if length(c)>0
    npointsi=find(c(2,:)==floor(c(2,:)));
    varhold=max(c(2,npointsi));
    i1=1+find(c(2,:)==varhold);
    iend=varhold+i1-1;
    inWAG2(:,:,i)=inpolygon(XC,YC,c(1,i1:iend),c(2,i1:iend));
    end
end

areaWAG2=squeeze(nansum(nansum(inWAG2.*repmat(RAC,[1 1 148]))));

%% changing area wind stress curl
inWAG2=logical(sign(inwagSareas));
for i=1:146
    holdvar=inwagSareas(:,:,i);
    
    if length(XC(inWAG2(:,:,i)))>0
   curlNow=griddata(double(lonCg).',double(latCg).',mean(curlNE(:,:,(4*i-3):4*i),3),XC(inWAG2(:,:,i)),YC(inWAG2(:,:,i)));
   lagrangeWAGwindInt(i)=nansum(curlNow.*holdvar(inWAG2(:,:,i)));
    else 
        lagrangeWAGwindInt(i)=nan;
    end
end

figure; plot(lagrangeWAGwindInt)
grid on
xlabel('simulation day'); ylabel('integrated wind stress curl, m^2/s')
title('Lagrangian WAG boundary wind stress forcing')
%% comparison plot
figure; plot(eulerWAGwindInt,'LineWidth',2); hold all; plot(lagrangeWAGwindInt,'r:','LineWidth',2)
grid on
xlabel('simulation day'); ylabel('integrated wind stress curl, m^2/s')
title('WAG wind stress forcing'); legend('eulerian','lagrangian')

%% load heat files
load('C:\Users\JayB\Desktop\climatologyIfremar\eraInterim200720082009.mat')
load('eraDataAndFall2007meansTry2.mat')
%12 hour frequency, hflux2 should be full heat flux, positive up (to atmos)
%hours since 1900 may be hours since 1990, check!
%% steady area heat flux, precip prep
%load('wagAreaAndFluxGate.mat', 'inWAG','XC','YC')
load('distancesAreas.mat','RAC')
xmin=min(min(XC)); ymin=min(min(YC));
xcm=111000*cosd(YC).*(XC-xmin*ones(size(XC)));
ycm=111000*(YC-ymin*ones(size(YC)));
%in1=sum(inWAG,3)>60;
iNov12007=608;
i148=904;%(march 28 2008)
h148=hflux2(:,:,iNov12007:i148);
ep148=EminusPms(:,:,iNov12007:i148);
%% steady area heat flux
for i=1:148
   hNow=griddata(double(lonEg).',double(latEg).',mean(h148(:,:,(2*i-1):2*i),3),XC(in1),YC(in1));
   eulerWAGheatInt(i)=-nansum(hNow.*RAC(in1));
end

figure; plot(eulerWAGheatInt)
grid on
xlabel('simulation day'); ylabel('heat flux')
title('Constant WAG boundary atmospheric heat flux, positive into ocean')
%% changing area heat flux

for i=1:146
    holdvar=inwagSareas(:,:,i);
    if length(XC(inWAG2(:,:,i)))>0
   hNow=griddata(double(lonEg).',double(latEg).',mean(h148(:,:,(2*i-1):2*i),3),XC(inWAG2(:,:,i)),YC(inWAG2(:,:,i)));
   lagrangeWAGheatInt(i)=-nansum(hNow.*holdvar(inWAG2(:,:,i)));
    else 
        lagrangeWAGheatInt(i)=nan;
    end
end

figure; plot(lagrangeWAGheatInt)
grid on
xlabel('simulation day'); ylabel('heat flux')
title('Lagrangian WAG boundary atmospheric heat flux, positive into ocean')
%% comparison plot
figure; plot(eulerWAGheatInt,'LineWidth',2); hold all; plot(lagrangeWAGheatInt,'r:','LineWidth',2)
grid on
xlabel('simulation day'); ylabel('integrated wind stress curl, m^2/s')
title('WAG atmospheric heat flux, positive into ocean'); legend('eulerian','lagrangian')




%% steady area precipitation- volume and mass flux
%EminusPms %evap-precip, m/s, loaded with heat flux files
%ep148 is prepped like h148, out of ocean is positive
for i=1:148
   vNow=griddata(double(lonEg).',double(latEg).',mean(ep148(:,:,(2*i-1):2*i),3),XC(in1),YC(in1));
   eulerWAGvolInt(i)=-nansum(vNow.*RAC(in1));%into ocean positive
   eulerWAGmassInt(i)=-1000*nansum(vNow.*RAC(in1));%into ocean positive
end

%% changing area precipitation- volume and mass flux
%load('diffFluxHorizS26.mat','inWAG2')
for i=1:146
    holdvar=inwagSareas(:,:,i);
    if length(XC(inWAG2(:,:,i)))>0
   vNow=griddata(double(lonEg).',double(latEg).',mean(ep148(:,:,(2*i-1):2*i),3),XC(inWAG2(:,:,i)),YC(inWAG2(:,:,i)));
   lagrangeWAGvolInt(i)=-nansum(vNow.*holdvar(inWAG2(:,:,i)));
   lagrangeWAGmassInt(i)=-1000*nansum(vNow.*RAC(inWAG2(:,:,i)));
    else 
        lagrangeWAGvolInt(i)=nan;
        lagrangeWAGmassInt(i)=nan;
    end
end

figure; plot(lagrangeWAGvolInt)
grid on
xlabel('simulation day'); ylabel('volume flux, P-E, m^3/s')
title('Lagrangian WAG boundary precip-evap, positive into ocean')

%%



