%blah lagrange Wag Budget cross-isopycnal velocity blah
% 

%%
%load('isoDepthsNFsnap.mat','isoDepth275')
%isoDepthS=isoDepth;
load('uvwNativeGridIsoDepth275rev.mat')
load('sigma275sigmagradientNFrev.mat')
load('inWagManifolds.mat')
load('geometrySpinupSteady.mat','Angle*','*C*')
 %load('iso275depthNRsnap.mat')
 %isoDepthS=isoDepth;
% load('uvwNativeGridIsoDepth275NR.mat','*Iso','isoDepth')
load('lagrangeWAGboundaryNF.mat')
% load('geometrySpinupSteady.mat','XC','YC','dInterface')
%load('sigma275sigmagradientNR.mat', '*I','dZi')
%load('geometrySpinupSteady','Angle*')
load('distancesAreas','RAC','dxg','dyg')
DXG=reshape(dxg,[700 200]);
DYG=reshape(dyg,[700 200]);
%I rotate uIso and vIso to eastward and northward to get correct
%product- orthogonal vectors in same directions
[nxg,nyg,ntg]=size(gsxI)
[~,~,ntv]=size(uIso)
nt=min(ntg,ntv)
Urot=uIso(:,:,1:nt).*repmat(AngleCS,[1 1 nt]) - vIso(:,:,1:nt).*repmat(AngleSN,[1 1 nt]);  
Vrot=uIso(:,:,1:nt).*repmat(AngleSN,[1 1 nt]) + vIso(:,:,1:nt).*repmat(AngleCS,[1 1 nt]); 
clear uIso vIso
%%
%dot product is magnitude of cross-isopycnal velocity
gsxI(700,200,ntg)=0; gsyI(700,200,ntg)=0; gszI(700,200,ntg)=0; dZi(700,200,ntg)=0;
mag1=sqrt(gsxI(:,:,1:nt).^2+gsyI(:,:,1:nt).^2+gszI(:,:,1:nt).^2);%magnitude of grad(sigma)
sizemag1=size(mag1)
vecU=[reshape(Urot,1,[]);reshape(Vrot,1,[]);reshape(wIso(:,:,1:nt),1,[])];
vecGI=[reshape(gsxI(:,:,1:nt)./mag1,1,[]);reshape(gsyI(:,:,1:nt)./mag1,1,[]);reshape(gszI(:,:,1:nt)./mag1,1,[])];
sizevecu=size(vecU)
sizevecgi=size(vecGI)
uDotGi=reshape(dot(vecU,vecGI),[700 200 nt]);%magnitude of projected vec(u)
sizeudot=size(uDotGi)
%east-north-up cross-iso vel components
udx=uDotGi.*gsxI(:,:,1:nt)./mag1;
sizeudx=size(udx)
udy=uDotGi.*gsyI(:,:,1:nt)./mag1;
udz=uDotGi.*gszI(:,:,1:nt)./mag1;
magCrossIsoVel=sqrt(udx.^2+udy.^2+udz.^2);
sizemagcross=size(magCrossIsoVel)
%native grid cross-iso vel components
udxN=udx.*repmat(AngleCS,[1 1 nt])+udy.*repmat(AngleSN,[1 1 nt]);
udyN=-udx.*repmat(AngleSN,[1 1 nt])+udy.*repmat(AngleCS,[1 1 nt]);
%uTest=Urot.*repmat(AngleCS,[1 1 131])+Vrot.*repmat(AngleSN,[1 1 131]);
%vTest=-Urot.*repmat(AngleSN,[1 1 131])+Vrot.*repmat(AngleCS,[1 1 131]);
% figure; quiver(XC(1:20:end,1:20:end),YC(1:20:end,1:20:end),uIso(1:20:end,1:20:end,20),vIso(1:20:end,1:20:end,20),'k');
% hold on;  quiver(XC(1:20:end,1:20:end),YC(1:20:end,1:20:end),uTest(1:20:end,1:20:end,20),vTest(1:20:end,1:20:end,20),'r');

[~,~,nti]=size(inWAG27)
nt2=min(nt,nti)
[~,~,nti2]=size(inWAG275)
nt3=min(nt,nti2)
%projecting onto grad(sigma) is outward direction, so negative is inward
%vertVol275=squeeze(nansum(nansum((-udz(:,:,1:nt2).*repmat(RAC,[1 1 nt2])-udxN(:,:,1:nt2).*dZi(:,:,1:nt2).*repmat(DYG,[1 1 nt2])-udyN(:,:,1:nt2).*dZi(:,:,1:nt2).*repmat(DXG,[1 1 nt2])).*double(inWAG27(:,:,1:nt2))))); %
% vecU=[reshape(uIso(:,:,1:131),1,[]);reshape(vIso(:,:,1:131),1,[]);reshape(wIso(:,:,1:131),1,[])];
% vecGI=[reshape(gsxI,1,[]);reshape(gsyI,1,[]);reshape(gszI,1,[])];
% uDotGi=reshape(dot(vecU,vecGI),[700 200 131]);
% vertVol275v6=-squeeze(nansum(nansum((uDotGi(:,:,1:nt2).*repmat(RAC,[1 1 nt2])).*double(inWAG27(:,:,1:nt2)))));
 
 %% new 10/31/2017
 snapDiff=diff(isoDepthS(:,:,1:nt2+1),1,3);%positive difference is deepening
 snapVel=snapDiff.*gszI(:,:,1:nt2)./(mag1(:,:,1:nt2).*86400);%positive down
 surfaceZ=repmat(RAC,[1 1 nt2]);
size(DYG)
size(DXG)
 surfaceXN=dZi(:,:,1:nt2).*repmat(DYG,[1 1 nt2]);
 surfaceYN=dZi(:,:,1:nt2).*repmat(DXG,[1 1 nt2]);
 surfaceX=surfaceXN.*repmat(AngleCS,[1 1 nt2]) - surfaceYN.*repmat(AngleSN,[1 1 nt2]);  
 surfaceY=surfaceXN.*repmat(AngleSN,[1 1 nt2]) + surfaceYN.*repmat(AngleCS,[1 1 nt2]); 
sizesY=size(surfaceY) 
vecSurf=[reshape(surfaceX,1,[]);reshape(surfaceY,1,[]);reshape(surfaceZ,1,[])];%all positive
vecGI2=[reshape(gsxI(:,:,1:nt2)./mag1(:,:,1:nt2),1,[]);reshape(gsyI(:,:,1:nt2)./mag1(:,:,1:nt2),1,[]);reshape(gszI(:,:,1:nt2)./mag1(:,:,1:nt2),1,[])];
surfaces=reshape(dot(vecSurf,vecGI2),[700 200 nt2]);%vecGI usually negative
 vertVelCross=uDotGi(:,:,1:nt2)-snapVel;%uDotGi positive out of gyre, generally down, total positive down
 vertVol=vertVelCross.*surfaces;%this is positive down*negative surface areas, giving Positive UP

 vertVol275inWag27=squeeze(nansum(nansum(vertVol(:,:,1:nt3).*double(inWAG27f(:,:,1:nt3)))));%mean -7e5
 %vertVol275inWag27v2=squeeze(nansum(nansum((uDotGi(:,:,1:nt2)+snapVel).*surfaces.*double(inWAG27f(:,:,1:nt3)))));
 vertVol275inWag275=squeeze(nansum(nansum(vertVol(:,:,1:nti).*double(inWAG275f(:,:,1:nti)))));
 save('vertVolIso275filled.mat','vertVelCross','vertVol','vertVol*')
 %% 12/20 indexing
   load('lagrangeWAGboundaryPortionFilled.mat', 'inWAG275f','inWAG27f')
   load('vertVolIso275rev.mat', 'vertVol')
 vertVol275inWag275i(2:147)=squeeze(nansum(nansum(vertVol(:,:,2:147).*inWAG275f(:,:,1:146)./5)));
 vertVol275inWag27i(2:147)=squeeze(nansum(nansum(vertVol(:,:,2:147).*inWAG27f(:,:,1:146)./5)));
 load('vertVolIso275filledNative.mat', 'crossIsoVol')
  vertVol275inWag275n(2:147)=squeeze(nansum(nansum(crossIsoVol(:,:,2:147).*inWAG275f(:,:,1:146)./5)));
 vertVol275inWag27n(2:147)=squeeze(nansum(nansum(crossIsoVol(:,:,2:147).*inWAG27f(:,:,1:146)./5)));
