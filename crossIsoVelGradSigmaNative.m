% %doing cross-iso volume flux in native grid to see differences
% 
% load('saTCtCpSigma162NF.mat', 'Sigma')
% Sigma=Sigma(:,:,1:22,1:150);
% [nx,ny,nz,nt]=size(Sigma);
% load('geometrySpinupSteady.mat')
% xmin=min(min(XC));
% ymin=min(min(YC));
% xcm=(XC-xmin).*111000.*cosd(YC);
% ycm=(YC-ymin).*111000;
% zcm=-0.5*dInterface(1:end-1)-0.5*dInterface(2:end);
% minZ=min(zcm)
% dZ=diff(zcm);
% dZ=dZ(1:22); zcm=zcm(1:22);
% load('distancesAreas.mat')
% DXC=reshape(dxc,[700 200]);
% DYC=reshape(dyc,[700 200]);
% disp('prep section done')
% %%
% dSigdx=diff(Sigma(:,1:end-1,1:end-1,:),1,1)./repmat(DXC(1:end-1,1:end-1),[1 1 nz-1 nt]);
% dSigdy=diff(Sigma(1:end-1,:,1:end-1,:),1,2)./repmat(DYC(1:end-1,1:end-1),[1 1 nz-1 nt]);
% dZ3(1,1,1:22)=dZ;
% dSigdz=diff(Sigma(1:end-1,1:end-1,:,:),1,3)./repmat(dZ3(1:end-1),[699 199 1 nt]);
% 
% save('varyingSigmaGradientNFnative.mat','dS*','-v7.3')
% disp('varyingSigmaGrad saved')
% clear Sigma
%%
isoStr='275'
fn=strcat('uvwNativeGridIsoDepth',isoStr,'rev.mat');
load(fn,'isoDepth')
size2=size(isoDepth)
maxDepth=max(isoDepth(:))
[~,~,nt]=size(isoDepth);
[~,~,~,nt2]=size(dSigdx);
dZ=diff(dInterface)
for k=1:min(nt,nt2)
   k
for i=1:699
   for j=1:199
       if isoDepth(i,j,k)>0
       gsxI(i,j,k)=interp1(zcm(1:20),squeeze(dSigdx(i,j,:,k)),-isoDepth(i,j,k));
       gsyI(i,j,k)=interp1(zcm(1:20),squeeze(dSigdy(i,j,:,k)),-isoDepth(i,j,k));
       gszI(i,j,k)=interp1(zcm(1:20),squeeze(dSigdz(i,j,:,k)),-isoDepth(i,j,k));
       index1=find(dInterface<=isoDepth(i,j,k),1,'last');
       dZi(i,j,k)=dZ(index1);
	end
   end
end
end
size3=size(gszI)

fns=strcat('sigma',isoStr,'sigmagradientNFnative.mat')
save(fns, '*I','dZi')
%I am projecting vec(u) onto grad(sigma)
clear dSig*
disp('sigmaGradMade')
%%
isoStr='275'
fn=strcat('uvwNativeGridIsoDepth',isoStr,'rev.mat');
load(fn,'*Iso')
%%
[~,~,nt0]=size(uIso);
[~,~,nt1]=size(gsxI);
nt2=min(nt0,nt1);
gsxI(700,200,nt1)=0;
gsyI(700,200,nt1)=0;
gszI(700,200,nt1)=0;
dZi(700,200,nt1)=0;
surfaceZ=repmat(RAC,[1 1 nt2]);
surfaceXN=dZi(:,:,1:nt2).*repmat(DYG,[1 1 nt2]);
surfaceYN=dZi(:,:,1:nt2).*repmat(DXG,[1 1 nt2]);
vecUvol=[reshape(uIso(:,:,1:nt2).*surfaceXN,1,[]);reshape(vIso(:,:,1:nt2).*surfaceYN,1,[]);reshape(wIso(:,:,1:nt2).*surfaceZ,1,[])];
mag1=sqrt(gsxI(:,:,1:nt2).^2+gsyI(:,:,1:nt2).^2+gszI(:,:,1:nt2).^2);%magnitude of grad(sigma)
vecGI=[reshape(gsxI(:,:,1:nt2)./mag1,1,[]);reshape(gsyI(:,:,1:nt2)./mag1,1,[]);reshape(gszI(:,:,1:nt2)./mag1,1,[])];
uvolDotGi=reshape(dot(vecUvol,vecGI),[700 200 nt2]);%magnitude of projected vec(u)
magUvol=sqrt((uIso(:,:,1:nt2).*surfaceXN).^2+(vIso(:,:,1:nt2).*surfaceYN).^2+(wIso(:,:,1:nt2).*surfaceZ).^2);
portionVolCross=uvolDotGi./magUvol;
magWvol=abs(wIso(:,:,1:nt2).*surfaceZ);
portionVolW=uvolDotGi./magWvol;
disp('cross-iso direction full vel made')
%%
load('iso275depthNFsnapRev.mat')
isoDepthS=isoDepth;

snapDiff=diff(isoDepthS(:,:,1:nt2+1),1,3);%positive means deepening means out of gyre
snapVol=snapDiff.*surfaceZ.*gszI(:,:,1:nt2)./(mag1(:,:,1:nt2).*86400);
crossIsoVol=uvolDotGi-snapVol;%grad sigma points to denser water, downward, out of gyre
%uvolDotGi is thus positive out of gyre (~down)
%snapDiff is also positive out of gyre
%difference is positive out of gyre
disp('cross-iso vel made')
%%
load('lagrangeWAGboundaryPortionFilled.mat','inWAG275f','inWAG27f')
[~,~,nt]=size(inWAG275f);
nt3=min(nt2,nt);
vertVol275inWag27=squeeze(nansum(nansum(crossIsoVol(:,:,1:nt3).*inWAG27f(:,:,1:nt3)./5)));
vertVol275inWag275=squeeze(nansum(nansum(crossIsoVol(:,:,1:nt3).*inWAG275f(:,:,1:nt3)./5)));
portioncrossVolW=crossIsoVol./magWvol;
portioncrossVol=crossIsoVol./magUvol;
save('vertVolIso275filledNative.mat','crossIsoVol','vertVol*','portion*','uvolDotGi','snapVol')
disp('done')
%%





