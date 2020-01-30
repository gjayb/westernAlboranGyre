addpath('/nobackup1/gbrett/mStuff/')
  load('geometrySpinupSteady.mat')
xmin=min(min(XC));
ymin=min(min(YC));
xcm=(XC-xmin).*111000.*cosd(YC);
ycm=(YC-ymin).*111000;
zcm=-0.5*dInterface(1:end-1)-0.5*dInterface(2:end);
minZ=min(zcm)
dZ=diff(dInterface);
size(dZ)
load('varyingSigmaGradientNF.mat','gsx')
gsx=gsx(:,:,1:20,:);
load('varyingSigmaGradientNF.mat','gsy')
gsy=gsy(:,:,1:20,:);
load('varyingSigmaGradientNF.mat','gsz')
gsz=gsz(:,:,1:20,:);
size1=size(gsz)
%%
isoStr='275'
fn=strcat('uvwNativeGridIsoDepth',isoStr,'.mat');
load(fn,'isoDepth')
size2=size(isoDepth)
maxDepth=max(isoDepth(:))
for k=1:150
   k
for i=1:700
   for j=1:200
       if isoDepth(i,j,k)>0
       gsxI(i,j,k)=interp1(zcm(1:20),squeeze(gsx(i,j,:,k)),-isoDepth(i,j,k));
       gsyI(i,j,k)=interp1(zcm(1:20),squeeze(gsy(i,j,:,k)),-isoDepth(i,j,k));
       gszI(i,j,k)=interp1(zcm(1:20),squeeze(gsz(i,j,:,k)),-isoDepth(i,j,k));
       index1=find(dInterface<=isoDepth(i,j,k),1,'last');
       dZi(i,j,k)=dZ(index1);
	end
   end
end
end
size3=size(gszI)
clear gsx gsy gsz
fns=strcat('sigma',isoStr,'sigmagradientNF.mat')
save(fns, '*I','dZi')
%I am projecting vec(u) onto grad(sigma)
disp('sigmaGradMade')