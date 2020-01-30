%% velocities rotated, moved to cell center, regridded
addpath('/nobackup1/gbrett/mStuff')
load('geometrySpinupSteady.mat');
dBin=0.5.*dInterface(1:end-1)+0.5.*dInterface(2:end);
isopycs=[26.3 26.5 26.75 27 27.5 28 28.5 28.9 29 29.1];
isopycStr=[263 265 2675 27 275 28 285 289 29 291];
%%
for iiso=1:length(isopycs)
    isopyc=isopycs(iiso)
    fnL=strcat('uvwNativeGridIsoDepth',num2str(isopycStr(iiso)),'.mat');
    load(fnL);
[nx,ny,nt]=size(uIso)
for i=1:nt
   i
   uIso(:,:,i) = griddata(XU,YU,uIso(:,:,i),XC,YC);
   vIso(:,:,i) = griddata(XV,YV,vIso(:,:,i),XC,YC);
end
   u=uIso.*repmat(AngleCS,[1 1 nt]) - vIso.*repmat(AngleSN,[1 1 nt]);  
   v=uIso.*repmat(AngleSN,[1 1 nt]) + vIso.*repmat(AngleCS,[1 1 nt]); 

clear uIso vIso
disp('section 1 done')
%% better grid

xmin=min(min(XC));
ymin=min(min(YC));
xinM=111000*cosd(YC).*(XC-xmin*ones(size(XC)));
yinM=111000*(YC-ymin.*ones(size(YC)));
xvel=0:1000:max(max(xinM));
yvel=0:1000:max(max(yinM));

[xvelg,yvelg]=meshgrid(xvel,yvel);
[ng,mg]=size(xvelg);
W=zeros([ng mg nt]);
U=W; V=W;
%clear u v
%V=U; W=U;
for k=1:nt
    k
%    for j=1:46
        U(:,:,k)=griddata(xinM,yinM,u(:,:,k),xvelg,yvelg);
        V(:,:,k)=griddata(xinM,yinM,v(:,:,k),xvelg,yvelg);
        W(:,:,k)=griddata(xinM,yinM,wIso(:,:,k),xvelg,yvelg);
%    end
end
clear u v wIso
%% save
disp('saving')
fn=strcat('uvwIso',num2str(isopycStr(iiso)),'InterpNF.mat');
save(fn,'-v7.3')
disp('done iso')
end
disp('all done')
