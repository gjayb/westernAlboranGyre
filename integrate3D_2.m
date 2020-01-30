addpath('../mStuff')
% load('uDaily175.mat','U')
% load('uDaily175.mat','xvel')
% load('uDaily175.mat','yvel')
% load('vDaily175.mat','V')
% load('wDaily175.mat','W')
%load('rhoVary175.mat')
%load('uvwDailyNativeNF.mat', 'U','V','W')
%[nx,ny,nz,nt]=size(U)
%size(V)
%size(W)
% load('uvwSSHDailyDepth1rotated148F.mat', 'Urot')
% load('uvwSSHDailyDepth1rotated148F.mat', 'Vrot')
% Urot=Urot(:,:,1:46);
% Vrot=Vrot(:,:,1:46);%using daily flows as depths bc this is a test
% load('uvwSSHDailyDepth1rotated148F.mat', 'XC')
% load('uvwSSHDailyDepth1rotated148F.mat', 'YC')
load('geometrySpinupSteady.mat')
dBin=0.5.*dInterface(1:end-1)+0.5.*dInterface(2:end);
zForW=-1*dInterface(1:end-1);

xmin=min(min(XC));
ymin=min(min(YC));

xcM=111000*cosd(YC).*(XC-xmin*ones(size(XC)));
ycM=111000*(YC-ymin.*ones(size(YC)));
xuM=111000*cosd(YU).*(XU-xmin*ones(size(XU)));
yuM=111000*(YU-ymin.*ones(size(YU)));
xvM=111000*cosd(YV).*(XV-xmin*ones(size(XV)));
yvM=111000*(YV-ymin.*ones(size(YV)));

xarray=0:2000:max(xcM(:));
yarray=0:2000:max(ycM(:));
[xgrid,ygrid]=meshgrid(xarray,yarray);

% Urot=U.*repmat(AngleCS,[1 1 46 nt])-V.*repmat(AngleSN,[1 1 46 nt]); 
% Vrot=U.*repmat(AngleSN,[1 1 46 nt])+V.*repmat(AngleCS,[1 1 46 nt]); clear U V
% 
% 
% for ti=1:nt
%     ti
% for k=1:46
%    k
%    U(:,:,k,ti)=griddata(xuM,yuM,Urot(:,:,k,ti),xgrid,ygrid);
%    V(:,:,k,ti)=griddata(xvM,yvM,Vrot(:,:,k,ti),xgrid,ygrid);
%    W2(:,:,k,ti)=griddata(xcM,ycM,W(:,:,k,ti),xgrid,ygrid);
% end
% end
% clear W; W=W2; clear W2
% disp('done interpolating')
% save('uvw3Dgrid.mat','U','V','W','-v7.3')
load('uvw3Dgrid.mat')
[nx,ny,nz,nt]=size(U)
tvel=0:86400:(nt-1)*86400;
disp('done loading')
% U=repmat(U,[1 1 1 2]);
% V=repmat(V,[1 1 1 2]);
% W=repmat(W,[1 1 1 2]);
% xin0=xcM(200,100):3000:xcM(450,100);
% yin0=ycM(400,1):3000:ycM(400,200);
% [xin1,yin1]=meshgrid(xin0,yin0);
% xin=repmat(xin1,[1 1 12]);
% yin=repmat(yin1,[1 1 12]);
% zin1(1,1,1:12)=[0 -5 -10 -15 -20 -25 -30 -35 -40 -70 -75 -80];
% zin=repmat(zin1,[size(xin1) 1]);
% 
% r0=[xin(:).'; yin(:).'; zin(:).'];

% xin0=[xcM(250,100)+100 xcM(450,100)-100];
% yin0=[ycM(250,100)+100 ycM(450,100)-100];
% [xin1,yin1]=meshgrid(xin0,yin0);
xin=xcM(200:50:500,1:50:200);
yin=ycM(200:50:500,1:50:200);
dpos=d(200:50:500,1:50:200)>0;
xin=xin(dpos);
yin=yin(dpos);
load('iso275depthNFrev.mat')




%%
for j=27
    j
    t1=86400*(j):7200:(j+1)*86400;
    t2=(j)*86400:-7200:(j-1)*86400;
   for k=6:8
	k
	o1=k; o2=k+3;
    options=odeset('RelTol',10^(-o1),'AbsTol',10^(-o2));
    zin=-isoDepth(200:50:500,1:50:200,j);
    zin=zin(dpos);
    r0=[xin(:).'; yin(:).'; zin(:).'];
    disp('int forward')
    [~,zz]=ode45(@HamEqSolver_TriLin,t1,r0(:),options,U,V,W,xarray,yarray,zForW,tvel);
    xf1=zz(:,1:3:end-2);
    yf1=zz(:,2:3:end-1);
    zf1=zz(:,3:3:end);
    clear zz
    
    xf=reshape(xf1(end,:),size(xin));
    yf=reshape(yf1(end,:),size(yin));
    zf=reshape(zf1(end,:),size(zin));
    
%     latF79(:,:,:,j-78)=ones(size(yf1)).*ymin+yf1./111000;
%     lonF79(:,:,:,j-78)=ones(size(xf1)).*xmin+xf1./(111000.*cosd(lattr));
%     zF79(:,:,:,j-78)=zf;
    
   % disp('int backward')
   % [~,zz]=ode45(@HamEqSolver_TriLin,t2,r0(:),options,U,V,W,xvel,yvel,zForW,tvel);
   % xb1=zz(:,1:3:end-2);
   % yb1=zz(:,2:3:end-1);
   % zb1=zz(:,3:3:end);
   % clear zz
    
   % xb=reshape(xf1(end,:),size(xin));
   % yb=reshape(yf1(end,:),size(yin));
   % zb=reshape(zf1(end,:),size(zin));
    
%     latB79(:,:,:,j-78)=ones(size(yf1)).*ymin+yf1./111000;
%     lonB79(:,:,:,j-78)=ones(size(xf1)).*xmin+xf1./(111000.*cosd(lattr));
%     zB79(:,:,:,j-78)=zf;
    
   fn=strcat('traj3Diso27sInt1day27option',num2str(k),'.mat');
   save(fn,'*f*')
   clear *f1
   end
end
%fn='results3Dint8day79to89.mat';
%save(fn,'-v7.3')
%%

