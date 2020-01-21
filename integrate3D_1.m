addpath('../mStuff')
load('uDaily175.mat','U')
size(U)
U=U(:,:,1:12,60:100);

load('uDaily175.mat','xvel')
size(xvel)
load('uDaily175.mat','yvel')
size(yvel)

%load('rhoVary175.mat')
load('geometrySpinupSteady.mat')
dBin=0.5.*dInterface(1:end-1)+0.5.*dInterface(2:end);
zForW=-1*dInterface(1:end-1);
tvel=0:86400:147*86400;
zForW=zForW(1:12);
tvel=tvel(60:100);

disp('loading v,w')
load('vDaily175.mat','V')
V=V(:,:,1:12,60:100);
load('wDaily175.mat','W')
W=W(:,:,1:12,60:100);

xmin=min(min(XC));
ymin=min(min(YC));

xcM=111000*cosd(YC).*(XC-xmin*ones(size(XC)));
ycM=111000*(YC-ymin.*ones(size(YC)));

xin0=xcM(200,100):100000:xcM(450,100);
yin0=ycM(400,1):100000:ycM(400,200);
length(xin0)
length(yin0)
[xin1,yin1]=meshgrid(xin0,yin0);
xin=repmat(xin1,[1 1 5]);
yin=repmat(yin1,[1 1 5]);
zin1(1,1,1:5)=[0 -5 -10 -15 -20];% -25 -30 -35 -40 -75 -80 -85];
zin=repmat(zin1,[size(xin1) 1]);

r0=[xin(:).'; yin(:).'; zin(:).'];

for j=79:89
    j
    t1=86400*(j):43200:(j+8)*86400;
    t2=(j)*86400:-43200:(j-8)*86400;
    options=odeset('RelTol',10^(-3),'AbsTol',10^(-6));
    
    disp('int forward')
    [~,zz]=ode45(@HamEqSolver_TriLin,t1,r0(:),options,U,V,W,xvel,yvel,zForW,tvel);
    xf1=zz(:,1:3:end-2);
    yf1=zz(:,2:3:end-1);
    zf1=zz(:,3:3:end);
    clear zz
    
    xf=reshape(xf1(end,:),size(xin));
    yf=reshape(yf1(end,:),size(yin));
    zf=reshape(zf1(end,:),size(zin));
    
    latF79(:,:,:,j-78)=ones(size(yf1)).*ymin+yf1./111000;
    lonF79(:,:,:,j-78)=ones(size(xf1)).*xmin+xf1./(111000.*cosd(latF79(:,:,:,j-78)));
    zF79(:,:,:,j-78)=zf;
    
    disp('int backward')
    [~,zz]=ode45(@HamEqSolver_TriLin,t2,r0(:),options,U,V,W,xvel,yvel,zForW,tvel);
    xb1=zz(:,1:3:end-2);
    yb1=zz(:,2:3:end-1);
    zb1=zz(:,3:3:end);
    clear zz
    
    xb=reshape(xf1(end,:),size(xin));
    yb=reshape(yf1(end,:),size(yin));
    zb=reshape(zf1(end,:),size(zin));
    
    latB79(:,:,:,j-78)=ones(size(yf1)).*ymin+yf1./111000;
    lonB79(:,:,:,j-78)=ones(size(xf1)).*xmin+xf1./(111000.*cosd(latB79(:,:,:,j-78)));
    zB79(:,:,:,j-78)=zf;
    
    fn=strcat('traj3Dint8day',num2str(j),'.mat');
    save(fn,'*f1','-v7.3')
    clear *f1
end
fn='results3Dint8day79to89.mat';
save(fn,'-v7.3')
