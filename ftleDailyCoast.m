addpath('../mStuff')
load('geometrySpinupSteady.mat')

xmin=min(min(XC));
%xmax=max(max(XC(225:430,1:110)));
ymin=min(min(YC));
%ymax=max(max(YC(225:430,1:110)));
lonCM=(lonCoast(597:840)-xmin).*111000.*cosd(latCoast(597:840));
latCM=(latCoast(597:840)-ymin).*111000;

xin=repmat(lonCM,[1 5]); yin(:,1)=latCM; yin(:,2:5)=repmat(latCM,[1 4])+1000*repmat((1:4),[244 1]);
xin=xin.'; yin=yin.';
z0=[xin(:)';yin(:)'];
NX=244; NY=5;

disp('found initial points')

addpath('../mStuff')
load('uvAdepth1interpolated.mat', 'XC')
load('uvAdepth1interpolated.mat', 'YC')
load('uvAdepth1interpolated.mat', 'u')
load('uvAdepth1interpolated.mat', 'v')
load('uvAdepth1interpolated.mat', 'yvel')
load('uvAdepth1interpolated.mat', 'xvel')
load('uvAdepth1interpolated.mat', 'tvel')

disp('load done')

disp('entering day loop')

for j=3:150; %9:1:30 %0:0.5:21
    j
    clear *tr*
    tmesh1=86400*(j-2):43200:(j)*86400;
    tmesh2=(j+2)*86400:-43200:(j)*86400;
    options=odeset('RelTol',10^(-10),'AbsTol',10^(-13));
    
    
    disp('int forward')
    [~,zz]=ode45(@HamEqSolver_BiLin_Irina,tmesh1,z0(:),options,u,v,xvel,yvel,tvel);%expects 1D xvel yvel tvel where u,v are defined on the meshgrid
    xtr(:,j-2)=zz(end,1:2:end-1);
    ytr(:,j-2)=zz(end,2:2:end);
    
    xend=reshape(squeeze(xtr(:,j-2)),NY,NX); yend=reshape(squeeze(ytr(:,j-2)),NY,NX);

size1=size(xtr)
%size2=size(xin)
%size3=size(xend)
    
    dx11=xend(2:end-1,3:end)-xend(2:end-1,1:end-2); 
    dx12=xend(3:end,2:end-1)-xend(1:end-2,2:end-1);
    dy21=yend(2:end-1,3:end)-yend(2:end-1,1:end-2);
    dy22=yend(3:end,2:end-1)-yend(1:end-2,2:end-1);
    dx0 =xin(2:end-1,3:end)-xin(2:end-1,1:end-2);
    dy0 =yin(3:end,2:end-1)-yin(1:end-2,2:end-1);

%size4=size(dx11)
%size5=size(dx0)
tIntegrate=tmesh1(end)-tmesh1(1);
    a11=(dx11./dx0).^2.0+(dy21./dx0).*(dy21./dx0);	
    a12=(dx11./dx0).*(dx12./dy0)+(dy21./dx0).*(dy22./dy0);
    a21=(dx12./dy0).*(dx11./dx0)+(dy22./dy0).*(dy21./dx0);
    a22=(dx12./dy0).*(dx12./dy0)+(dy22./dy0).^2.0;
    lambda1=(a22+a11)/2.0+sqrt(((a22+a11)/2.0).^2.0-(a11.*a22-a12.*a21));
    lambda2=(a22+a11)/2.0-sqrt(((a22+a11)/2.0).^2.0-(a11.*a22-a12.*a21));
    ftleP(:,:,j-2)=log(sqrt(max(lambda1,lambda2)))./abs(tIntegrate);
    
    
    disp('int backward')
    [~,zz]=ode45(@HamEqSolver_BiLin_Irina,tmesh2,z0(:),options,u,v,xvel,yvel,tvel);%expects 1D xvel yvel tvel where u,v are defined on the meshgrid
    xtr2(:,j-2)=zz(end,1:2:end-1);
    ytr2(:,j-2)=zz(end,2:2:end);
disp('check 1')    
    xend=reshape(xtr2(:,j-2),NY,NX); yend2=reshape(ytr(:,j-2),NY,NX);
    
    dx11=xend(2:end-1,3:end)-xend(2:end-1,1:end-2); 
    dx12=xend(3:end,2:end-1)-xend(1:end-2,2:end-1);
    dy21=yend(2:end-1,3:end)-yend(2:end-1,1:end-2);
    dy22=yend(3:end,2:end-1)-yend(1:end-2,2:end-1);
    dx0 =xin(2:end-1,3:end)-xin(2:end-1,1:end-2);
    dy0 =yin(3:end,2:end-1)-yin(1:end-2,2:end-1);
disp('check 2')
    a11=(dx11./dx0).^2.0+(dy21./dx0).*(dy21./dx0);	
    a12=(dx11./dx0).*(dx12./dy0)+(dy21./dx0).*(dy22./dy0);
    a21=(dx12./dy0).*(dx11./dx0)+(dy22./dy0).*(dy21./dx0);
    a22=(dx12./dy0).*(dx12./dy0)+(dy22./dy0).^2.0;
    lambda1=(a22+a11)/2.0+sqrt(((a22+a11)/2.0).^2.0-(a11.*a22-a12.*a21));
    lambda2=(a22+a11)/2.0-sqrt(((a22+a11)/2.0).^2.0-(a11.*a22-a12.*a21));
    ftleN(:,:,j-2)=log(sqrt(max(lambda1,lambda2)))./abs(tIntegrate);
end%day loop
disp('loop done, saving')
save('ftleCoast3to150.mat','-v7.3')
disp('done')
