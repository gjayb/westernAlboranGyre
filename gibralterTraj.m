%addpath('../mStuff')
load('geometrySpinupSteady.mat')

xmin=min(min(XC));
%xmax=max(max(XC(225:430,1:110)));
ymin=min(min(YC));
%ymax=max(max(YC(225:430,1:110)));

%lonCM=(lonCoast(597:840)-xmin).*111000.*cosd(latCoast(597:840));
%latCM=(latCoast(597:840)-ymin).*111000;
xcm=(XC-xmin).*111000.*cosd(YC);
ycm=(YC-ymin).*111000;
%xin=repmat(lonCM,[1 5]); yin(:,1)=latCM; yin(:,2:5)=repmat(latCM,[1 4])+1000*repmat((1:4),[244 1]);
xinmin=min(min(xcm(200:230,:)));
xinmax=max(max(xcm(200:230,:)));
yinmin=min(min(ycm(200:230,:)));
yinmax=max(max(ycm(200:230,:)));
xin1=xinmin:3000:xinmax; yin1=yinmin:3000:yinmax;
[xin,yin]=meshgrid(xin1,yin1);
z0=[xin(:)';yin(:)'];
[NY,NX]=size(xin)

disp('found initial points')

addpath('../mStuff')
load('uvAdepth1interpolated.mat', 'XC')
load('uvAdepth1interpolated.mat', 'YC')

load('uvAdepth1interpolated.mat', 'yvel')
load('uvAdepth1interpolated.mat', 'xvel')
load('uvAdepth1interpolated.mat', 'tvel')

disp('load 1 done')

disp('entering day loop')

for j=101:150; 
    j
    load('uvAdepth1interpolated.mat', 'u')
load('uvAdepth1interpolated.mat', 'v')
   disp('load 2 done')
    clear *tr*
    tmesh1=86400*(j):43200:(j+12)*86400;
    tmesh2=(j)*86400:-43200:(j-2)*86400;
    options=odeset('RelTol',10^(-10),'AbsTol',10^(-13));
 tIntegrate=2*86400;   
    
 %if j==35
 %    u(:,:,34)=u(:,:,50); u(:,:,33)=u(:,:,49);
 %    v(:,:,34)=v(:,:,50); v(:,:,33)=v(:,:,49);
 %elseif j==36
 %    u(:,:,34)=u(:,:,50); 
 %    v(:,:,34)=v(:,:,50);
 %elseif j==49
 %    u(:,:,51)=u(:,:,35); 
 %    v(:,:,51)=v(:,:,35);
 %elseif j==50
 %    u(:,:,51)=u(:,:,35); u(:,:,52)=u(:,:,36);
 %    v(:,:,51)=v(:,:,35); v(:,:,52)=v(:,:,36);
 %end
 
 
    disp('int forward')
    [~,zz]=ode45(@HamEqSolver_BiLin_Irina,tmesh1,z0(:),options,u,v,xvel,yvel,tvel);%expects 1D xvel yvel tvel where u,v are defined on the meshgrid
    xtr=zz(:,1:2:end-1);
    ytr=zz(:,2:2:end);
    
%    xend=reshape(xtr(end,:,:),NY,NX); yend=reshape(ytr(end,:,:),NY,NX);
%s1=size(xend);
%s2=size(xin);    
%    dx11=xend(2:end-1,3:end)-xend(2:end-1,1:end-2); 
%    dx12=xend(3:end,2:end-1)-xend(1:end-2,2:end-1);
%    dy21=yend(2:end-1,3:end)-yend(2:end-1,1:end-2);
%    dy22=yend(3:end,2:end-1)-yend(1:end-2,2:end-1);
%    dx0 =xin(2:end-1,3:end)-xin(2:end-1,1:end-2);
%    dy0 =yin(3:end,2:end-1)-yin(1:end-2,2:end-1);%
%
%    a11=(dx11./dx0).^2.0+(dy21./dx0).*(dy21./dx0);	
%    a12=(dx11./dx0).*(dx12./dy0)+(dy21./dx0).*(dy22./dy0);
%    a21=(dx12./dy0).*(dx11./dx0)+(dy22./dy0).*(dy21./dx0);
%    a22=(dx12./dy0).*(dx12./dy0)+(dy22./dy0).^2.0;
%    lambda1=(a22+a11)/2.0+sqrt(((a22+a11)/2.0).^2.0-(a11.*a22-a12.*a21));
%    lambda2=(a22+a11)/2.0-sqrt(((a22+a11)/2.0).^2.0-(a11.*a22-a12.*a21));
%    ftleP=log(sqrt(max(lambda1,lambda2)))./abs(tIntegrate);
%    
%    disp('int backward')
%    [~,zz]=ode45(@HamEqSolver_BiLin_Irina,tmesh2,z0(:),options,u,v,xvel,yvel,tvel);%expects 1D xvel yvel tvel where u,v are defined on the meshgrid
%    xtr2=zz(:,1:2:end-1);
%    ytr2=zz(:,2:2:end);
%    
%    xend=reshape(xtr2(end,:,:),NY,NX); yend2=reshape(ytr(end,:,:),NY,NX);
%    
%    dx11=xend(2:end-1,3:end)-xend(2:end-1,1:end-2); 
%    dx12=xend(3:end,2:end-1)-xend(1:end-2,2:end-1);
%    dy21=yend(2:end-1,3:end)-yend(2:end-1,1:end-2);
%    dy22=yend(3:end,2:end-1)-yend(1:end-2,2:end-1);
 %   dx0 =xin(2:end-1,3:end)-xin(2:end-1,1:end-2);
%%    dy0 =yin(3:end,2:end-1)-yin(1:end-2,2:end-1);
%
%    a11=(dx11./dx0).^2.0+(dy21./dx0).*(dy21./dx0);	
%    a12=(dx11./dx0).*(dx12./dy0)+(dy21./dx0).*(dy22./dy0);
%    a21=(dx12./dy0).*(dx11./dx0)+(dy22./dy0).*(dy21./dx0);
%    a22=(dx12./dy0).*(dx12./dy0)+(dy22./dy0).^2.0;
%    lambda1=(a22+a11)/2.0+sqrt(((a22+a11)/2.0).^2.0-(a11.*a22-a12.*a21));
%    lambda2=(a22+a11)/2.0-sqrt(((a22+a11)/2.0).^2.0-(a11.*a22-a12.*a21));
%    ftleN=log(sqrt(max(lambda1,lambda2)))./abs(tIntegrate);
    
    clear u v
    disp('saving')
    fn=strcat('gibralterDay',num2str(j),'.mat');
    save(fn)
end%day loop
disp('done loop 1')

for i=1:150
    fn=strcat('gibralterDay',num2str(i),'.mat');
    load(fn,'xtr')
    xtrA(:,:,i)=xtr;
    load(fn,'ytr')
    ytrA(:,:,i)=ytr;
end
load(fn)
clear u v
fn=strcat('gibralterAll1to150.mat');
save(fn)
