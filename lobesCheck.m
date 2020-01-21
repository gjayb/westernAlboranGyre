%uses output from lobes2
addpath('../mStuff')
load('lobes8day148DailyVSurfaceV2.mat')
%load('uvAdepth1interpolated.mat', 'u')
%load('uvAdepth1interpolated.mat', 'v')
%load('uvAdepth1interpolated.mat', 'yvel')
%load('uvAdepth1interpolated.mat', 'xvel')
%load('uvAdepth1interpolated.mat', 'tvel')

load('uvwIso27interpolated.mat')

disp('load done')

tvel=0:86400:147*86400;

lobedayi=[5 13 46 49];
lobeday=lobedayi+8;

for i=lobeday
    i
    xini=XC(inside(:,:,i-8)==1); 
    yini=YC(inside(:,:,i-8)==1); 
    xinM1=(xini-min(min(XC))*ones(size(xini))).*111000.*cosd(yini); 
    yinM1=(yini-min(min(YC))*ones(size(yini))).*111000;
    z0=[xinM1(:)';yinM1(:)']; 
    
    tmesh1=86400*(i):86400:(i+12)*86400;
    tmesh2=(i)*86400:-86400:(i-12)*86400;
    options=odeset('RelTol',10^(-10),'AbsTol',10^(-13));

    disp('int forward')
    [~,zz]=ode45(@HamEqSolver_BiLin_Irina,tmesh1,z0(:),options,u,v,xvel,yvel,tvel);%expects 1D xvel yvel tvel where u,v are defined on the meshgrid
    xtr=zz(:,1:2:end-1);
    ytr=zz(:,2:2:end);
    
     lattrF=ones(size(ytr)).*min(min(YC))+ytr./111000;
     lontrF=ones(size(xtr)).*min(min(XC))+xtr./(111000.*cosd(lattrF));
    disp('int backward')
    [~,zz]=ode45(@HamEqSolver_BiLin_Irina,tmesh2,z0(:),options,u,v,xvel,yvel,tvel);%expects 1D xvel yvel tvel where u,v are defined on the meshgrid
    xtr2=zz(:,1:2:end-1);
    ytr2=zz(:,2:2:end);
    
     lattrB=ones(size(ytr2)).*min(min(YC))+ytr2./111000;
     lontrB=ones(size(xtr2)).*min(min(XC))+xtr2./(111000.*cosd(lattrB));

disp('saving')
fn=strcat('lobeIso27Day',num2str(i),'int12.mat');
save(fn,'-v7.3')
end
disp('done')
