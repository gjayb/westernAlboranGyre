addpath('../mStuff')
load('uvAdepth1interpolated.mat', 'XC')
load('uvAdepth1interpolated.mat', 'YC')
load('uvAdepth1interpolated.mat', 'u')
load('uvAdepth1interpolated.mat', 'v')
load('uvAdepth1interpolated.mat', 'yvel')
load('uvAdepth1interpolated.mat', 'xvel')
load('uvAdepth1interpolated.mat', 'tvel')

load('wagCore1.mat')

disp('load done')

z0=[xinm(:).';yinm(:).'];
disp('entering day loop')


for j=37;%13:12:133;%109:12:133; %9:1:30 %0:0.5:21
    j
    clear *tr*
    tmesh1=86400*(j):86400:(j+12)*86400;
    tmesh2=(j)*86400:-86400:(j-12)*86400;
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
    fnsave=strcat('wagCore1Day',num2str(j),'.mat');
    save(fnsave,'-v7.3')
end%day loop
disp('done advecting')
