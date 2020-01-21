%% advects lobe forward in time
% uses lobes found 'by hand', saves out daily positions of points

addpath('../mStuff')
load('geometrySpinupSteady.mat')
xmin=min(min(XC));
ymin=min(min(YC));

%load lobes
load('lobesByHandSigma275.mat')
lobeNames275=lobeNames;
xlobeH275=xlobeH;
ylobeH275=ylobeH;
lobeAhand275=lobeAhand;
clear xlobeHm ylobeHm *LobeH
%% load surface vel
load('uvAdepth1interpolated.mat','YC')
load('uvAdepth1interpolated.mat','u')
load('uvAdepth1interpolated.mat','v')
load('uvAdepth1interpolated.mat','xvel')
load('uvAdepth1interpolated.mat','yvel')
load('uvAdepth1interpolated.mat','tvel')
load('uvAdepth1interpolated.mat','latCoast')
load('uvAdepth1interpolated.mat','lonCoast')
%% find and advect surface lobes
% n1=find(ylobeHS(7,1,:)>0,1,'last'); %day 15
% xL016=squeeze(xlobeHS(7,1,1:n1)); 
% yL016=squeeze(ylobeHS(7,1,1:n1)); 
% xIn=(xL016-xmin).*111000.*cosd(yL016);
% yIn=(yL016-ymin).*111000;
% z0=[xIn(:)';yIn(:)'];
% timesWanted=86400*15:86400:86400*19;
%     options=odeset('RelTol',10^(-5),'AbsTol',10^(-8));
%     disp('entering integration')
%     [ttr,positions]=ode45(@HamEqSolver_BiLin_Irina,timesWanted,z0(:),options,u,v,xvel,yvel,tvel);
% 
%     xtr=positions(:,1:2:end-1);
%     ytr=positions(:,2:2:end);
%     disp('integration 1 done')
% lattrL0S=ones(size(ytr)).*ymin+ytr./111000;
% lontrL0S=ones(size(xtr)).*xmin+xtr./(111000.*cosd(lattrL0S));
% %b
% n1=find(ylobeHS(7,3,:)>0,1,'last'); %day 15
% xB=squeeze(xlobeHS(7,3,1:n1)); 
% yB=squeeze(ylobeHS(7,3,1:n1)); 
% xIn=(xB-xmin).*111000.*cosd(yB);
% yIn=(yB-ymin).*111000;
% z0=[xIn(:)';yIn(:)'];
% timesWanted=86400*15:86400:86400*19;
%     options=odeset('RelTol',10^(-5),'AbsTol',10^(-8));
%     disp('entering integration')
%     [ttr,positions]=ode45(@HamEqSolver_BiLin_Irina,timesWanted,z0(:),options,u,v,xvel,yvel,tvel);
% 
%     xtr=positions(:,1:2:end-1);
%     ytr=positions(:,2:2:end);
%     disp('integration 1 done')
% lattrBS=ones(size(ytr)).*ymin+ytr./111000;
% lontrBS=ones(size(xtr)).*xmin+xtr./(111000.*cosd(lattrBS));
% %c
% n1=find(ylobeHS(9,4,:)>0,1,'last'); %day 17
% xC=squeeze(xlobeHS(9,4,1:n1)); 
% yC=squeeze(ylobeHS(9,4,1:n1)); 
% xIn=(xC-xmin).*111000.*cosd(yC);
% yIn=(yC-ymin).*111000;
% z0=[xIn(:)';yIn(:)'];
% timesWanted=86400*17:86400:86400*20;
%     options=odeset('RelTol',10^(-5),'AbsTol',10^(-8));
%     disp('entering integration')
%     [ttr,positions]=ode45(@HamEqSolver_BiLin_Irina,timesWanted,z0(:),options,u,v,xvel,yvel,tvel);
% 
%     xtr=positions(:,1:2:end-1);
%     ytr=positions(:,2:2:end);
%     disp('integration 1 done')
% lattrCS=ones(size(ytr)).*ymin+ytr./111000;
% lontrCS=ones(size(xtr)).*xmin+xtr./(111000.*cosd(lattrCS));
% %L1
% n1=find(ylobeHS(12,1,:)>0,1,'last'); %day 20
% xL1=squeeze(xlobeHS(12,1,1:n1)); 
% yL1=squeeze(ylobeHS(12,1,1:n1)); 
% xIn=(xL1-xmin).*111000.*cosd(yL1);
% yIn=(yL1-ymin).*111000;
% z0=[xIn(:)';yIn(:)'];
% timesWanted=86400*20:86400:86400*24;
%     options=odeset('RelTol',10^(-5),'AbsTol',10^(-8));
%     disp('entering integration')
%     [ttr,positions]=ode45(@HamEqSolver_BiLin_Irina,timesWanted,z0(:),options,u,v,xvel,yvel,tvel);
% 
%     xtr=positions(:,1:2:end-1);
%     ytr=positions(:,2:2:end);
%     disp('integration last S done')
% lattrL1S=ones(size(ytr)).*ymin+ytr./111000;
% lontrL1S=ones(size(xtr)).*xmin+xtr./(111000.*cosd(lattrL1S));
% clear xtr ytr
% %save
% fn='surfaceAdection4lobes.mat'
% save(fn,'x*','y*','*tr*')
% %%
% clear u v
% load('uvwIso265interpolated.mat', 'u');
% load('uvwIso265interpolated.mat', 'v');
% %% sigma=26.5 advection
% % find and advect lobes
% n1=find(ylobeH265(15,2,:)>0,1,'last'); %day 15
% xL016=squeeze(xlobeH265(15,2,1:n1)); 
% yL016=squeeze(ylobeH265(15,2,1:n1)); 
% xIn=(xL016-xmin).*111000.*cosd(yL016);
% yIn=(yL016-ymin).*111000;
% z0=[xIn(:)';yIn(:)'];
% timesWanted=86400*15:86400:86400*19;
%     options=odeset('RelTol',10^(-5),'AbsTol',10^(-8));
%     disp('entering integration')
%     [ttr,positions]=ode45(@HamEqSolver_BiLin_Irina,timesWanted,z0(:),options,u,v,xvel,yvel,tvel);
% 
%     xtr=positions(:,1:2:end-1);
%     ytr=positions(:,2:2:end);
%     disp('integration 1 26.5 done')
% lattrL0265=ones(size(ytr)).*ymin+ytr./111000;
% lontrL0265=ones(size(xtr)).*xmin+xtr./(111000.*cosd(lattrL0265));
% %b
% n1=find(ylobeH265(15,4,:)>0,1,'last'); %day 15
% xB=squeeze(xlobeH265(15,4,1:n1)); 
% yB=squeeze(ylobeH265(15,4,1:n1)); 
% xIn=(xB-xmin).*111000.*cosd(yB);
% yIn=(yB-ymin).*111000;
% z0=[xIn(:)';yIn(:)'];
% timesWanted=86400*15:86400:86400*19;
%     options=odeset('RelTol',10^(-5),'AbsTol',10^(-8));
%     disp('entering integration')
%     [ttr,positions]=ode45(@HamEqSolver_BiLin_Irina,timesWanted,z0(:),options,u,v,xvel,yvel,tvel);
% 
%     xtr=positions(:,1:2:end-1);
%     ytr=positions(:,2:2:end);
%     disp('integration 1 done')
% lattrB265=ones(size(ytr)).*ymin+ytr./111000;
% lontrB265=ones(size(xtr)).*xmin+xtr./(111000.*cosd(lattrB265));
% %c
% n1=find(ylobeH265(17,3,:)>0,1,'last'); %day 17
% xC=squeeze(xlobeH265(17,3,1:n1)); 
% yC=squeeze(ylobeH265(17,3,1:n1)); 
% xIn=(xC-xmin).*111000.*cosd(yC);
% yIn=(yC-ymin).*111000;
% z0=[xIn(:)';yIn(:)'];
% timesWanted=86400*17:86400:86400*20;
%     options=odeset('RelTol',10^(-5),'AbsTol',10^(-8));
%     disp('entering integration')
%     [ttr,positions]=ode45(@HamEqSolver_BiLin_Irina,timesWanted,z0(:),options,u,v,xvel,yvel,tvel);
% 
%     xtr=positions(:,1:2:end-1);
%     ytr=positions(:,2:2:end);
%     disp('integration 1 done')
% lattrC265=ones(size(ytr)).*ymin+ytr./111000;
% lontrC265=ones(size(xtr)).*xmin+xtr./(111000.*cosd(lattrC265));
% %L1
% n1=find(ylobeH265(20,1,:)>0,1,'last'); %day 20
% xL1=squeeze(xlobeH265(20,1,1:n1)); 
% yL1=squeeze(ylobeH265(20,1,1:n1)); 
% xIn=(xL1-xmin).*111000.*cosd(yL1);
% yIn=(yL1-ymin).*111000;
% z0=[xIn(:)';yIn(:)'];
% timesWanted=86400*20:86400:86400*24;
%     options=odeset('RelTol',10^(-5),'AbsTol',10^(-8));
%     disp('entering integration')
%     [ttr,positions]=ode45(@HamEqSolver_BiLin_Irina,timesWanted,z0(:),options,u,v,xvel,yvel,tvel);
% 
%     xtr=positions(:,1:2:end-1);
%     ytr=positions(:,2:2:end);
%     disp('integration last 26.5 done')
% lattrL1265=ones(size(ytr)).*ymin+ytr./111000;
% lontrL1265=ones(size(xtr)).*xmin+xtr./(111000.*cosd(lattrL1265));
% clear xtr ytr
% %save
% fn='sigma265Adection4lobes.mat'
% save(fn,'x*','y*','*tr*')
%%
clear u v
load('uvwIso27interpolated.mat', 'u');
load('uvwIso27interpolated.mat', 'v');
%% sigma=27 advection
% find and advect lobes
n1=find(ylobeH27(15,2,:)>0,1,'last'); %day 15
xL016=squeeze(xlobeH27(15,2,1:n1)); 
yL016=squeeze(ylobeH27(15,2,1:n1)); 
xIn=(xL016-xmin).*111000.*cosd(yL016);
yIn=(yL016-ymin).*111000;
z0=[xIn(:)';yIn(:)'];
timesWanted=86400*15:86400:86400*19;
    options=odeset('RelTol',10^(-5),'AbsTol',10^(-8));
    disp('entering integration')
    [ttr,positions]=ode45(@HamEqSolver_BiLin_Irina,timesWanted,z0(:),options,u,v,xvel,yvel,tvel);

    xtr=positions(:,1:2:end-1);
    ytr=positions(:,2:2:end);
    disp('integration 1 26.5 done')
lattrL027=ones(size(ytr)).*ymin+ytr./111000;
lontrL027=ones(size(xtr)).*xmin+xtr./(111000.*cosd(lattrL027));
%b
n1=find(ylobeH27(15,4,:)>0,1,'last'); %day 15
xB=squeeze(xlobeH27(15,4,1:n1)); 
yB=squeeze(ylobeH27(15,4,1:n1)); 
xIn=(xB-xmin).*111000.*cosd(yB);
yIn=(yB-ymin).*111000;
z0=[xIn(:)';yIn(:)'];
timesWanted=86400*15:86400:86400*19;
    options=odeset('RelTol',10^(-5),'AbsTol',10^(-8));
    disp('entering integration')
    [ttr,positions]=ode45(@HamEqSolver_BiLin_Irina,timesWanted,z0(:),options,u,v,xvel,yvel,tvel);

    xtr=positions(:,1:2:end-1);
    ytr=positions(:,2:2:end);
    disp('integration 1 done')
lattrB27=ones(size(ytr)).*ymin+ytr./111000;
lontrB27=ones(size(xtr)).*xmin+xtr./(111000.*cosd(lattrB27));
%c
n1=find(ylobeH27(17,3,:)>0,1,'last'); %day 17
xC=squeeze(xlobeH27(17,3,1:n1)); 
yC=squeeze(ylobeH27(17,3,1:n1)); 
xIn=(xC-xmin).*111000.*cosd(yC);
yIn=(yC-ymin).*111000;
z0=[xIn(:)';yIn(:)'];
timesWanted=86400*17:86400:86400*20;
    options=odeset('RelTol',10^(-5),'AbsTol',10^(-8));
    disp('entering integration')
    [ttr,positions]=ode45(@HamEqSolver_BiLin_Irina,timesWanted,z0(:),options,u,v,xvel,yvel,tvel);

    xtr=positions(:,1:2:end-1);
    ytr=positions(:,2:2:end);
    disp('integration 1 done')
lattrC27=ones(size(ytr)).*ymin+ytr./111000;
lontrC27=ones(size(xtr)).*xmin+xtr./(111000.*cosd(lattrC27));
%L1
n1=find(ylobeH27(20,1,:)>0,1,'last'); %day 20
xL1=squeeze(xlobeH27(20,1,1:n1)); 
yL1=squeeze(ylobeH27(20,1,1:n1)); 
xIn=(xL1-xmin).*111000.*cosd(yL1);
yIn=(yL1-ymin).*111000;
z0=[xIn(:)';yIn(:)'];
timesWanted=86400*20:86400:86400*24;
    options=odeset('RelTol',10^(-5),'AbsTol',10^(-8));
    disp('entering integration')
    [ttr,positions]=ode45(@HamEqSolver_BiLin_Irina,timesWanted,z0(:),options,u,v,xvel,yvel,tvel);

    xtr=positions(:,1:2:end-1);
    ytr=positions(:,2:2:end);
    disp('integration last 26.5 done')
lattrL127=ones(size(ytr)).*ymin+ytr./111000;
lontrL127=ones(size(xtr)).*xmin+xtr./(111000.*cosd(lattrL127));
clear xtr ytr
%save
fn='sigma27Adection4lobes.mat'
save(fn,'x*','y*','*tr*')
