%1-week integration of trajectories centered on even days

%load('uvAdepth1steadyNov.mat');
load('uvwIso29interpolated.mat');
%load('uvAdepth1interpolated.mat');
addpath('../mStuff')
clear w
[~,~,nt]=size(u);

for intLoop=115 %:1:34
    tmesh=(intLoop)*86400:7200:(intLoop+8)*86400;
    tmesh2=(intLoop)*86400:-7200:(intLoop-8)*86400;
    intLoop
    
	tvel=0:86400:((nt-1)*86400);

    options=odeset('RelTol',10^(-10),'AbsTol',10^(-13));
disp('integration 1')
    %[~,zz]=ode45(@HamEqSolver_Duffing_periodic,tmesh,z0(:),options,eps,a,nu1);
    [~,zz]=ode45(@HamEqSolver_BiLin_Irina,tmesh,z0(:),options,u,v,xvel,yvel,tvel);%expects 1D xvel yvel tvel where u,v are defined on the meshgrid
    %[~,zz]=ode45(@HamEqSolver_BiLin_Irina,tmesh,z0(:),options,u,v,w,xvel,yvel,zvel,tvel);
    xtr=zz(:,1:2:end-1);ytr=zz(:,2:2:end);
    clear zz
disp('done 1')
     
     %convert back to lon/lat
     lattr=ones(size(ytr)).*ymin+ytr./111000;
     lontr=ones(size(xtr)).*xmin+xtr./(111000.*cosd(lattr));
     figure; plot(lonCoast,latCoast); hold on; plot(lontr(:,1:5:end),lattr(:,1:5:end)); title('forward trajectories')
axis([-6 0 34 38])
save2pdf('forwardtraj8dayIso28Day115.pdf')
close all     
%save('forwardtraj8dayIso26Day115.mat','-v7.3')
%%axis([-6 0 34 38])
    % lattrF(:,:,(intLoop-3))=lattr;
    % lontrF(:,:,(intLoop-3))=lontr;
    % clear lattr lontr
     disp('entering integration')
     [~,zz]=ode45(@HamEqSolver_BiLin_Irina,tmesh2,z0(:),options,u,v,xvel,yvel,tvel);%expects 1D xvel yvel tvel where u,v are defined on the meshgrid
    %[~,zz]=ode45(@HamEqSolver_BiLin_Irina,tmesh,z0(:),options,u,v,w,xvel,yvel,zvel,tvel);
  disp('done 2')
    xtr2=zz(:,1:2:end-1); ytr2=zz(:,2:2:end);
    clear zz
    lattr2=ones(size(ytr2)).*ymin+ytr2./111000;
     lontr2=ones(size(xtr2)).*xmin+xtr2./(111000.*cosd(lattr2));
     figure; plot(lonCoast,latCoast); hold on; plot(lontr2(:,1:5:end),lattr2(:,1:5:end)); title('backward trajectories')
axis([-6 0 34 38])
save2pdf('backwardtraj8dayIso28Day115.pdf')
close all    
%%axis([-6 0 34 38])
   % lattrB(:,:,(intLoop-3))=lattr2;
   % lontrB(:,:,(intLoop-3))=lontr2;
%   %  clear lattr2 lontr2
end
disp('saving')
fname='/nobackup1/gbrett/traj8dayIso29Day115c.mat';
save(fname,'-v7.3')
disp('done all')

