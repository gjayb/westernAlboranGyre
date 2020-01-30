%1-week integration of trajectories centered on even days

load('uvAdepth1steadyNov.mat');
addpath('../mStuff')

for intLoop=1 %:1:34
    tmesh=(intLoop-1)*86400:3600:(intLoop+13)*86400;
    tmesh2=(intLoop+13)*86400:-3600:(intLoop-1)*86400;
    intLoop
    
	tvel=-1:86401:(147*86401);

    options=odeset('RelTol',10^(-11),'AbsTol',10^(-14));
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
save2pdf('novF1.pdf')
close all     
%%axis([-6 0 34 38])
    % lattrF(:,:,(intLoop-3))=lattr;
    % lontrF(:,:,(intLoop-3))=lontr;
    % clear lattr lontr
%     disp('entering integration')
%     [~,zz]=ode45(@HamEqSolver_BiLin_Irina,tmesh2,z0(:),options,u,v,xvel,yvel,tvel);%expects 1D xvel yvel tvel where u,v are defined on the meshgrid
    %[~,zz]=ode45(@HamEqSolver_BiLin_Irina,tmesh,z0(:),options,u,v,w,xvel,yvel,zvel,tvel);
%    disp('done 2')
%    xtr2=zz(:,1:2:end-1); ytr2=zz(:,2:2:end);
%    clear zz
%    lattr2=ones(size(ytr2)).*ymin+ytr2./111000;
%     lontr2=ones(size(xtr2)).*xmin+xtr2./(111000.*cosd(lattr2));
%     figure; plot(lonCoast,latCoast); hold on; plot(lontr2(:,1:5:end),lattr2(:,1:5:end)); title('backward trajectories')
%save2pdf('novB1.pdf')
%close all    
%%axis([-6 0 34 38])
   % lattrB(:,:,(intLoop-3))=lattr2;
   % lontrB(:,:,(intLoop-3))=lontr2;
%   %  clear lattr2 lontr2
end
disp('saving')
fname='/nobackup1/gbrett/trajMITgcmSteadyNovF.mat';
save(fname,'-v7.3')
disp('done all')

