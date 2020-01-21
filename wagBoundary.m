
addpath('../mStuff')
load('uvHdepth1interpolated.mat')

circleCenters=zeros(2,2);
    %148 steady for WAG
    circleCenters(1:2,:)=[-5.35,35.86;-3.05,35.4];
    
r=8/110;
    ang=0:0.01:2*pi;
    yin=[]; xin=[];
   % for i=1:2
    xc1=circleCenters(1,1)+r.*cos(ang);
    yc1=circleCenters(1,2)+r.*sin(ang);
    %yin=cat(1,yin,yc);
    %xin=cat(1,xin,xc);
   xc2=circleCenters(2,1)+r.*cos(ang);
    yc2=circleCenters(2,2)+r.*sin(ang);
    % end
% r=2/110;
% for i=1:2
%     xc=circleCenters(i,1)+r.*cos(ang);
%     yc=circleCenters(i,2)+r.*sin(ang);
%     yin=cat(1,yin,yc);
%     xin=cat(1,xin,xc);
% end

%xin=xin(:);
%    yin=yin(:);
    %xmin=min(xin);
    %ymin=min(yin);
    xinM1=(xc1-min(min(XC))*ones(size(xc1))).*111000.*cosd(yc1); 
    yinM1=(yc1-min(min(YC))*ones(size(yc1))).*111000;
    xinM2=(xc2-min(min(XC))*ones(size(xc2))).*111000.*cosd(yc2); 
    yinM2=(yc2-min(min(YC))*ones(size(yc2))).*111000;
    %z0=[xinM(:).';yinM(:).'];
    xmin=min(min(XC));
    ymin=min(min(YC));
    
    z0f=[xinM1(:)';yinM1(:)'];
    z0b=[xinM2(:)';yinM2(:)'];
    
%%
disp('entering day loop')

for j=10:12 %0:0.5:21
    j
    clear *tr*
    tmesh=86400*(j-8):1440:(j-6)*86400;
    tmesh2=(j+8)*86400:-14400:(j+6)*86400;
    options=odeset('RelTol',10^(-10),'AbsTol',10^(-13));
    [~,zz]=ode45(@HamEqSolver_BiLin_Irina,tmesh,z0f(:),options,u,v,xvel,yvel,tvel);%expects 1D xvel yvel tvel where u,v are defined on the meshgrid
    xtrA=zz(:,1:2:end-1);ytrA=zz(:,2:2:end);
    clear zz
     
     %convert back to lon/lat
     lattrA=ones(size(ytrA)).*ymin+ytrA./111000;
     lontrA=ones(size(xtrA)).*xmin+xtrA./(111000.*cosd(lattrA));
     %figure; plot(lonCoast,latCoast); hold on; plot(lontr(:,1:5:end),lattr(:,1:5:end)); title('forward trajectories')
     %axis([-6 0 34 38])
     
     disp('entering integration 2a')
     [~,zz]=ode45(@HamEqSolver_BiLin_Irina,tmesh2,z0b(:),options,u,v,xvel,yvel,tvel);%expects 1D xvel yvel tvel where u,v are defined on the meshgrid
    sizezz=size(zz)
     xtr2A=zz(:,1:2:end-1); ytr2A=zz(:,2:2:end);
    clear zz
    lattr2A=ones(size(ytr2A)).*ymin+ytr2A./111000;
     lontr2A=ones(size(xtr2A)).*xmin+xtr2A./(111000.*cosd(lattr2A));
    
     disp('save A')
%   fn=strcat('wagBoundaryday',num2str(j),'A2day.mat');
%   save(fn,'-v7.3')
    figure
    plot(lonCoast,latCoast,'k')
    hold on
    plot(lontrA(end,1:2:end),lattrA(end,1:2:end),'ro')
    plot(lontr2A(end,2:2:end),lattr2A(end,2:2:end),'bo')
    title(strcat('Manifolds day 2 of 8 day integration day',num2str(j)))
    fn=strcat('wagBoundaryday',num2str(j),'A2dayH.pdf');
    save2pdf(fn)
    
    %%
    tmesh=86400*(j-6):14400:(j-4)*86400;
    tmesh2=(j+6)*86400:-14400:(j+4)*86400;
    %options=odeset('RelTol',10^(-10),'AbsTol',10^(-13));
    
    xin=xtrA(end,:); yin=ytrA(end,:);
    xinb=0.5.*xin(2:end)+0.5.*xin(1:end-1);
    yinb=0.5.*yin(2:end)+0.5.*yin(1:end-1);
    %xinbb=0.5.*xin(2:2:end-2)+0.5.*xin(4:2:end);
    %yinbb=0.5.*yin(2:2:end-2)+0.5.*yin(4:2:end);
    %xinB=cat(1,xin.',xinb.');%,xinbb.');
    %yinB=cat(1,yin.',yinb.');%,yinbb.');
    
    xin2=xtr2A(end,:); yin2=ytr2A(end,:);
    xinb2=0.5.*xin2(2:end)+0.5.*xin2(1:end-1);
    yinb2=0.5.*yin2(2:end)+0.5.*yin2(1:end-1);
    
    xinB=zeros(size(xin)+size(xinb));
    yinB=xinB;
    yin2B=xinB;
    xin2B=xinB;
    for i=1:length(xinb)
       xinB(2*i-1)= xin(i);
       xinB(2*i)= xinb(i);
       yinB(2*i-1)= yin(i);
       yinB(2*i)= yinb(i);
       xin2B(2*i-1)= xin2(i);
       xin2B(2*i)= xinb2(i);
       yin2B(2*i-1)= yin2(i);
       yin2B(2*i)= yinb2(i);
    end
    xinB(end)=xin(end);
    yinB(end)=yin(end);
    z0a=[xinB(:)';yinB(:)'];
    xin2B(end)=xin2(end);
    yin2B(end)=yin2(end);
    z02a=[xin2B(:)';yin2B(:)'];
    
    disp('entering integration 1b')
    [~,zz]=ode45(@HamEqSolver_BiLin_Irina,tmesh,z0a(:),options,u,v,xvel,yvel,tvel);%expects 1D xvel yvel tvel where u,v are defined on the meshgrid
    xtrB=zz(:,1:2:end-1);ytrB=zz(:,2:2:end);
    clear zz
     
     %convert back to lon/lat
     lattrB=ones(size(ytrB)).*ymin+ytrB./111000;
     lontrB=ones(size(xtrB)).*xmin+xtrB./(111000.*cosd(lattrB));
     %figure; plot(lonCoast,latCoast); hold on; plot(lontr(:,1:5:end),lattr(:,1:5:end)); title('forward trajectories')
     %axis([-6 0 34 38])
     
     disp('entering integration 2b')
     [~,zz]=ode45(@HamEqSolver_BiLin_Irina,tmesh2,z02a(:),options,u,v,xvel,yvel,tvel);%expects 1D xvel yvel tvel where u,v are defined on the meshgrid
    xtr2B=zz(:,1:2:end-1); ytr2B=zz(:,2:2:end);
    clear zz
    lattr2B=ones(size(ytr2B)).*ymin+ytr2B./111000;
     lontr2B=ones(size(xtr2B)).*xmin+xtr2B./(111000.*cosd(lattr2B));
    
   %fn=strcat('wagBoundaryday',num2str(j),'B2day.mat');
   %save(fn,'-v7.3')
    figure; plot(lontrB(end,:),lattrB(end,:),'ro')
    hold on; plot(lontr2B(end,:),lattr2B(end,:),'bo')
    hold on; plot(lontr2B(end,:),lattr2B(end,:),'bo')
    hold on; plot(lonCoast,latCoast,'k')
    title(strcat('Manifolds day 4 of 8 day integration day',num2str(j)))
    fn=strcat('wagBoundaryday',num2str(j),'B2dayH.pdf');
    save2pdf(fn)
    %%
        tmesh=86400*(j-4):14400:(j-2)*86400;
    tmesh2=(j+4)*86400:-14400:(j+2)*86400;
%     %options=odeset('RelTol',10^(-10),'AbsTol',10^(-13));
%     
        xin=xtrB(end,:); yin=ytrB(end,:);
        xin2=xtr2B(end,:); yin2=ytr2B(end,:);
        xinc=0.5.*xin(2:end)+0.5.*xin(1:end-1);
        yinc=0.5.*yin(2:end)+0.5.*yin(1:end-1);
        xinc2=0.5.*xin2(2:end)+0.5.*xin2(1:end-1);
        yinc2=0.5.*yin2(2:end)+0.5.*yin2(1:end-1);
        
    xinC=zeros(size(xin)+size(xinc));
    yinC=xinC;
    yin2C=xinC;
    xin2C=xinC;
    for i=1:length(xinc)
       xinC(2*i-1)= xin(i);
       xinC(2*i)= xinc(i);
       yinC(2*i-1)= yin(i);
       yinC(2*i)= yinc(i);
       xin2C(2*i-1)= xin2(i);
       xin2C(2*i)= xinc2(i);
       yin2C(2*i-1)= yin2(i);
       yin2C(2*i)= yinc2(i);
    end
    xinC(end)=xin(end);
    yinC(end)=yin(end);
    z0bb=[xinC(:)';yinC(:)'];
    xin2C(end)=xin2(end);
    yin2C(end)=yin2(end);
    z02b=[xin2C(:)';yin2C(:)'];
%     xinc=0.5.*xin(2:end)+0.5.*xin(1:end-1);
%     yinc=0.5.*yin(2:end)+0.5.*yin(1:end-1);
%     xinC=cat(1,xin,xinc);
%     yinC=cat(1,yin,yinc);
%     z0b=[xinC(:)';yinC(:)'];
%     

%     xin2c=0.5.*xin2(2:end)+0.5.*xin2(1:end-1);
%     yin2c=0.5.*yin2(2:end)+0.5.*yin2(1:end-1);
%     xin2C=cat(1,xin,xin2c);
%     yin2C=cat(1,yinC,yin2c);
%     z02b=[xin2C(:)';yin2C(:)'];
%     
    disp('entering integration 1c')
    [~,zz]=ode45(@HamEqSolver_BiLin_Irina,tmesh,z0bb(:),options,u,v,xvel,yvel,tvel);%expects 1D xvel yvel tvel where u,v are defined on the meshgrid
    xtrC=zz(:,1:2:end-1);ytrC=zz(:,2:2:end);
    clear zz
     
     %convert back to lon/lat
     lattrC=ones(size(ytrC)).*ymin+ytrC./111000;
     lontrC=ones(size(xtrC)).*xmin+xtrC./(111000.*cosd(lattrC));
     %figure; plot(lonCoast,latCoast); hold on; plot(lontr(:,1:5:end),lattr(:,1:5:end)); title('forward trajectories')
     %axis([-6 0 34 38])
     
     disp('entering integration 2c')
     [~,zz]=ode45(@HamEqSolver_BiLin_Irina,tmesh2,z02b(:),options,u,v,xvel,yvel,tvel);%expects 1D xvel yvel tvel where u,v are defined on the meshgrid
    xtr2C=zz(:,1:2:end-1); ytr2C=zz(:,2:2:end);
    clear zz
    lattr2C=ones(size(ytr2C)).*ymin+ytr2C./111000;
     lontr2C=ones(size(xtr2C)).*xmin+xtr2C./(111000.*cosd(lattr2C));
    
   %fn=strcat('wagBoundaryday',num2str(j),'C2day.mat');
   %save(fn,'-v7.3')
    figure
    plot(lonCoast,latCoast,'k')
    hold on
    plot(lontrC(end,:),lattrC(end,:),'ro')
    plot(lontr2C(end,:),lattr2C(end,:),'bo')
    title(strcat('Manifolds day 6 of 8 day integration day',num2str(j)))
    fn=strcat('wagBoundaryday',num2str(j),'C2dayH.pdf');
    save2pdf(fn)
    
 %%   
        tmesh=86400*(j-2):14400:(j)*86400;
    tmesh2=(j+2)*86400:-14400:(j)*86400;
%     %options=odeset('RelTol',10^(-10),'AbsTol',10^(-13));
%     
        xin=xtrC(end,:); yin=ytrC(end,:);
        xin2=xtr2C(end,:); yin2=ytr2C(end,:);
        xind=0.5.*xin(2:end)+0.5.*xin(1:end-1);
        yind=0.5.*yin(2:end)+0.5.*yin(1:end-1);
        xind2=0.5.*xin2(2:end)+0.5.*xin2(1:end-1);
        yind2=0.5.*yin2(2:end)+0.5.*yin2(1:end-1);
        
    xinD=zeros(size(xin)+size(xind));
    yinD=xinD;
    yin2D=xinD;
    xin2D=xinD;
    for i=1:length(xind)
       xinD(2*i-1)= xin(i);
       xinD(2*i)= xind(i);
       yinD(2*i-1)= yin(i);
       yinD(2*i)= yind(i);
       xin2D(2*i-1)= xin2(i);
       xin2D(2*i)= xind2(i);
       yin2D(2*i-1)= yin2(i);
       yin2D(2*i)= yind2(i);
    end
    xinD(end)=xin(end);
    yinD(end)=yin(end);
    z0c=[xinC(:)';yinC(:)'];
    xin2D(end)=xin2(end);
    yin2D(end)=yin2(end);
    z02c=[xin2D(:)';yin2D(:)'];
    
    
    disp('entering integration 1d')
    [~,zz]=ode45(@HamEqSolver_BiLin_Irina,tmesh,z0c(:),options,u,v,xvel,yvel,tvel);%expects 1D xvel yvel tvel where u,v are defined on the meshgrid
    xtrD=zz(:,1:2:end-1);ytrD=zz(:,2:2:end);
    clear zz
     
     %convert back to lon/lat
     lattrD=ones(size(ytrD)).*ymin+ytrD./111000;
     lontrD=ones(size(xtrD)).*xmin+xtrD./(111000.*cosd(lattrD));
     %figure; plot(lonCoast,latCoast); hold on; plot(lontr(:,1:5:end),lattr(:,1:5:end)); title('forward trajectories')
     %axis([-6 0 34 38])
     
     disp('entering integration 2d')
     [~,zz]=ode45(@HamEqSolver_BiLin_Irina,tmesh2,z02c(:),options,u,v,xvel,yvel,tvel);%expects 1D xvel yvel tvel where u,v are defined on the meshgrid
    xtr2D=zz(:,1:2:end-1); ytr2D=zz(:,2:2:end);
    clear zz
    lattr2D=ones(size(ytr2C)).*ymin+ytr2C./111000;
     lontr2D=ones(size(xtr2C)).*xmin+xtr2C./(111000.*cosd(lattr2C));
    
   fn=strcat('wagBoundaryday',num2str(j),'D2dayH.mat');
   save(fn,'-v7.3')
    figure
    plot(lonCoast,latCoast,'k')
    hold on
    plot(lontrD(end,:),lattrD(end,:),'ro')
    plot(lontr2D(end,:),lattr2D(end,:),'bo')
    title(strcat('Manifolds 8 day integration day',num2str(j)))
    fn=strcat('wagBoundaryday',num2str(j),'D2dayH.pdf');
    save2pdf(fn)
    
   disp('done j')
end
disp('done')
