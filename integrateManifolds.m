%drawing manifolds

%load('C:\Users\JayB\Documents\MATLAB\MITgcm\integrated1weekSurfaceSteady.mat')

xcM=(XC-(min(XC(:)))*ones(size(XC))).*111000.*cosd(YC); 
ycM=(YC-(min(YC(:)))*ones(size(YC))).*111000;
xv=min(min(xcM)):1000:max(max(xcM));
yv=min(min(ycM)):1000:max(max(ycM));
[xvel,yvel]=meshgrid(xv,yv);
latg=ones(size(yvel)).*min(min(YC))+yvel./111000;
long=ones(size(xvel)).*min(min(XC))+xvel./(111000.*cosd(latg));
% Urot200s=mean(Urot200,3);
% Vrot200s=mean(Vrot200,3);
% Ug=griddata(xcM,ycM,Urot200s(:,:,1),xvel,yvel);
% Vg=griddata(xcM,ycM,Vrot200s(:,:,1),xvel,yvel);
% Ug(isnan(Ug))=0;
% Vg(isnan(Vg))=0;
% figure; quiver(long,latg,Ug,Vg)
% hold on; plot(lonCoast,latCoast,'k')
% grid on
% ke=Ug.^2+Vg.^2;
% figure; pcolorJ(long,latg,ke); shading 'flat'; colorbar

circleCenters=zeros(20,2);
%surface:
%circleCenters(1:10,:)=[-5.33000000000000,36.1300000000000;-5.3000000000000,35.9500000000000;-2.98000000000000,35.4400000000000;-3.4100,36.0800;-4.28,35.7500000;-4.55000000000000,36.5300000000000;-2.12000000000000,36.7200000000000;-1.2000000000000,35.5800000000000;-1.33000000000000,36.4300000000000;-2.100, 35.8800];
%100m depth:
%circleCenters(1:11,:)=[-5.33000000000000,36.130000000000;-5.27000000000000,35.9200000000000;-2.96000000000000,35.4500000000000;-3.200,35.900;-4.3000000000000,35.7500000000000;-4.6000000000000,36.4500000000000;-3.6250,35.3500;-2.12000000000000,36.7200000000000;-1.25000000000000,35.650000000000;-1.5000000000,36.60000000000;-2.100, 35.8600];
%200m depth
circleCenters(1:13,:)=[-5.33000000000000,36.090000000000;-5.25000000000000,35.900000000000;-2.94000000000000,35.5200000000000;-4.20,36.10;-4.17000000000,35.800000000000;-4.8000000000000,35.500000000000;-3.550,35.400;-2.95, 36.1;-2.1000000000000,36.4900000000;-0.6,36.0000000000;-2.0700, 36.150; -1.44,36.92; -1.4, 36.3];

r=5/110;
ang=0:0.1:2*pi;
yin=[]; xin=[];
for i=1:13
xc=circleCenters(i,1)+r.*cos(ang);
yc=circleCenters(i,2)+r.*sin(ang);
yin=cat(1,yin,yc);
xin=cat(1,xin,xc);
end
xin=xin(:);
yin=yin(:);
xmin=min(xin);
ymin=min(yin);
xinM=(xin-min(min(XC))*ones(size(xin))).*111000.*cosd(yin); 
yinM=(yin-min(min(YC))*ones(size(yin))).*111000;
qin=[xinM(:).';yinM(:).'];
xmin=min(min(XC));
ymin=min(min(YC));



uInt=repmat(Ug,[1 1 2]);
vInt=repmat(Vg,[1 1 2]);

i=1;
tmesh=((i-1)*86400):3600:((i+69)*86400); % times for trajectory integration in seconds; 1h=3600
    tmesh2=((i+69)*86400):-3600:((i-1)*86400);
    options=odeset('RelTol',10^(-6),'AbsTol',10^(-9));
    %[~,zz]=ode45(@HamEqSolver_Duffing_periodic,tmesh,z0(:),options,eps,a,nu1);
    [~,zz]=ode45(@HamEqSolver_BiLin_Irina,tmesh,qin(:),options,uInt,vInt,xv,yv,[min(tmesh) max(tmesh)+1]);%expects 1D xvel yvel tvel where u,v are defined on the meshgrid
    %[~,zz]=ode45(@HamEqSolver_BiLin_Irina,tmesh,z0(:),options,u,v,w,xvel,yvel,zvel,tvel);
    xtr=zz(:,1:2:end-1);ytr=zz(:,2:2:end);
    clear zz
     
     %convert back to lon/lat
     lattr=ones(size(ytr)).*ymin+ytr./111000;
     lontr=ones(size(xtr)).*xmin+xtr./(111000.*cosd(lattr));
     
     disp('entering integration')
     [~,zz]=ode45(@HamEqSolver_BiLin_Irina,tmesh2,qin(:),options,uInt,vInt,xv,yv,[min(tmesh) max(tmesh)+1]);%expects 1D xvel yvel tvel where u,v are defined on the meshgrid
    %[~,zz]=ode45(@HamEqSolver_BiLin_Irina,tmesh,z0(:),options,u,v,w,xvel,yvel,zvel,tvel);
    disp('done')
    xtr2=zz(:,1:2:end-1); ytr2=zz(:,2:2:end);
    clear zz
    lattr2=ones(size(ytr2)).*ymin+ytr2./111000;
     lontr2=ones(size(xtr2)).*xmin+xtr2./(111000.*cosd(lattr2));
     figure; plot(lonCoast,latCoast,'k'); hold on; plot(lontr2,lattr2,'b'); plot(lontr,lattr,'r'); plot(xin,yin,'ko'); title('forward and backward trajectories')
    axis([-6 0 34 38])

    figure; plot(lonCoast,latCoast,'k'); hold on; plot(lontr,lattr,'r'); plot(lontr2,lattr2,'b'); plot(xin,yin,'ko'); title('forward and backward trajectories')
    axis([-6 0 34 38])