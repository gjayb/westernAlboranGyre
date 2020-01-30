
%load('geometrySpinupSteady.mat');

% uh=rdmds('U',NaN);
% uh=squeeze(uh(:,:,1,1:48));
% vh=rdmds('V',NaN);
% vh=squeeze(vh(:,:,1,1:48));
% uhr =uh.*repmat(AngleCS,[1 1 48]) - vh.*repmat(AngleSN,[1 1 48]);  
% vhr =uh.*repmat(AngleSN,[1 1 48]) + vh.*repmat(AngleCS,[1 1 48]);
% uhv=repmat(uhr,[1 1 6]);
% vhv=repmat(vhr,[1 1 6]);
 tvelh=0:3600:11*86400;
% 
% ud=rdmds('Uave',NaN);
% vd=rdmds('Vave',NaN);
% ud=squeeze(ud(:,:,1,:));
% vd=squeeze(vd(:,:,1,:));
% udr =ud.*repmat(AngleCS,[1 1 2]) - vd.*repmat(AngleSN,[1 1 2]);  
% vdr =ud.*repmat(AngleSN,[1 1 2]) + vd.*repmat(AngleCS,[1 1 2]);
% udv=repmat(udr,[1 1 6]);
% vdv=repmat(vdr,[1 1 6]);
 tveld=0:86400:11*86400;
%load('dailyhourly.mat');


% NX=120; NY=80;
%     xmin=-5.5;
%     xmax=1;
%     ymin=35;
%     ymax=37;
% 
%     xin=linspace(xmin,xmax,NX);
%     yin=linspace(ymin,ymax,NY);
% 
%     [xin2, yin2]=meshgrid(xin,yin);
%     %convert to meters from lat/lon; lat is y, lon is x
%     xinM=(xin2-xmin*ones(size(xin2))).*111000.*cosd(yin2); 
%     yinM=(yin2-ymin*ones(size(yin2))).*111000;
%     
% xvel2=(XC-xmin*ones(size(XC))).*111000.*cosd(YC); 
% yvel2=(YC-ymin*ones(size(YC))).*111000;
% 
% xvel=min(min(xvel2)):1000:max(max(xvel2));
% yvel=min(min(yvel2)):1000:max(max(yvel2));
% 
% [xvelg,yvelg]=meshgrid(xvel,yvel);
% for k=1:length(tveld)
%     udvel(:,:,k)=griddata(xvel2,yvel2,udv(:,:,k),xvelg,yvelg);
%     vdvel(:,:,k)=griddata(xvel2,yvel2,vdv(:,:,k),xvelg,yvelg);
% end
% 
% for k=1:length(tvelh)
%     uhvel(:,:,k)=griddata(xvel2,yvel2,uhv(:,:,k),xvelg,yvelg);
%     vhvel(:,:,k)=griddata(xvel2,yvel2,vhv(:,:,k),xvelg,yvelg);
% end
% 
% z0=[xinM(:)';yinM(:)'];

%load('dailyhourlyGridded.mat');

disp('entering integration')

tmesh=0:1800:8*86400-1;
    tmesh2=8*86400-1:-1800:0;
    options=odeset('RelTol',10^(-7),'AbsTol',10^(-10));
    %[~,zz]=ode45(@HamEqSolver_Duffing_periodic,tmesh,z0(:),options,eps,a,nu1);
    [~,zz]=ode45(@HamEqSolver_BiLin_Irina,tmesh,z0(:),options,udvel,vdvel,xvel,yvel,tveld);%expects 1D xvel yvel tvel where u,v are defined on the meshgrid
    %[~,zz]=ode45(@HamEqSolver_BiLin_Irina,tmesh,z0(:),options,u,v,w,xvel,yvel,zvel,tvel);
    xtrd=zz(:,1:2:end-1);ytrd=zz(:,2:2:end);
    clear zz
     
     %convert back to lon/lat
     lattrd=ones(size(ytrd)).*ymin+ytrd./111000;
     lontrd=ones(size(xtrd)).*xmin+xtrd./(111000.*cosd(lattrd));
     figure; plot(lonCoast,latCoast); hold on; plot(lontrd(:,1:5:end),lattrd(:,1:5:end)); title('forward daily field trajectories')
     axis([-6 0 34 38])
     
     disp('entering integration')
     [~,zz]=ode45(@HamEqSolver_BiLin_Irina,tmesh2,z0(:),options,udvel,vdvel,xvel,yvel,tveld);%expects 1D xvel yvel tvel where u,v are defined on the meshgrid
    %[~,zz]=ode45(@HamEqSolver_BiLin_Irina,tmesh,z0(:),options,u,v,w,xvel,yvel,zvel,tvel);
    disp('done')
    xtr2d=zz(:,1:2:end-1); ytr2d=zz(:,2:2:end);
    clear zz
    lattr2d=ones(size(ytr2d)).*ymin+ytr2d./111000;
     lontr2d=ones(size(xtr2d)).*xmin+xtr2d./(111000.*cosd(lattr2d));
     figure; plot(lonCoast,latCoast); hold on; plot(lontr2d(:,1:5:end),lattr2d(:,1:5:end)); title('backward daily field trajectories')
    axis([-6 0 34 38])
    
%%hourly    

tmeshH=0:900:8*86400-1;
    tmesh2H=8*86400-1:-1800:0;
    
    options=odeset('RelTol',10^(-7),'AbsTol',10^(-10));
    %[~,zz]=ode45(@HamEqSolver_Duffing_periodic,tmesh,z0(:),options,eps,a,nu1);
    [~,zz]=ode45(@HamEqSolver_BiLin_Irina,tmeshH,z0(:),options,uhvel,vhvel,xvel,yvel,tvelh);%expects 1D xvel yvel tvel where u,v are defined on the meshgrid
    %[~,zz]=ode45(@HamEqSolver_BiLin_Irina,tmesh,z0(:),options,u,v,w,xvel,yvel,zvel,tvel);
    xtrh=zz(:,1:2:end-1);ytrh=zz(:,2:2:end);
    clear zz
     
     %convert back to lon/lat
     lattrh=ones(size(ytrh)).*ymin+ytrh./111000;
     lontrh=ones(size(xtrh)).*xmin+xtrh./(111000.*cosd(lattrh));
     figure; plot(lonCoast,latCoast); hold on; plot(lontrh(:,1:5:end),lattrh(:,1:5:end)); title('forward hourly field trajectories')
     axis([-6 0 34 38])
     
     disp('entering integration')
     [~,zz]=ode45(@HamEqSolver_BiLin_Irina,tmesh2H,z0(:),options,uhvel,vhvel,xvel,yvel,tvelh);%expects 1D xvel yvel tvel where u,v are defined on the meshgrid
    %[~,zz]=ode45(@HamEqSolver_BiLin_Irina,tmesh,z0(:),options,u,v,w,xvel,yvel,zvel,tvel);
    disp('done')
    xtr2h=zz(:,1:2:end-1); ytr2h=zz(:,2:2:end);
    clear zz
    lattr2h=ones(size(ytr2h)).*ymin+ytr2h./111000;
     lontr2h=ones(size(xtr2h)).*xmin+xtr2h./(111000.*cosd(lattr2h));
     figure; plot(lonCoast,latCoast); hold on; plot(lontr2h(:,1:5:end),lattr2h(:,1:5:end)); title('backward hourly field trajectories')
    axis([-6 0 34 38])
    
    disp('saving')
    %save('dailyhourlyIntegrated8days2.mat')
    %%
    %calculate arclength
    %shape of xtr,ytr is time x trajectories, so xtr(2,:)-xtr(1,:) gives dx for
    %t=1 to 2
      
    dxdth=diff(xtrh);
    dydth=diff(ytrh);
    dsh=sqrt(dxdth.^2 +dydth.^2);
    lengthsh=nansum(dsh);
    
    dxdt2h=diff(xtr2h);
    dydt2h=diff(ytr2h);
    ds2h=sqrt(dxdt2h.^2 +dydt2h.^2);
    lengths2h=nansum(ds2h);
    
    
    dxdtd=diff(xtrd);
    dydtd=diff(ytrd);
    dsd=sqrt(dxdtd.^2 +dydtd.^2);
    lengthsd=nansum(dsd);
    
    dxdt2d=diff(xtr2d);
    dydt2d=diff(ytr2d);
    ds2d=sqrt(dxdt2d.^2 +dydt2d.^2);
    lengths2d=nansum(ds2d);
    
    figure; plot(lonCoast,latCoast,'k'); hold on; pcolorJ(xin2,yin2,reshape(lengthsh,NY,NX)); shading 'flat'; colorbar; title('Hourly Field Arclength')
    colormap(cbrewer('seq','Reds',100)); axis([-6 0 34 38])
    figure; plot(lonCoast,latCoast,'k'); hold on; pcolorJ(xin2,yin2,reshape(lengths2h,NY,NX)); shading 'flat'; colorbar; title('Hourly Field Arclength Negative Time')
    colormap(cbrewer('seq','Reds',100)); axis([-6 0 34 38])
    
    
    figure; plot(lonCoast,latCoast,'k'); hold on; pcolorJ(xin2,yin2,reshape(lengthsd,NY,NX)); shading 'flat'; colorbar; title('Daily Field Arclength')
    colormap(cbrewer('seq','Reds',100)); axis([-6 0 34 38])
    figure; plot(lonCoast,latCoast,'k'); hold on; pcolorJ(xin2,yin2,reshape(lengths2d,NY,NX)); shading 'flat'; colorbar; title('Daily Field Arclength Negative Time')
    colormap(cbrewer('seq','Reds',100)); axis([-6 0 34 38])
    
    figure; plot(lonCoast,latCoast,'k'); hold on; plot(lontrh(:,lengthsh>2e5),lattrh(:,lengthsh>2e5));
%     figure; contour(xin2,yin2,reshape(lengths,NY,NX),[1e4 2e4])
%     colormap(darkb2r(-3e5,3e5)); %axis([-6 0 34 38])
%     hold on; contour(xin2,yin2,-reshape(lengths2,NY,NX),[-2e4 -1e4]); title('Arclength both, 1km')

%%
disp('ftle')
xend=reshape(xtrh(end,:),NY,NX); yend=reshape(ytrh(end,:),NY,NX); 

    dx11=xend(2:end-1,3:end)-xend(2:end-1,1:end-2); 
    dx12=xend(3:end,2:end-1)-xend(1:end-2,2:end-1);
    dy21=yend(2:end-1,3:end)-yend(2:end-1,1:end-2);
    dy22=yend(3:end,2:end-1)-yend(1:end-2,2:end-1);
    dx0 =xinM(2:end-1,3:end)-xinM(2:end-1,1:end-2);
    dy0 =yinM(3:end,2:end-1)-yinM(1:end-2,2:end-1);

    a11=(dx11./dx0).^2.0+(dy21./dx0).*(dy21./dx0);	
    a12=(dx11./dx0).*(dx12./dy0)+(dy21./dx0).*(dy22./dy0);
    a21=(dx12./dy0).*(dx11./dx0)+(dy22./dy0).*(dy21./dx0);
    a22=(dx12./dy0).*(dx12./dy0)+(dy22./dy0).^2.0;
    lambda1=(a22+a11)/2.0+sqrt(((a22+a11)/2.0).^2.0-(a11.*a22-a12.*a21));
    lambda2=(a22+a11)/2.0-sqrt(((a22+a11)/2.0).^2.0-(a11.*a22-a12.*a21));
    ftlePh=log(sqrt(max(lambda1,lambda2)))./abs((tmesh(end)-tmesh(1)));
    
    xend2=reshape(xtr2h(end,:),NY,NX); yend2=reshape(ytr2h(end,:),NY,NX); 

    dx11=xend2(2:end-1,3:end)-xend2(2:end-1,1:end-2); 
    dx12=xend2(3:end,2:end-1)-xend2(1:end-2,2:end-1);
    dy21=yend2(2:end-1,3:end)-yend2(2:end-1,1:end-2);
    dy22=yend2(3:end,2:end-1)-yend2(1:end-2,2:end-1);
    dx0 =xinM(2:end-1,3:end)-xinM(2:end-1,1:end-2);
    dy0 =yinM(3:end,2:end-1)-yinM(1:end-2,2:end-1);

    a11=(dx11./dx0).^2.0+(dy21./dx0).*(dy21./dx0);	
    a12=(dx11./dx0).*(dx12./dy0)+(dy21./dx0).*(dy22./dy0);
    a21=(dx12./dy0).*(dx11./dx0)+(dy22./dy0).*(dy21./dx0);
    a22=(dx12./dy0).*(dx12./dy0)+(dy22./dy0).^2.0;
    lambda1=(a22+a11)/2.0+sqrt(((a22+a11)/2.0).^2.0-(a11.*a22-a12.*a21));
    lambda2=(a22+a11)/2.0-sqrt(((a22+a11)/2.0).^2.0-(a11.*a22-a12.*a21));
    ftleNh=log(sqrt(max(lambda1,lambda2)))./abs((tmesh(end)-tmesh(1)));


xend=reshape(xtrd(end,:),NY,NX); yend=reshape(ytrd(end,:),NY,NX); 

    dx11=xend(2:end-1,3:end)-xend(2:end-1,1:end-2); 
    dx12=xend(3:end,2:end-1)-xend(1:end-2,2:end-1);
    dy21=yend(2:end-1,3:end)-yend(2:end-1,1:end-2);
    dy22=yend(3:end,2:end-1)-yend(1:end-2,2:end-1);
    dx0 =xinM(2:end-1,3:end)-xinM(2:end-1,1:end-2);
    dy0 =yinM(3:end,2:end-1)-yinM(1:end-2,2:end-1);

    a11=(dx11./dx0).^2.0+(dy21./dx0).*(dy21./dx0);	
    a12=(dx11./dx0).*(dx12./dy0)+(dy21./dx0).*(dy22./dy0);
    a21=(dx12./dy0).*(dx11./dx0)+(dy22./dy0).*(dy21./dx0);
    a22=(dx12./dy0).*(dx12./dy0)+(dy22./dy0).^2.0;
    lambda1=(a22+a11)/2.0+sqrt(((a22+a11)/2.0).^2.0-(a11.*a22-a12.*a21));
    lambda2=(a22+a11)/2.0-sqrt(((a22+a11)/2.0).^2.0-(a11.*a22-a12.*a21));
    ftlePd=log(sqrt(max(lambda1,lambda2)))./abs((tmesh(end)-tmesh(1)));
    
    xend2d=reshape(xtr2d(end,:),NY,NX); yend2d=reshape(ytr2d(end,:),NY,NX); 

    dx11=xend2(2:end-1,3:end)-xend2(2:end-1,1:end-2); 
    dx12=xend2(3:end,2:end-1)-xend2(1:end-2,2:end-1);
    dy21=yend2(2:end-1,3:end)-yend2(2:end-1,1:end-2);
    dy22=yend2(3:end,2:end-1)-yend2(1:end-2,2:end-1);
    dx0 =xinM(2:end-1,3:end)-xinM(2:end-1,1:end-2);
    dy0 =yinM(3:end,2:end-1)-yinM(1:end-2,2:end-1);

    a11=(dx11./dx0).^2.0+(dy21./dx0).*(dy21./dx0);	
    a12=(dx11./dx0).*(dx12./dy0)+(dy21./dx0).*(dy22./dy0);
    a21=(dx12./dy0).*(dx11./dx0)+(dy22./dy0).*(dy21./dx0);
    a22=(dx12./dy0).*(dx12./dy0)+(dy22./dy0).^2.0;
    lambda1=(a22+a11)/2.0+sqrt(((a22+a11)/2.0).^2.0-(a11.*a22-a12.*a21));
    lambda2=(a22+a11)/2.0-sqrt(((a22+a11)/2.0).^2.0-(a11.*a22-a12.*a21));
    ftleNd=log(sqrt(max(lambda1,lambda2)))./abs((tmesh(end)-tmesh(1)));
    
    
figure; plot(lonCoast,latCoast,'k')
hold on; pcolorJ(xin2,yin2,ftlePd);
shading 'flat'; title('ftleP daily')
figure; plot(lonCoast,latCoast,'k')
hold on; pcolorJ(xin2,yin2,ftleNd);
shading 'flat'; title('ftleN daily')

figure; plot(lonCoast,latCoast,'k')
hold on; pcolorJ(xin2,yin2,ftlePh);
shading 'flat'; title('ftleP hourly')
figure; plot(lonCoast,latCoast,'k')
hold on; pcolorJ(xin2,yin2,ftleNh);
shading 'flat'; title('ftleN hourly')

disp('saving')
save('dailyhourlyIntegrated8daysFtleEnd.mat')
%%
circleCenters=zeros(20,2);
    %surface:
    circleCenters(1:10,:)=[-5.300000000000,36.080000000000;-4.9000000000000,35.85000000000;-2.9000000000000,35.60000000000;-3.2500,36.1500;-4.38000000000000,35.9300000000000;-3.700000000000,35.450000000000;-2.12000000000000,36.7200000000000;-1.5500000000000,35.30000000000;-1.200000000000,36.60000000000;-2.1200, 35.8400];
    %100m depth:
    %circleCenters(1:12,:)=[-5.33000000000000,36.100000000000;-5.27000000000000,35.9200000000000;-2.94000000000000,35.5200000000000;-3.600,36.300;-4.38000000000000,35.9300000000000;-4.3000000000000,35.400000000000;-4.6000000000000,36.4500000000000;-3.6750,35.400;-2.12000000000000,36.7200000000000;-1.25000000000000,35.700000000000;-1.5000000000000,36.400000000000;-2.1200, 35.8400];
    %200m depth
    %circleCenters(1:12,:)=[-5.33000000000000,36.090000000000;-5.25000000000000,35.900000000000;-2.94000000000000,35.5200000000000;-4.00,36.300;-4.1000000000000,35.800000000000;-4.8000000000000,35.500000000000;-4.6000000000000,36.300000000000;-3.550,35.4500;-2.1000000000000,36.4900000000;-1.000000000000,35.800000000000;-2.1500, 36.100; -1.44,36.92];

    r=4/110;
    ang=0:0.2:2*pi;
    yin=[]; xin=[];
    for i=1:10
    xc=circleCenters(i,1)+r.*cos(ang);
    yc=circleCenters(i,2)+r.*sin(ang);
    yin=cat(1,yin,yc);
    xin=cat(1,xin,xc);
    end
    xin=xin(:);
    yin=yin(:);
    %xmin=min(xin);
    %ymin=min(yin);
    xmin=-5.5;
    xmax=1;
    ymin=35;
    ymax=37;
    xinM=(xin-xmin*ones(size(xin))).*111000.*cosd(yin); 
    yinM=(yin-ymin*ones(size(yin))).*111000;
    z0=[xinM(:).';yinM(:).'];
    %xmin=min(min(XC));
    %ymin=min(min(YC));
    
    disp('entering integration')

tmesh=0:1800:max(tveld)-1;
    tmesh2=max(tveld)-1:-10800:0;
    options=odeset('RelTol',10^(-6),'AbsTol',10^(-9));
    %[~,zz]=ode45(@HamEqSolver_Duffing_periodic,tmesh,z0(:),options,eps,a,nu1);
    [~,zz]=ode45(@HamEqSolver_BiLin_Irina,tmesh,z0(:),options,udvel,vdvel,xvel,yvel,tveld);%expects 1D xvel yvel tvel where u,v are defined on the meshgrid
    %[~,zz]=ode45(@HamEqSolver_BiLin_Irina,tmesh,z0(:),options,u,v,w,xvel,yvel,zvel,tvel);
    xtrdc=zz(:,1:2:end-1);ytrdc=zz(:,2:2:end);
    clear zz
     
     %convert back to lon/lat
     lattrdc=ones(size(ytrdc)).*ymin+ytrdc./111000;
     lontrdc=ones(size(xtrdc)).*xmin+xtrdc./(111000.*cosd(lattrdc));
     figure; plot(lonCoast,latCoast); hold on; plot(lontrdc(:,1:5:end),lattrdc(:,1:5:end)); title('forward daily field trajectories')
     axis([-6 0 34 38])
     
     disp('entering integration')
     [~,zz]=ode45(@HamEqSolver_BiLin_Irina,tmesh2,z0(:),options,udvel,vdvel,xvel,yvel,tveld);%expects 1D xvel yvel tvel where u,v are defined on the meshgrid
    %[~,zz]=ode45(@HamEqSolver_BiLin_Irina,tmesh,z0(:),options,u,v,w,xvel,yvel,zvel,tvel);
    disp('done')
    xtr2dc=zz(:,1:2:end-1); ytr2dc=zz(:,2:2:end);
    clear zz
    lattr2dc=ones(size(ytr2dc)).*ymin+ytr2dc./111000;
     lontr2dc=ones(size(xtr2dc)).*xmin+xtr2dc./(111000.*cosd(lattr2dc));
     figure; plot(lonCoast,latCoast); hold on; plot(lontr2dc(:,1:5:end),lattr2dc(:,1:5:end)); title('backward daily field trajectories')
    axis([-6 0 34 38])
    
%%hourly    

tmeshH=0:900:max(tvelh)-1;
    tmesh2H=max(tvelh)-1:-900:0;
    
    options=odeset('RelTol',10^(-6),'AbsTol',10^(-9));
    %[~,zz]=ode45(@HamEqSolver_Duffing_periodic,tmesh,z0(:),options,eps,a,nu1);
    [~,zz]=ode45(@HamEqSolver_BiLin_Irina,tmesh,z0(:),options,uhvel,vhvel,xvel,yvel,tvelh);%expects 1D xvel yvel tvel where u,v are defined on the meshgrid
    %[~,zz]=ode45(@HamEqSolver_BiLin_Irina,tmesh,z0(:),options,u,v,w,xvel,yvel,zvel,tvel);
    xtrhc=zz(:,1:2:end-1);ytrhc=zz(:,2:2:end);
    clear zz
     
     %convert back to lon/lat
     lattrhc=ones(size(ytrhc)).*ymin+ytrhc./111000;
     lontrhc=ones(size(xtrhc)).*xmin+xtrhc./(111000.*cosd(lattrhc));
     figure; plot(lonCoast,latCoast); hold on; plot(lontrhc(:,1:5:end),lattrhc(:,1:5:end)); title('forward hourly field trajectories')
     axis([-6 0 34 38])
     
     disp('entering integration')
     [~,zz]=ode45(@HamEqSolver_BiLin_Irina,tmesh2,z0(:),options,uhvel,vhvel,xvel,yvel,tvelh);%expects 1D xvel yvel tvel where u,v are defined on the meshgrid
    %[~,zz]=ode45(@HamEqSolver_BiLin_Irina,tmesh,z0(:),options,u,v,w,xvel,yvel,zvel,tvel);
    disp('done')
    xtr2hc=zz(:,1:2:end-1); ytr2hc=zz(:,2:2:end);
    clear zz
    lattr2hc=ones(size(ytr2hc)).*ymin+ytr2hc./111000;
     lontr2hc=ones(size(xtr2hc)).*xmin+xtr2hc./(111000.*cosd(lattr2hc));
     figure; plot(lonCoast,latCoast); hold on; plot(lontr2hc(:,1:5:end),lattr2hc(:,1:5:end)); title('backward hourly field trajectories')
    axis([-6 0 34 38])
    
    
    figure; plot(lonCoast,latCoast); hold on; plot(lontr2hc,lattr2hc,'b'); plot(lontrhc,lattrhc,'r'); title('Hourly')
    axis([-6 0 34 38])
    
        figure; plot(lonCoast,latCoast); hold on; plot(lontr2dc,lattr2dc,'b'); plot(lontrdc,lattrdc,'r'); title('Daily')
    axis([-6 0 34 38])