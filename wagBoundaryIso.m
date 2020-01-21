addpath('../mStuff')
load('uDaily175.mat')
load('rhoVary175.mat')
load('geometrySpinupSteady.mat')
dBin=0.5.*dInterface(1:end-1)+0.5.*dInterface(2:end);
load('vDaily175.mat','V')
load('wDaily175.mat','W')

circleCenters=zeros(2,2);
    %148 steady for WAG
    circleCenters(1:2,:)=[-5.25,35.85;-3.15,35.35];
    
r=15/110;
    ang=0:0.25:2*pi;%was 0.1 steps for manifolds %0.02:2*pi;
    yin=[]; xin=[];
   % for i=1:2
    xc1=circleCenters(1,1)+r.*cos(ang);
    yc1=circleCenters(1,2)+r.*sin(ang);
    %yin=cat(1,yin,yc);
    %xin=cat(1,xin,xc);
   xc2=circleCenters(2,1)+r.*cos(ang);
    yc2=circleCenters(2,2)+r.*sin(ang);
    
    xinM1=(xc1-min(min(XC))*ones(size(xc1))).*111000.*cosd(yc1); 
    yinM1=(yc1-min(min(YC))*ones(size(yc1))).*111000;
    xinM2=(xc2-min(min(XC))*ones(size(xc2))).*111000.*cosd(yc2); 
    yinM2=(yc2-min(min(YC))*ones(size(yc2))).*111000;
    %z0=[xinM(:).';yinM(:).'];
    xmin=min(min(XC));
    ymin=min(min(YC));
    
    %z0f=[xinM1(:)';yinM1(:)'];
    %z0b=[xinM2(:)';yinM2(:)'];
    x0f=xinM1(:)';
    y0f=yinM1(:)';
    x0b=xinM2(:)';
    y0b=yinM2(:)';
    
    tvel=10*times-times(1);
    xmin=min(min(XC));
ymin=min(min(YC));
xinM=111000*cosd(YC).*(XC-xmin*ones(size(XC)));
yinM=111000*(YC-ymin.*ones(size(YC)));
%%
disp('entering day loop')

for j=35;%:10:145; %9:1:30 %0:0.5:21
    j
    clear *tr* rho1* rho2*
    %set times
    t1=86400*(j-8):8640:(j)*86400;
    t2=(j+8)*86400:-8640:(j)*86400;
    options=odeset('RelTol',10^(-3),'AbsTol',10^(-6));%was -10,-13 for manifolds
    %find z0f,z0b
    isopyc=1025.5;
    for depthi=1:36
    rho1(:,depthi)=griddata(xinM(230:500,1:130),yinM(230:500,1:130),Rho(:,:,depthi,j-8),x0f,y0f);
    rho2(:,depthi)=griddata(xinM(230:500,1:130),yinM(230:500,1:130),Rho(:,:,depthi,j+8),x0b,y0b);
    end
    depths1=3:dBin(36)-1;
    for ii=1:length(x0f)
        rho11=interp1(dBin(1:36),rho1(ii,:),depths1);
        rho22=interp1(dBin(1:36),rho2(ii,:),depths1);
        %figure; plot(rho11,-depths1,'g'); hold on; plot(rho22,-depths1,'r')
        %di11=find(rho11<isopyc,1,'last');
        di12=find(abs(rho11-isopyc)==min(abs(rho11-isopyc)),1,'first');       
        %di21=find(rho22<isopyc,1,'last');
        di22=find(abs(rho22-isopyc)==min(abs(rho22-isopyc)),1,'first');
        if isempty(di12)
            z0f(ii)=0;
        else
            z0f(ii)=depths1(di12);
        end
        if isempty(di22)
            z0b(ii)=0;
        else
            z0b(ii)=depths1(di22);
        end
    end
    
    %integrations
    for k=1:length(t1)-1 % different integration time
        k
    r0f=[x0f;y0f;z0f];
    r0b=[x0b;y0b;z0b];
    tmesh1=[t1(k) (t1(k+1)+t1(k))/2 t1(k+1)];
    tmesh2=[t2(k) (t2(k+1)+t2(k))/2 t2(k+1)];
    
    disp('int forward')
    [~,zz]=ode45(@HamEqSolver_TriLin,tmesh1,r0f(:),options,U,V,W,xvel,yvel,dBin,tvel);%expects 1D xvel yvel tvel where u,v are defined on the meshgrid
        if k<length(t1)-1 %if not the last integration
           
            dlmin=2000; %RR/3; %10^(-3);

            xunst=[zz(end,1:3:end-2) zz(end,1)]; %end of previous integration, make a loop
            yunst=[zz(end,2:3:end-1) zz(end,2)];
            zunst=[zz(end,3:3:end)   zz(end,3)];

            tunst=(1:1:length(xunst)); t1unst=(1:0.005:length(xunst));
            x1unst=interp1(tunst,xunst,t1unst); y1unst=interp1(tunst,yunst,t1unst); z1unst=interp1(tunst,zunst,t1unst);%interpolate to 200x as many points

            dxunst=diff(x1unst); dyunst=diff(y1unst); dzunst=diff(z1unst); dlunst=sqrt(dxunst.^2+dyunst.^2+dzunst.^2); %distances between points

            nunst=[];

            ii=1; 
            while ii<(length(dlunst)-1) %while index ii is less than the number of points
                if dlunst(ii)>dlmin
                    nunst=[nunst ii]; ii=ii+1; %if this distance is bigger than the minimum, add ii to list nunst
                else 
                    jj=1;
                    while ((dlunst(ii)+dlunst(ii+jj)<dlmin)&&((ii+jj)<length(dlunst)))
                          dlunst(ii)=dlunst(ii)+dlunst(ii+jj); jj=jj+1;  %otherwise, find how many points forward you can move until it gets too far
                     end
                    nunst=[nunst ii+jj]; ii=ii+jj; %add the next index, where the distance just goes over dlmin, to nunst
                end
                %ii
            end
            x0f=x1unst(nunst); y0f=y1unst(nunst); z0f=z1unst(nunst);%set input to next integration to all the interpolated points that are needed to have separations no larger than dlmin
            %%%%%%%% end of resampling 1 %%%%%%%

        else
                 x0f=zz(end,1:3:end-2);
                 y0f=zz(end,2:3:end-1);
                 z0f=zz(end,3:3:end);
        end %if not last integration
 clear zz
     disp('entering integration 2')
     [~,zz]=ode45(@HamEqSolver_TriLin,tmesh2,r0b(:),options,U,V,W,xvel,yvel,dBin,tvel);%expects 1D xvel yvel tvel where u,v are defined on the meshgrid
         if k<length(t1)-1 %if not the last integration
           
            dlmin=2000; %RR/3; %10^(-3);

            xunst=[zz(end,1:3:end-2) zz(end,1)]; %end of previous integration, make a loop
            yunst=[zz(end,2:3:end-1) zz(end,2)];
            zunst=[zz(end,3:3:end)   zz(end,3)];

            tunst=(1:1:length(xunst)); t1unst=(1:0.005:length(xunst));
            x1unst=interp1(tunst,xunst,t1unst); y1unst=interp1(tunst,yunst,t1unst); z1unst=interp1(tunst,zunst,t1unst);%interpolate to 200x as many points

            dxunst=diff(x1unst); dyunst=diff(y1unst); dzunst=diff(z1unst); dlunst=sqrt(dxunst.^2+dyunst.^2+dzunst.^2); %distances between points

            nunst=[];

            ii=1; 
            while ii<(length(dlunst)-1) %while index ii is less than the number of points
                if dlunst(ii)>dlmin
                    nunst=[nunst ii]; ii=ii+1; %if this distance is bigger than the minimum, add ii to list nunst
                else 
                    jj=1;
                    while ((dlunst(ii)+dlunst(ii+jj)<dlmin)&&((ii+jj)<length(dlunst)))
                          dlunst(ii)=dlunst(ii)+dlunst(ii+jj); jj=jj+1;  %otherwise, find how many points forward you can move until it gets too far
                     end
                    nunst=[nunst ii+jj]; ii=ii+jj; %add the next index, where the distance just goes over dlmin, to nunst
                end
                %ii
            end
            x0b=x1unst(nunst); y0b=y1unst(nunst); z0b=z1unst(nunst);%set input to next integration to all the interpolated points that are needed to have separations no larger than dlmin
            %%%%%%%% end of resampling 1 %%%%%%%

        else
                 x0b=zz(end,1:3:end-2);
                 y0b=zz(end,2:3:end-1);
                 z0b=zz(end,3:3:end);
        end %if not last integration

    
    end %different integration times
    
    %convert back to lon/lat
     lattrF=ones(size(y0f)).*ymin+y0f./111000;
     lontrF=ones(size(x0f)).*xmin+x0f./(111000.*cosd(lattrF));
     %figure; plot(lonCoast,latCoast); hold on; plot(lontr(:,1:5:end),lattr(:,1:5:end)); title('forward trajectories')
     %axis([-6 0 34 38])
     
    lattrB=ones(size(y0b)).*ymin+y0b./111000;
     lontrB=ones(size(x0b)).*xmin+x0b./(111000.*cosd(lattrB));

    zF=z0f; zB=z0b;
    
   fn=strcat('wagBoundaryday',num2str(j),'Iso',num2str(isopyc),'int8.mat');
   save(fn,'-v7.3')
    figure
    plot(lonCoast,latCoast,'k')
    hold on
    plot(lontrF(end,:),lattrF(end,:),'ro')
    plot(lontrB(end,:),lattrB(end,:),'bo')
    title(strcat('Manifolds 8 day integration day',num2str(j)))
    axis([-7 -1 34 38])
    fn=strcat('wagBoundaryday',num2str(j),'Iso',num2str(isopyc),'int8.pdf');
    save2pdf(fn)
    
   disp('done j')
end %day loop
