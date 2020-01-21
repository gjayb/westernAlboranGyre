
addpath('../mStuff')
%load('uvAdepth1interpolated.mat')
load('uvwIso265InterpNF.mat', 'XC','YC','U','V','xvel','yvel','tvel','*Coast')
%load('uvAdepth1interpolated.mat', 'YC')
%load('uvAdepth1interpolated.mat', 'u');
%load('uvAdepth1interpolated.mat', 'v');
%load('uvAdepth1interpolated.mat', 'yvel')
%load('uvAdepth1interpolated.mat', 'xvel')
%load('uvAdepth1interpolated.mat', 'tvel')
%load('uvAdepth1interpolated.mat', 'latCoast')
%load('uvAdepth1interpolated.mat', 'lonCoast')

%tvel
u=U; clear U; v=V; clear V; %for isopycnalsNF
[nx,ny,nt]=size(u)
tvel=0:86400:(nt-1)*86400;

circleCenters=zeros(2,2);
    %148 steady for WAG
    circleCenters(1:2,:)=[-5.35,35.86;-3.05,35.4];


    
r=8/110;

    ang=0:0.2:2*pi;%0.02:2*pi;
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
    
    %z0f=[xinM1(:)';yinM1(:)'];
    %z0b=[xinM2(:)';yinM2(:)'];
    x0f=xinM1(:)';
    y0f=yinM1(:)';
    x0b=xinM2(:)';
    y0b=yinM2(:)';
%%
disp('entering day loop')

%notaBene='day 67 was skipped, could not run in under 12 hours';

%notaBene1='dlmin changes from 1km to 2km at day 68';

for j=118:-1:90;%[9:18 26:57 65:102 107:112 118:127 134:139];%[19:25 58:64 103:106 113:117  128:133];%113:117; %[10 31:1:61] %0:0.5:21
    j
    clear *tr*
    t1=86400*(j-14):86400:(j)*86400;
    t2=(j+14)*86400:-86400:(j)*86400;
    options=odeset('RelTol',10^(-10),'AbsTol',10^(-13));

x0f=xinM1(:)';
y0f=yinM1(:)';
x0b=xinM2(:)';
y0b=yinM2(:)';
    
    for k=1:length(t1)-1 % different integration time
	j
        k
    z0f=[x0f;y0f];
    z0b=[x0b;y0b];
    tmesh1=[t1(k) (t1(k+1)+t1(k))/2 t1(k+1)];
    tmesh2=[t2(k) (t2(k+1)+t2(k))/2 t2(k+1)];

    
    disp('int forward')
    [~,zz]=ode45(@HamEqSolver_BiLin_Irina,tmesh1,z0f(:),options,u,v,xvel,yvel,tvel);%expects 1D xvel yvel tvel where u,v are defined on the meshgrid
disp('int done')   
     if k<length(t1)-1 %if not the last integration
           
            dlmin=2000  %RR/3; %10^(-3);

            xunst=[zz(end,1:2:end-1) zz(end,1)]; %end of previous integration, make a loop
            yunst=[zz(end,2:2:end) zz(end,2)];

            tunst=(1:1:length(xunst)); t1unst=(1:0.005:length(xunst));
            x1unst=interp1(tunst,xunst,t1unst); y1unst=interp1(tunst,yunst,t1unst); %interpolate to 200x as many points

            dxunst=diff(x1unst); dyunst=diff(y1unst); dlunst=sqrt(dxunst.^2+dyunst.^2); %distances between points

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
            x0f=x1unst(nunst); y0f=y1unst(nunst); %set input to next integration to all the interpolated points that are needed to have separations no larger than dlmin
            %%%%%%%% end of resampling 1 %%%%%%%

        else
                 x0f=zz(end,1:2:end-1);
                 y0f=zz(end,2:2:end);
            
        end
 clear zz
     disp('entering integration 2')
     [~,zz]=ode45(@HamEqSolver_BiLin_Irina,tmesh2,z0b(:),options,u,v,xvel,yvel,tvel);%expects 1D xvel yvel tvel where u,v are defined on the meshgrid
disp('int done')      
   if k<length(t1)-1 %if not the last integration
           
            dlmin=2000; %RR/3; %10^(-3);

            xunst=[zz(end,1:2:end-1) zz(end,1)]; %end of previous integration, make a loop
            yunst=[zz(end,2:2:end) zz(end,2)];

            tunst=(1:1:length(xunst)); t1unst=(1:0.005:length(xunst));
            x1unst=interp1(tunst,xunst,t1unst); y1unst=interp1(tunst,yunst,t1unst); %interpolate to 200x as many points

            dxunst=diff(x1unst); dyunst=diff(y1unst); dlunst=sqrt(dxunst.^2+dyunst.^2); %distances between points

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
            x0b=x1unst(nunst); y0b=y1unst(nunst); %set input to next integration to all the interpolated points that are needed to have separations no larger than dlmin
            %%%%%%%% end of resampling 1 %%%%%%%

        else
                 x0b=zz(end,1:2:end-1);
                 y0b=zz(end,2:2:end);
            
        end
  
    
     
    end
    
    %%

    %%

    
%%
     
     %convert back to lon/lat
disp('convert to lon/lat')
     lattrF=ones(size(y0f)).*ymin+y0f./111000;
     lontrF=ones(size(x0f)).*xmin+x0f./(111000.*cosd(lattrF));
     %figure; plot(lonCoast,latCoast); hold on; plot(lontr(:,1:5:end),lattr(:,1:5:end)); title('forward trajectories')
     %axis([-6 0 34 38])
     
    lattrB=ones(size(y0b)).*ymin+y0b./111000;
     lontrB=ones(size(x0b)).*xmin+x0b./(111000.*cosd(lattrB));
    
   fn=strcat('wagBoundaryday',num2str(j),'Iso265Int14.mat');
%fn='wagBoundarySteadySurfaceInt28.mat';   
save(fn,'-v7.3')
   % figure
   % plot(lonCoast,latCoast,'k')
   % hold on
   % plot(lontrF(end,:),lattrF(end,:),'r')
   % plot(lontrB(end,:),lattrB(end,:),'b')
   % title(strcat('Manifolds 8 day integration hourly vel day',num2str(j)))
   % fn=strcat('wagBoundaryday',num2str(j),'Iso275int8.pdf');
   % save2pdf(fn)
    
   disp('done j')
end
disp('done')

