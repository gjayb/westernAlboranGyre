%vorticity budget, euler WAG, surface
%uses circulationInt.m and crossCurveFlux.m
%% load for circulation
load('edgesWAGeuler2017NF.mat', 'inWag','XC','YC')
load('geometrySpinupSteady.mat','*Coast')
xmin=min(XC(:)); ymin=min(YC(:));
xm=111111*cosd(YC).*(XC-xmin*ones(size(XC)));%was using 111000
ym=111111*(YC-ymin*ones(size(YC)));
xcoast=111111*cosd(latCoast).*(lonCoast-xmin*ones(size(lonCoast)));
ycoast=111111*(latCoast-ymin*ones(size(latCoast)));
%load('uvwSSHDailyDepth1rotated148F.mat','Urot','Vrot','SSHa')

%% edge
inWag2=inWag(:,:,1).*double(d>300);
inWag3=inWag(:,:,1).*double(d>200);
figure; [c,h]=contour(xm,ym,inWag(:,:,1),[1 1]);
edgeWagXo=c(1,2:511);
edgeWagYo=c(2,2:511);
figure; [c,h]=contour(xm,ym,inWag2(:,:,1),[1 1]);
edgeWagX=c(1,2:459);
edgeWagY=c(2,2:459);
figure; [c,h]=contour(xm,ym,inWag3(:,:,1),[1 1]);
figure; [c,h]=contour(xm,ym,inWag2(:,:,1),[1 1]);
edgeWagX=c(1,2:444);
edgeWagY=c(2,2:444);
hold on; plot(edgeWagX,edgeWagY,'b');plot(edgeWagXo,edgeWagYo,'r')
%%
%dEdge=griddata(xm,ym,d,edgeWagX,edgeWagY); figure; plot(dEdge)
%%
%edgeWagX=edgeWagX(98:350);
%edgeWagY=edgeWagY(98:350);
%% interpolate edge
xunst=[xSquare.' xSquare(1)];%[edgeWagX edgeWagX(1)];
yunst=[ySquare.' ySquare(1)];%[edgeWagY edgeWagY(1)];
tunst=(1:1:length(xunst)); t1unst=(1:0.0001:length(xunst));
x1unst=interp1(tunst,xunst,t1unst); y1unst=interp1(tunst,yunst,t1unst); %interpolate to 200x as many points
dxunst=diff(x1unst); dyunst=diff(y1unst); dlunst=sqrt(dxunst.^2+dyunst.^2); %distances between points

nunst=[];
ii=1; 
dlminset=[10 51 101 250 501 1001];
for nsi=1:6
    dlmin=dlminset(nsi)
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
%     if mod(ii,10000)==0
%         disp(num2str(ii))
%     end
end
switch nsi
    case 1
        x10=x1unst(nunst); y10=y1unst(nunst); 
    case 2
        x50=x1unst(nunst); y50=y1unst(nunst); 
    case 3
         x100=x1unst(nunst); y100=y1unst(nunst); 
    case 4
        x250=x1unst(nunst); y250=y1unst(nunst); 
    case 5
        x500=x1unst(nunst); y500=y1unst(nunst); 
    case 6
        x1000=x1unst(nunst); y1000=y1unst(nunst); 
%     case 5
%         x2000=x1unst(nunst); y2000=y1unst(nunst); 
end
end
%figure; plot(edgeWagX,edgeWagY); hold all; plot(x100,y100)
% edgeOldX=edgeWagX;
% edgeWagX=x0b;
% edgeOldY=edgeWagY;
% edgeWagY=y0b;
%% edge using edges of cells
% load('edgesWAGeuler2017NF.mat','open*')
% openE=openE(:,:,1);
% openW=openW(:,:,1);
% openN=openN(:,:,1);
% openS=openS(:,:,1);
% edgeX(1,1)=XG(247,103); edgeX(1,2)=XG(247,104);
%% edge through cells
load('geometrySpinupSteady')
xmin=min(XC(:)); ymin=min(YC(:));
xm=111111*cosd(YC).*(XC-xmin*ones(size(XC)));
ym=111111*(YC-ymin*ones(size(YC)));
xcoast=111111*cosd(latCoast).*(lonCoast-xmin*ones(size(lonCoast)));
ycoast=111111*(latCoast-ymin*ones(size(latCoast)));
xum=111111*cosd(YU).*(XU-xmin*ones(size(XU)));
yum=111111*(YU-ymin*ones(size(YU)));
xvm=111111*cosd(YV).*(XV-xmin*ones(size(XV)));
yvm=111111*(YV-ymin*ones(size(YV)));
xi1=375;%300
xi2=425;%350
yi1=100;%60
yi2=150;%100
x1=xm(xi1,yi1:yi2);%left bottom to top, !corners!
y1=ym(xi1,yi1:yi2);
x2=xm(xi1:xi2,yi2);%top left to right
y2=ym(xi1:xi2,yi2);
x3=xm(xi2,yi1:yi2);%right bottom to top
y3=ym(xi2,yi1:yi2);
x4=xm(xi1:xi2,yi1);%bottom left to right
y4=ym(xi1:xi2,yi1);
x2=x2(:); y2=y2(:); x1=x1(:); y1=y1(:); 
x4=x4(:); y4=y4(:); x3=x3(:); y3=y3(:); 
xSquare=[x4(1:end-1); x3(1:end-1);x2(end:-1:2);x1(end:-1:1)];
ySquare=[y4(1:end-1); y3(1:end-1);y2(end:-1:2);y1(end:-1:1)];
ulogic2=false([700 200]); ulogic4=ulogic2; ulogic4(xi1+1:xi2,yi1)=true; ulogic2(xi1+1:xi2,yi2)=true;
vlogic1=false([700 200]); vlogic3=vlogic1; vlogic1(xi1,yi1+1:yi2)=true; vlogic3(xi2,yi1+1:yi2)=true;
load('distancesAreas')
DYU=reshape(dyu,[700 200]);
DXV=reshape(dxv,[700 200]);
DYC=reshape(dyc,[700 200]);
DXC=reshape(dxc,[700 200]);
%earthR=6371.0*1e3;
ds1=DYC(xi1,yi1+1:yi2);%sqrt(diff([x4(1); x1(:)]).^2 +diff([y4(1); y1(:)]).^2);%
ds2=DXC(xi1+1:xi2,yi2);%sqrt(diff([x1(end); x2(:)]).^2 +diff([y1(end);y2(:)]).^2);%
ds3=DYC(xi2,yi1+1:yi2);%sqrt(diff([x3(:); x2(end)]).^2 +diff([y3(:); y2(end)]).^2);%
ds4=DXC(xi1+1:xi2,yi1);%sqrt(diff([x4(:); x3(1)]).^2 +diff([y4(:); y3(1)]).^2);%
ds1=ds1(:);ds2=ds2(:);ds3=ds3(:);ds4=ds4(:);
%[~,~,nt]=size(Uext);
dx=diff([xSquare(:); xSquare(1)]); dy=diff([ySquare(:); ySquare(1)]);
vortLogic=false([700 200]); vortLogic(xi1+1:xi2,yi1+1:yi2)=true;
%%
load('uvwDailyNativeNF.mat','U');
U=squeeze(U(:,:,1,:));
load('uvwDailyNativeNF.mat','V');
V=squeeze(V(:,:,1,:));
load('vorticitySurface')
load('distancesAreas')
%%
rAz=reshape(raz,[700 200]);
for i=1:148%nt
    hold1=Vtend(:,:,i);
    hold2=Utend(:,:,i);
    tendS(i)=sum(-hold1(vlogic1).*ds1)+sum(hold1(vlogic3).*ds3)+sum(-hold2(ulogic2).*ds2)+sum(hold2(ulogic4).*ds4);
    %tendVals(:,i)=[-hold1(vlogic1);hold2(ulogic2);hold1(vlogic3);-hold2(ulogic4)];
    hold3=vort(:,:,i);
    vortint(i)=squeeze(nansum(nansum(hold3.*rAz.*double(vortLogic))));
    hold4=U(:,:,i);
    hold5=V(:,:,i);
    circ(i)=sum(-hold5(vlogic1).*ds1)+sum(hold5(vlogic3).*ds3)+sum(-hold4(ulogic2).*ds2)+sum(hold4(ulogic4).*ds4);
    hold3=vorta(:,:,i);
    vortAInt(i)=squeeze(nansum(nansum(hold3.*rAz.*double(vortLogic))));
end
dvortdT=tendS; dvortdT(2:end-1)=(vortint(3:end)-vortint(1:end-2));
dvortadT=tendS; dvortadT(2:end-1)=(vortAInt(3:end)-vortAInt(1:end-2));
figure; plot(tendS); hold all; plot(dvortdT); plot(dvortadT);
figure; plot(circ); hold all; plot(vortAInt);
%figure; scatter(xSquare,ySquare,36,tendVals(:,1),'filled');
%%
nt=148;
tendE1=Utend.*repmat(AngleCS,[1 1 nt]) - Vtend.*repmat(AngleSN,[1 1 nt]);  
tendN1=Utend.*repmat(AngleSN,[1 1 nt]) + Vtend.*repmat(AngleCS,[1 1 nt]); 
[ tendS1, sign0 ] = circulationInt( xum,yum,tendE1,xvm,yvm,tendN1,xSquare.',ySquare.');
figure; plot(tendS1); hold all; plot(diff(vortInt));
[a3,a4]=corrcoef(tendS1(1:end-1),diff(vortInt))
[ circ1, sign01 ] = circulationInt( xm,ym,Urot,xm,ym,Vrot,xSquare.',ySquare.');
figure; plot(circ1); hold all; plot(circ); plot(vortAInt);
figure; plot(diff(circ)); hold all; plot(diff(vortAInt))
%%
[ circ10c, ~ ] = circulationInt( xm,ym,Urot,xm,ym,Vrot,x10,y10);
[ circ50c, ~ ] = circulationInt( xm,ym,Urot,xm,ym,Vrot,x50,y50);
%[ circ100c, ~ ] = circulationInt3( xm,ym,Urot,xm,ym,Vrot,x100,y100);
%[ circ250c, ~ ] = circulationInt3( xm,ym,Urot,xm,ym,Vrot,x250,y250);
%[ circ500c, ~ ] = circulationInt3( xm,ym,Urot,xm,ym,Vrot,x500,y500);
[ circ1000c, ~ ] = circulationInt( xm,ym,Urot,xm,ym,Vrot,x1000,y1000);
figure; plot(vortAInt,'linewidth',2); hold all; plot(circ1); plot(circ10); plot(circ50); plot(circ10c); plot(circ50c);% plot(circ500); plot(circ1000);
legend('vorticity integral','native grid','10m','50m','10m c','50m c')%,'500m','1000m')
%%
[~,~,nt]=size(Uext);
% windE=Uext.*repmat(AngleCS,[1 1 nt]) - Vext.*repmat(AngleSN,[1 1 nt]);  
% windN=Uext.*repmat(AngleSN,[1 1 nt]) + Vext.*repmat(AngleCS,[1 1 nt]); 
% diffE=UDif2a.*repmat(AngleCS,[1 1 nt]) - VDif2a.*repmat(AngleSN,[1 1 nt]);  
% diffN=UDif2a.*repmat(AngleSN,[1 1 nt]) + VDif2a.*repmat(AngleCS,[1 1 nt]); 
% diffE2=UDif2b.*repmat(AngleCS,[1 1 nt]) - VDif2b.*repmat(AngleSN,[1 1 nt]);  
% diffN2=UDif2b.*repmat(AngleSN,[1 1 nt]) + VDif2b.*repmat(AngleCS,[1 1 nt]);
% dissE=UDiss.*repmat(AngleCS,[1 1 nt]) - VDiss.*repmat(AngleSN,[1 1 nt]);  
% dissN=UDiss.*repmat(AngleSN,[1 1 nt]) + VDiss.*repmat(AngleCS,[1 1 nt]); 
% abE=AbU.*repmat(AngleCS,[1 1 nt]) - AbV.*repmat(AngleSN,[1 1 nt]);  
% abN=AbU.*repmat(AngleSN,[1 1 nt]) + AbV.*repmat(AngleCS,[1 1 nt]); 
% tendE=Utend.*repmat(AngleCS,[1 1 nt])./86400 - Vtend.*repmat(AngleSN,[1 1 nt])./86400;  
% tendN=Utend.*repmat(AngleSN,[1 1 nt])./86400 + Vtend.*repmat(AngleCS,[1 1 nt])./86400; 
% press1E=UdPdx.*repmat(AngleCS,[1 1 nt]) - VdPdy.*repmat(AngleSN,[1 1 nt]);  
% press1N=UdPdx.*repmat(AngleSN,[1 1 nt]) + VdPdy.*repmat(AngleCS,[1 1 nt]); 
% advE=AdvU.*repmat(AngleCS,[1 1 nt]) - AdvV.*repmat(AngleSN,[1 1 nt]);  
% advN=AdvU.*repmat(AngleSN,[1 1 nt]) + AdvV.*repmat(AngleCS,[1 1 nt]);
% 
% [ winds, ~ ] = circulationInt( xum,yum,windE,xvm,yvm,windN,xSquare.',ySquare.' );
% [ pclins, ~ ] = circulationInt( xum,yum,press1E,xvm,yvm,press1N,xSquare.',ySquare.');
% [ tends, ~ ] = circulationInt( xum,yum,tendE,xvm,yvm,tendN,xSquare.',ySquare.');
% [ diffs, ~ ] = circulationInt( xum,yum,diffE2,xvm,yvm,diffN2,xSquare.',ySquare.' );
% [ disss, ~ ] = circulationInt( xum,yum,dissE,xvm,yvm,dissN,xSquare.',ySquare.' );
% [ times, ~ ] = circulationInt( xum,yum,abE,xvm,yvm,abN,xSquare.',ySquare.' );
% [ advs, ~ ] = circulationInt( xum,yum,advE,xvm,yvm,advN,xSquare.',ySquare.' );
figure; plot(advs+times+pclins+winds+diffs+disss-tends,'linewidth',3); hold all; plot(advs); plot(times); plot(pclins); plot(winds); plot(diffs); plot(disss); plot(-tends)
%% calculate circulation for each day
[ integral1, sign1 ] = circulationInt( xum,yum,Urot,xvm,yvm,Vrot,edgeWagX,edgeWagY );
%[ integralc, ~ ] = circulationInt3( xum,yum,Urot,xvm,yvm,Vrot,edgeWagX,edgeWagY );
disp('1')
[ integral1o, sign1o ] = circulationInt( xum,yum,Urot,xvm,yvm,Vrot,edgeWagXo,edgeWagYo );
% [ integral1000, ~ ] = circulationInt( xm,ym,Urot,Vrot,x1000,y1000 );
% [ integral2000, ~ ] = circulationInt( xm,ym,Urot,Vrot,x2000,y2000 );
% disp('first 3 done')
%[ integral100, ~] = circulationInt( xm,ym,Urot,Vrot,x100,y100 );
%[ integral250, ~ ] = circulationInt( xm,ym,Urot,Vrot,x250,y250 );
%[ integral100c, ~] = circulationInt3( xm,ym,Urot,Vrot,x100,y100 );
%[ integral250c, ~ ] = circulationInt3( xm,ym,Urot,Vrot,x250,y250 );
%disp('first 5 done')
%[ integral500, ~ ] = circulationInt( xm,ym,Urot,Vrot,x500,y500 );
%[ integral500c, ~ ] = circulationInt3( xm,ym,Urot,Vrot,x500,y500 );

% figure; plot(integral2000); hold all; plot(integral1000); plot(integral500); plot(integral250); plot(integral100); plot(integral1)
% figure; plot(integral2000-integral1000); hold all; plot(integral1000-integral500); plot(integral500-integral250); plot(integral250-integral100)
dCircdt=diff(integral1)./86400;
% vorticity area integral
load('vorticitySurface.mat', 'vorta','vort')
load('distancesAreas')
rAz=reshape(raz,[700 200]);
vortInt=squeeze(nansum(nansum(vort(:,:,1:149).*repmat(inWag(:,:,1).*rAz.*hFacC(:,:,1),[1 1 149]))));
vortaInt=squeeze(nansum(nansum(vorta(:,:,1:149).*repmat(inWag(:,:,1).*rAz.*hFacC(:,:,1),[1 1 149]))));
dVortdtA=diff(vortaInt)./86400;%daily averages
dVortdtS=diff(vortInt)./86400;%snapshots

vortInt3=squeeze(nansum(nansum(vort(:,:,1:149).*repmat(inWag3(:,:,1).*rAz.*hFacC(:,:,1),[1 1 149]))));
vortaInt3=squeeze(nansum(nansum(vorta(:,:,1:149).*repmat(inWag3(:,:,1).*rAz.*hFacC(:,:,1),[1 1 149]))));
dVortdtA2=diff(vortaInt2)./86400;%daily averages
dVortdtS2=diff(vortInt2)./86400;%snapshots
%% calculate wind stress component
load('momentumWindForcing.mat','*ext')
%rotate to east-north
load('geometrySpinupSteady','Angle*')
[~,~,nt]=size(Uext);
windE=Uext.*repmat(AngleCS,[1 1 nt]) - Vext.*repmat(AngleSN,[1 1 nt]);  
windN=Uext.*repmat(AngleSN,[1 1 nt]) + Vext.*repmat(AngleCS,[1 1 nt]); 
[ wind, sign2 ] = circulationInt( xum,yum,windE,xvm,yvm,windN,edgeWagX,edgeWagY );
[ windc, ~ ] = circulationInt3( xum,yum,windE,xvm,yvm,windN,edgeWagX,edgeWagY );

% edgeWagXH=0.5*(edgeWagX(2:end)+edgeWagX(1:end-1));
% edgeWagYH=0.5*(edgeWagY(2:end)+edgeWagY(1:end-1));
% [ windstressLess, ~ ] = circulationInt( xm,ym,windE,windN,edgeWagXH,edgeWagYH );
% edgeWagX2=[edgeWagX(1:end-1);edgeWagXH]; edgeWagX2=edgeWagX2(:);
% edgeWagY2=[edgeWagY(1:end-1);edgeWagYH]; edgeWagY2=edgeWagY2(:);
% [ windstressTwice, ~ ] = circulationInt( xm,ym,windE,windN,edgeWagX2.',edgeWagY2.' );
% [ wind100, ~] = circulationInt( xm,ym,windE,windN,x100,y100 );
% [ wind250, ~] = circulationInt( xm,ym,windE,windN,x250,y250 );
% [ wind500, ~] = circulationInt( xm,ym,windE,windN,x500,y500 );
% [ wind100c, ~] = circulationInt3( xm,ym,windE,windN,x100,y100 );
% [ wind250c, ~] = circulationInt3( xm,ym,windE,windN,x250,y250 );
% [ wind500c, ~] = circulationInt3( xm,ym,windE,windN,x500,y500 );
%[corr1,p1]=corrcoef(dCircdt,windstress(1:end-1))
%% calculate diffusion (viscosity) component
load('distancesAreas')
dZ=5;
rAw=reshape(raw,[700 200]);
rAs=reshape(ras,[700 200]);
cellVolU=rAw.*hFacW(:,:,1).*dZ;
cellVolV=rAs.*hFacS(:,:,1).*dZ;
% load('momentumPressureSurface.mat','SSHa')
% [~,~,nt]=size(SSHa);
% cellVolU2=repmat(rAw.*hFacW(:,:,1),[1 1 nt]).*(dZ+SSHa);
% cellVolV2=repmat(rAs.*hFacS(:,:,1),[1 1 nt]).*(dZ+SSHa);
%  load('momentumDiagnostics148dayNF2.mat','VisZU')
%  VisZU=VisZU(:,:,1:2,:);
%  [~,~,~,nt]=size(VisZU);
%  UDif2a=(VisZU(:,:,2:end,:)-VisZU(:,:,1:end-1,:))./repmat(cellVolU(:,:,1),[1 1 1 nt]);
%  UDif2a=squeeze(UDif2a(:,:,1,:)); UDif2b=squeeze(VisZU(:,:,2:end,:)-VisZU(:,:,1:end-1,:))./cellVolU2;
% clear VisZU
%  disp('UDif2a done')
%  load('momentumDiagnostics148dayNF2.mat','VisZV')
%  VisZV=VisZV(:,:,1:2,:);
%  VDif2a=(VisZV(:,:,2:end,:)-VisZV(:,:,1:end-1,:))./repmat(cellVolV(:,:,1),[1 1 1 nt]);
%  VDif2a=squeeze(VDif2a(:,:,1,:));VDif2b=squeeze(VisZV(:,:,2:end,:)-VisZV(:,:,1:end-1,:))./cellVolV2;
%  clear VisZU
%  disp('VDif2a done')
%  save('momentumDiffusionSurface.mat','*Dif2*','cellVol*')
%  disp('save done')
%% calculate diffusion (viscosity) component
load('momentumDiffusionSurface.mat')
%remove nans (where there is land)
UDif2a(isnan(UDif2a))=0;
VDif2a(isnan(VDif2a))=0;
UDif2b(isnan(UDif2b))=0;
VDif2b(isnan(VDif2b))=0;
%rotate to east-north
diffE=UDif2a.*repmat(AngleCS,[1 1 nt]) - VDif2a.*repmat(AngleSN,[1 1 nt]);  
diffN=UDif2a.*repmat(AngleSN,[1 1 nt]) + VDif2a.*repmat(AngleCS,[1 1 nt]); 
[ diffusion, sign3 ] = circulationInt( xum,yum,diffE,xvm,yvm,diffN,edgeWagX,edgeWagY );
%[ diff100, ~] = circulationInt( xm,ym,diffE,diffN,x100,y100 );
[ diffusionc, sign3 ] = circulationInt3( xum,yum,diffE,xvm,yvm,diffN,edgeWagX,edgeWagY );
%[ diff100c, ~] = circulationInt3( xm,ym,diffE,diffN,x100,y100 );
%[ diff500, ~] = circulationInt( xm,ym,diffE,diffN,x500,y500 );
%[ diff500c, ~] = circulationInt3( xm,ym,diffE,diffN,x500,y500 );
%[ diff250, ~] = circulationInt( xm,ym,diffE,diffN,x250,y250 );
%[ diff250c, ~] = circulationInt3( xm,ym,diffE,diffN,x250,y250 );

diffE2=UDif2b.*repmat(AngleCS,[1 1 nt]) - VDif2b.*repmat(AngleSN,[1 1 nt]);  
diffN2=UDif2b.*repmat(AngleSN,[1 1 nt]) + VDif2b.*repmat(AngleCS,[1 1 nt]); 
[ diffusion2, sign30 ] = circulationInt( xum,yum,diffE2,xvm,yvm,diffN2,edgeWagX,edgeWagY );
%[ diff100b, ~] = circulationInt( xm,ym,diffE2,diffN2,x100,y100 );
[ diffusion2c, sign30 ] = circulationInt3( xum,yum,diffE2,xvm,yvm,diffN2,edgeWagX,edgeWagY );
%[ diff100bc, ~] = circulationInt3( xm,ym,diffE2,diffN2,x100,y100 );
%[ diff500b, ~] = circulationInt( xm,ym,diffE2,diffN2,x500,y500 );
%[ diff500bc, ~] = circulationInt3( xm,ym,diffE2,diffN2,x500,y500 );
%[ diff250b, ~] = circulationInt( xm,ym,diffE2,diffN2,x250,y250 );
%[ diff250bc, ~] = circulationInt3( xm,ym,diffE2,diffN2,x250,y250 );

%[corr2,p2]=corrcoef(dCircdt,diffusion(1:end-1))
%[corr3,p3]=corrcoef(windstress,diffusion)
%[corr4,p4]=corrcoef(dCircdt,diffusion(1:end-1)+windstress(1:end-1))
%figure; plot(dCircdt); hold all; plot(windstress); plot(diffusion); plot(diffusion2)
%% dissipation
load('momentumAdvDissSurface.mat','*Diss')
%rotate to east-north
dissE=UDiss.*repmat(AngleCS,[1 1 nt]) - VDiss.*repmat(AngleSN,[1 1 nt]);  
dissN=UDiss.*repmat(AngleSN,[1 1 nt]) + VDiss.*repmat(AngleCS,[1 1 nt]); 
[ dissipation, sign4 ] = circulationInt( xum,yum,dissE,xvm,yvm,dissN,edgeWagX,edgeWagY );
[ dissipationc, ~ ] = circulationInt3( xum,yum,dissE,xvm,yvm,dissN,edgeWagX,edgeWagY );
%[ diss100, ~] = circulationInt( xm,ym,dissE,dissN,x100,y100 );
%[ diss250, ~] = circulationInt( xm,ym,dissE,dissN,x250,y250 );
%[ diss500, ~] = circulationInt( xm,ym,dissE,dissN,x500,y500 );
%[ diss100c, ~] = circulationInt3( xm,ym,dissE,dissN,x100,y100 );
%[ diss250c, ~] = circulationInt3( xm,ym,dissE,dissN,x250,y250 );
%[ diss500c, ~] = circulationInt3( xm,ym,dissE,dissN,x500,y500 );

%% timestep
load('momentumTendAbSurface.mat','Ab*')
%rotate to east-north
abE=AbU.*repmat(AngleCS,[1 1 nt]) - AbV.*repmat(AngleSN,[1 1 nt]);  
abN=AbU.*repmat(AngleSN,[1 1 nt]) + AbV.*repmat(AngleCS,[1 1 nt]); 
[ timestep, sign5 ] = circulationInt( xum,yum,abE,xvm,yvm,abN,edgeWagX,edgeWagY );
%[ time100, ~] = circulationInt( xm,ym,abE,abN,x100,y100 );
%[ time250, ~] = circulationInt( xm,ym,abE,abN,x250,y500 );
%[ time500, ~] = circulationInt( xm,ym,abE,abN,x250,y500 );
[ timestepc, ~ ] = circulationInt3( xum,yum,abE,xvm,yvm,abN,edgeWagX,edgeWagY );
% [ time100c, ~] = circulationInt3( xm,ym,abE,abN,x100,y100 );
% [ time250c, ~] = circulationInt3( xm,ym,abE,abN,x250,y500 );
% [ time500c, ~] = circulationInt3( xm,ym,abE,abN,x250,y500 );
%% advection
load('momentumAdvDissSurface.mat','Adv*')
%rotate to east-north
advE=AdvU.*repmat(AngleCS,[1 1 nt]) - AdvV.*repmat(AngleSN,[1 1 nt]);  
advN=AdvU.*repmat(AngleSN,[1 1 nt]) + AdvV.*repmat(AngleCS,[1 1 nt]); 
[ advection, sign6 ] = circulationInt( xum,yum,advE,xvm,yvm,advN,edgeWagX,edgeWagY );
[ advectionc, ~ ] = circulationInt3( xum,yum,advE,xvm,yvm,advN,edgeWagX,edgeWagY );
%[ adv100, ~] = circulationInt( xm,ym,advE,advN,x100,y100 );
% [ adv250, ~] = circulationInt( xm,ym,advE,advN,x250,y250 );
% [ adv500, ~] = circulationInt( xm,ym,advE,advN,x500,y500 );
% [ adv100c, ~] = circulationInt3( xm,ym,advE,advN,x100,y100 );
% [ adv250c, ~] = circulationInt3( xm,ym,advE,advN,x250,y250 );
% [ adv500c, ~] = circulationInt3( xm,ym,advE,advN,x500,y500 );

%% calculate presure field
load('momentumPressureSurface.mat')
pressUssh=Utend/86400-AdvU-AbU-UDiss-UDif2a-Uext-UdPdx;
pressVssh=Vtend/86400-AdvV-AbV-VDiss-VDif2a-Vext-VdPdy;

pressUssh2=Utend/86400-AdvU-AbU-UDiss-UDif2b-Uext-UdPdx;
pressVssh2=Vtend/86400-AdvV-AbV-VDiss-VDif2b-Vext-VdPdy;
%rotate to east-north
press1E=UdPdx.*repmat(AngleCS,[1 1 nt]) - VdPdy.*repmat(AngleSN,[1 1 nt]);  
press1N=UdPdx.*repmat(AngleSN,[1 1 nt]) + VdPdy.*repmat(AngleCS,[1 1 nt]); 

press2E=pressUssh.*repmat(AngleCS,[1 1 nt]) - pressVssh.*repmat(AngleSN,[1 1 nt]);  
press2N=pressUssh.*repmat(AngleSN,[1 1 nt]) + pressVssh.*repmat(AngleCS,[1 1 nt]); 

pressE=(pressUssh+UdPdx).*repmat(AngleCS,[1 1 nt]) - (pressVssh+VdPdy).*repmat(AngleSN,[1 1 nt]);  
pressN=(pressUssh+UdPdx).*repmat(AngleSN,[1 1 nt]) + (pressVssh+VdPdy).*repmat(AngleCS,[1 1 nt]); 

pressE2=(pressUssh2+UdPdx).*repmat(AngleCS,[1 1 nt]) - (pressVssh2+VdPdy).*repmat(AngleSN,[1 1 nt]);  
pressN2=(pressUssh2+UdPdx).*repmat(AngleSN,[1 1 nt]) + (pressVssh2+VdPdy).*repmat(AngleCS,[1 1 nt]); 
% pressure terms
[ pressureTrop, sign7 ] = circulationInt( xum,yum,press2E,xvm,yvm,press2N,edgeWagX,edgeWagY );
[ pressureClin, sign9 ] = circulationInt( xum,yum,press1E,xvm,yvm,press1N,edgeWagX,edgeWagY );
%[ pclin100, ~] = circulationInt( xm,ym,press1E,press1N,x100,y100 );
%[ pclin250, ~] = circulationInt( xm,ym,press1E,press1N,x250,y250 );
%[ pclin500, ~] = circulationInt( xm,ym,press1E,press1N,x500,y500 );
[ pressureClinc, ~ ] = circulationInt3( xum,yum,press1E,xvm,yvm,press1N,edgeWagX,edgeWagY );
%[ pclin100c, ~] = circulationInt3( xm,ym,press1E,press1N,x100,y100 );
%[ pclin250c, ~] = circulationInt3( xm,ym,press1E,press1N,x250,y250 );
%[ pclin500c, ~] = circulationInt3( xm,ym,press1E,press1N,x500,y500 );

[ pressure, sign10 ] = circulationInt( xum,yum,pressE,xvm,yvm,pressN,edgeWagX,edgeWagY );
[ pressure2, sign20 ] = circulationInt( xum,yum,pressE2,xvm,yvm,pressN2,edgeWagX,edgeWagY );%negligible differences!
%%


%% other pressure calculations
% DXC=reshape(dxc,[700 200]);
% DYC=reshape(dyc,[700 200]);
% [~,~,nt]=size(sshSnap)
% dhdx=zeros(size(sshSnap));
% dhdy=dhdx;
% g=9.81;
% dhdx(2:end,:,:)=-g*(sshSnap(2:end,:,:)-sshSnap(1:end-1,:,:))./repmat(DXC(1:end-1,:),[1 1 nt]);
% dhdy(:,2:end,:)=-g*(sshSnap(:,2:end,:)-sshSnap(:,1:end-1,:))./repmat(DYC(:,1:end-1),[1 1 nt]);
% press3E=(UdPdx+dhdx).*repmat(AngleCS,[1 1 nt]) - (VdPdy+dhdy).*repmat(AngleSN,[1 1 nt]);  
% press3N=(UdPdx+dhdx).*repmat(AngleSN,[1 1 nt]) + (VdPdy+dhdy).*repmat(AngleCS,[1 1 nt]); 
% 
% [ pressure3, sign11 ] = circulationInt( xm,ym,press3E,press3N,edgeWagX,edgeWagY );%much too big

%% tendency (rather than change in daily-average circ, this would give daily-average change in circ)
load('momentumTendAbSurface.mat','*tend')
nt=148;
%rotate to east-north
tendE=Utend.*repmat(AngleCS,[1 1 nt])./86400 - Vtend.*repmat(AngleSN,[1 1 nt])./86400;  
tendN=Utend.*repmat(AngleSN,[1 1 nt])./86400 + Vtend.*repmat(AngleCS,[1 1 nt])./86400; 
[ tendency, sign8 ] = circulationInt( xum,yum,tendE,xvm,yvm,tendN,edgeWagX,edgeWagY );
disp('1')
[ tendencyC, sign8 ] = circulationInt3( xum,yum,tendE,xvm,yvm,tendN,edgeWagX,edgeWagY );
% [ tend1000, ~ ] = circulationInt( xm,ym,tendE,tendN,x1000,y1000 );
% [ tend2000, ~ ] = circulationInt( xm,ym,tendE,tendN,x2000,y2000 );
disp('first 3 done')
%[ tend100, ~] = circulationInt( xm,ym,tendE,tendN,x100,y100 );
%[ tend250, ~ ] = circulationInt( xm,ym,tendE,tendN,x250,y250 );
disp('first 5 done')
%[ tend500, ~ ] = circulationInt( xm,ym,tendE,tendN,x500,y500 );

%figure; plot(tendency-tend2000); hold all; plot(tend2000-tend1000); plot(tend1000-tend500); plot(tend500-tend250); plot(tend250-tend100)
%%
% [ tendencyN, sign8 ] = circulationInt( xm,ym,tendE,tendN,edgeWagX,edgeWagY );
% disp('1')
% [ tend1000N, ~ ] = circulationInt( xm,ym,tendE,tendN,x1000,y1000 );
% [ tend2000N, ~ ] = circulationInt( xm,ym,tendE,tendN,x2000,y2000 );
% disp('first 3 done')
% [ tend100N, ~] = circulationInt( xm,ym,tendE,tendN,x100,y100 );
% [ tend250N, ~ ] = circulationInt( xm,ym,tendE,tendN,x250,y250 );
% disp('first 5 done')
% [ tend500N, ~ ] = circulationInt( xm,ym,tendE,tendN,x500,y500 );

% figure; plot(tendency-tend2000); hold all; plot(tend2000-tend1000); plot(tend1000-tend500); plot(tend500-tend250); plot(tend250-tend100)
%%
% [ pressureClinN, sign9 ] = circulationInt( xm,ym,press1E,press1N,edgeWagX,edgeWagY );
% [ pclin100N, ~] = circulationInt( xm,ym,press1E,press1N,x100,y100 );
% [ pclin500N, ~] = circulationInt( xm,ym,press1E,press1N,x500,y500 );
% disp('nearest done')
% [ tendencyNa, ~ ] = circulationInt2( xm,ym,tendE,tendN,edgeWagX,edgeWagY );
% disp('1')
%%
% [ tend1000Na, ~ ] = circulationInt2( xm,ym,tendE,tendN,x1000,y1000 );
% [ tend2000Na, ~ ] = circulationInt2( xm,ym,tendE,tendN,x2000,y2000 );
% disp('first 3 done')
% [ tend100Na, ~] = circulationInt2( xm,ym,tendE,tendN,x100,y100 );
% [ tend250Na, ~ ] = circulationInt2( xm,ym,tendE,tendN,x250,y250 );
% disp('first 5 done')
% [ tend500Na, ~ ] = circulationInt2( xm,ym,tendE,tendN,x500,y500 );
% 
% [ pressureClinNa, ~ ] = circulationInt2( xm,ym,press1E,press1N,edgeWagX,edgeWagY );
% [ pclin100Na, ~] = circulationInt2( xm,ym,press1E,press1N,x100,y100 );
% [ pclin500Na, ~] = circulationInt2( xm,ym,press1E,press1N,x500,y500 );
% disp('natural done')
%[ pressureClinC, sign9 ] = circulationInt3( xm,ym,press1E,press1N,edgeWagX,edgeWagY );
%[ pclin100C, ~] = circulationInt3( xm,ym,press1E,press1N,x100,y100 );
%[ pclin250C, ~] = circulationInt3( xm,ym,press1E,press1N,x250,y250 );
%[ pclin500C, ~] = circulationInt3( xm,ym,press1E,press1N,x500,y500 );

%[ tendencyC, ~ ] = circulationInt3( xm,ym,tendE,tendN,edgeWagX,edgeWagY );
disp('1')
%[ tend1000C, ~ ] = circulationInt3( xm,ym,tendE,tendN,x1000,y1000 );
%[ tend2000C, ~ ] = circulationInt3( xm,ym,tendE,tendN,x2000,y2000 );
disp('first 3 done')
%[ tend100C, ~] = circulationInt3( xm,ym,tendE,tendN,x100,y100 );
%[ tend250C, ~ ] = circulationInt3( xm,ym,tendE,tendN,x250,y250 );
disp('first 5 done')
%[ tend500C, ~ ] = circulationInt3( xm,ym,tendE,tendN,x500,y500 );
disp('cubic done, saving')
%save('circulationIntegralTesting.mat')
%% coriolis term (part of advection term above)
load('momentumCoriSurface.mat')
nt=148;
UCori=UCori(:,:,1:148);
VCori=VCori(:,:,1:148);
coriE=UCori.*repmat(AngleCS,[1 1 nt])./86400 - VCori.*repmat(AngleSN,[1 1 nt])./86400;  
coriN=UCori.*repmat(AngleSN,[1 1 nt])./86400 + VCori.*repmat(AngleCS,[1 1 nt])./86400; 
[ coriolis, sign12 ] = circulationInt( xum,yum,coriE,xvm,yvm,coriN,edgeWagX,edgeWagY );
[ coriolisc, ~ ] = circulationInt3( xum,yum,coriE,xvm,yvm,coriN,edgeWagX,edgeWagY );
% [ cori100, ~ ] = circulationInt( xm,ym,coriE,coriN,x100,y100 );
% [ cori250, ~ ] = circulationInt( xm,ym,coriE,coriN,x250,y250 );
% [ cori500, ~ ] = circulationInt( xm,ym,coriE,coriN,x500,y500 );
% [ cori100c, ~ ] = circulationInt3( xm,ym,coriE,coriN,x100,y100 );
% [ cori250c, ~ ] = circulationInt3( xm,ym,coriE,coriN,x250,y250 );
% [ cori500c, ~ ] = circulationInt3( xm,ym,coriE,coriN,x500,y500 );
%% does it close with area-integrated vorticity changes instead?
nt=148;
load('vorticitySurface.mat')
load('distancesAreas')
rAz=reshape(raz,[700 200]);
dvortdt=diff(squeeze(nansum(nansum(vort(:,:,1:149).*repmat(rAz.*inWag(:,:,1),[1 1 149])))))./86400;
dvortadt=diff(squeeze(nansum(nansum(vorta(:,:,1:149).*repmat(rAz.*inWag(:,:,1),[1 1 149])))))./86400;
figure; plot(tendency); hold all; plot(dvortdt); plot(dvortadt)
%% plots
% figure; plot(tendency,'linewidth',3); hold all; plot(windstress,'linewidth',2);
% plot(diffusion,'linewidth',2);plot(dissipation,'linewidth',2);plot(advection,'linewidth',2);plot(pressureClin,'linewidth',2);plot(timestep);
% legend('d/dt','wind','vertical diffusion','dissipation','advection','pressure','timestep')
% 
% figure; plot(tend100,'linewidth',3); hold all; plot(windstress100,'linewidth',2);
% plot(diff100,'linewidth',2);plot(diss100,'linewidth',2);plot(adv100,'linewidth',2);plot(pclin100,'linewidth',2);plot(time100);
% legend('d/dt','wind','vertical diffusion','dissipation','advection','pressure','timestep'); title('New')
% 
% figure; plot(-tendency+windstress+dissipation+advection+timestep+pressureClin+diffusion,'--','linewidth',2)
% hold all; plot(-tend100+windstress100+diss100+adv100+time100+pclin100+diff100)
% 
% 
figure; plot(-tendency,'linewidth',2); hold all; plot(wind,'linewidth',2);plot(diffusion,'linewidth',2);plot(dissipation,'linewidth',2);
plot(advection,'linewidth',2);plot(pressureClin,'linewidth',2);plot(timestep,'linewidth',2);
plot(-tendency+wind+dissipation+advection+timestep+pressureClin+diffusion,'--','linewidth',2)
legend('-d/dt','wind','vertical diffusion','dissipation','advection','pressure','timestep','total')

figure; plot(-tendency+wind+dissipation+advection+timestep+pressureClin+diffusion,'--','linewidth',2)
hold all; plot(-dvortdt+wind+dissipation+advection+timestep+pressureClin+diffusion,'linewidth',2)
plot(-dvortadt+windstress+dissipation+advection+timestep+pressureClin+diffusion,'linewidth',2)
% 
% figure; plot(-tendency,'linewidth',2); hold all; plot(windstress+dissipation,'linewidth',2);plot(diffusion+advection,'linewidth',2);
% plot(pressure,'linewidth',2);plot(timestep,'linewidth',2);
% plot(-tendency+windstress+dissipation+advection+timestep+pressure+diffusion,'--','linewidth',2)
% legend('-d/dt','wind+dissipation','vertical diffusion+advection','pressure','timestep','total')
% 
% figure; plot(-tendency,'linewidth',2); hold all; plot(windstress+diffusion,'linewidth',2);plot(dissipation+advection,'linewidth',2);
% plot(pressure,'linewidth',2);plot(timestep,'linewidth',2);
% plot(-tendency+windstress+dissipation+advection+timestep+pressure+diffusion,'--','linewidth',2)
% legend('-d/dt','wind+diffusion','dissipation+advection','pressure','timestep','total')
% 
% figure; plot(tendency,'linewidth',2); hold all; plot(windstress+dissipation+diffusion+advection,'linewidth',2);
% plot(pressure,'linewidth',2);plot(timestep,'linewidth',2);
% legend('d/dt','wind+dissipation+vertical+diffusion+advection','pressure','timestep')

%figure; plot(tendency,'linewidth',2); hold all; plot(dCircdt)
%plot(windstress+dissipation+advection+timestep+pressure+diffusion)
%plot(windstress+dissipation+advection+timestep+diffusion)
%legend('curl of d/dt','d/dt of circ','forces sum','except pressure')

%figure; plot(tendency-(windstress+dissipation+advection+timestep+pressure)); hold all; plot(diffusion)
%plot(tendency-(windstress+dissipation+advection+timestep+pressure+diffusion));
% %%
% figure; plot(windstress./max(windstress)); hold all; 
% plot(diffusion./max(diffusion))
% plot(diff(integral1)./max(diff(integral1)));
%% means and correlations
disp('last section')
allTerms=[tendency wind diffusion dissipation advection pressureClin timestep];

[rhos,ps]=corr(allTerms);
means=mean(allTerms,1);
meanMag=mean(abs(allTerms));

% allTerms100=[tend100 wind100 diff100 diss100 adv100 pclin100 time100];
% 
% [rhos100,ps100]=corr(allTerms100);
% means100=mean(allTerms100,1);
% meanMag100=mean(abs(allTerms100));

clear *E *N *dPd* *rot *ext *Dif* Adv* Ab* *Cori Utend Vtend
%save('circulationIntegralForErrbars.mat')
disp('done')
%% plots of different methods etc
figure; plot(wind500); hold all; plot(wind500c); plot(wind100); plot(wind100c)
%% spread of each value- 500 to 100, linear and cubic
windset=[wind windc];%wind100 windc wind100c];
windM=mean(windset,2);
windS=std(windset,0,2);
windU=max(windset,[],2)-windM;
windL=windM-min(windset,[],2);
tendset=[tendency tendencyC diff(vortaInt3)/86400 diff(vortInt3)/86400];%tend100 tendencyC tend100C];
tendM=mean(tendset,2);
tendS=std(tendset,0,2);
tendU=max(-tendset,[],2)-(-tendM);
tendL=-tendM-min(-tendset,[],2);
dissset=[dissipation dissipationc];
dissM=mean(dissset,2);
dissS=std(dissset,0,2);
dissU=max(dissset,[],2)-dissM;
dissL=dissM-min(dissset,[],2);
advset=[advection advectionc];
advM=mean(advset,2);
advS=std(advset,0,2);
advU=max(advset,[],2)-advM;
advL=advM-min(advset,[],2);
timeset=[timestep timestepc];
timeM=mean(timeset,2);
timeS=std(timeset,0,2);
timeU=max(timeset,[],2)-timeM;
timeL=timeM-min(timeset,[],2);
clinset=[pressureClin pressureClinc];
clinM=mean(clinset,2);
clinS=std(clinset,0,2);
clinU=max(clinset,[],2)-clinM;
clinL=clinM-min(clinset,[],2);
diffset=[diffusion diffusionc diffusion2 diffusion2c];
diffM=mean(diffset,2);
diffS=std(diffset,0,2);
diffU=max(diffset,[],2)-diffM;
diffL=diffM-min(diffset,[],2);
%% does it close
totM=-tendM+windM+dissM+advM+timeM+clinM+diffM;
totS=tendS+windS+dissS+advS+timeS+clinS+diffS;
figure; plot(totS./totM)
figure; plot(-tendM+windM+dissM+advM+timeM+clinM+diffM)
hold on; errorbar(-tendM+windM+dissM+advM+timeM+clinM+diffM,tendS+windS+dissS+advS+timeS+clinS+diffS)
hold all; errorbar(-tendM+windM+dissM+advM+timeM+clinM+diffM,abs(tendency-diff(vortaInt3)/86400))

figure; errorbar(-tendM+windM+dissM+advM+timeM+clinM+diffM,tendS+windS+dissS+advS+timeS+clinS+diffS)
hold all; errorbar(-diff(vortaInt3)/86400+windM+dissM+advM+timeM+clinM+diffM,tendS+windS+dissS+advS+timeS+clinS+diffS)

err1=-tendM+windM+dissM+advM+timeM+clinM+diffM;
errU=tendS+windS+dissS+advS+timeS+clinS+diffS;
errL=(tendS+windS+dissS+advS+timeS+clinS+diffS);
errU2=tendU+windU+dissU+advU+timeU+clinU+diffU;
errL2=(tendL+windL+dissL+advL+timeL+clinL+diffL);
%nClosed=sum(errU2>0 &errL2<0)%35, 38
%%
L1=tendL+windL+dissL+advL+timeL+clinL+diffL;
U1=tendU+windU+dissU+advU+timeU+clinU+diffU;
figure; plot(tendM); hold all; plot(windM); plot(diffM); plot(dissM); plot(advM); plot(clinM); plot(timeM);
plot(windM+dissM+advM+timeM+clinM+diffM,'--','linewidth',2)
hold on; errorbar(1:148,windM+dissM+advM+timeM+clinM+diffM,windU+dissU+advU+timeU+clinU+diffU,windL+dissL+advL+timeL+clinL+diffL)
legend('d/dt','wind','vertical diff','horiz diff and drag','adv and cori','baroclinic','timestep','sum of forces','error')

[rho1,p1]=corrcoef(tendM,windM+dissM+advM+timeM+clinM+diffM)
[rho2,p2]=corrcoef(tendency,wind+dissipation+advection+timestep+pressureClin+diffusion)%best!
[rho3,p3]=corrcoef(tendencyC,windc+dissipationc+advectionc+timestepc+pressureClinC+diffusionc)
%% does it close
figure; plot(-tendency+wind+dissipation+advection+timestep+pressureClin+diffusion2)
hold on; errorbar(-tendency+wind+dissipation+advection+timestep+pressureClin+diffusion2,tendS+windS+dissS+advS+timeS+clinS+diffS,'.')

errU3=-tendency+wind+dissipation+advection+timestep+pressureClin+diffusion2+tendU+windU+dissU+advU+timeU+clinU+diffU;
errL3=-tendency+wind+dissipation+advection+timestep+pressureClin+diffusion2-(tendL+windL+dissL+advL+timeL+clinL+diffL);
nClosed=sum(errU3>0 &errL3<0)%50
%%
isClosed=errU3>0 &errL3<0;
err3=-tendency+wind+dissipation+advection+timestep+pressureClin+diffusion;
days=1:148;
figure; plot(err3)
hold on; errorbar(err3,tendS+windS+dissS+advS+timeS+clinS+diffS,'.')
hold on; plot(days(isClosed),err3(isClosed),'go')

%%
figure; plot(windL./abs(wind)); hold all; plot(windU./abs(wind)); plot(windS./abs(wind))
figure; plot(diffL./abs(diffusion)); hold all; plot(diffU./abs(diffusion)); plot(diffS./abs(diffusion))

errEst=0.05*abs(tendency)+0.05*abs(wind)+0.05*abs(dissipation)+0.05*abs(advection)+0.05*abs(timestep)+0.05*abs(pressureClin)+0.05*abs(diffusion2);
errU4=-tendency+wind+dissipation+advection+timestep+pressureClin+diffusion2+errEst;
errL4=-tendency+wind+dissipation+advection+timestep+pressureClin+diffusion2-(errEst);
nClosed=sum(errU4>0 &errL4<0)%49
errU5=max(errU3,errU4);
errL5=min(errL3,errL4);
nClosed=sum(errU5>0 &errL5<0)%62
%%
isClosed=errU5>0 &errL5<0;
days=1:148;
tot1=-tendency+wind+dissipation+advection+timestep+pressureClin+diffusion;
figure; plot(-tendency,'linewidth',2); hold all; plot(wind,'linewidth',2);plot(diffusion,'linewidth',2);plot(dissipation,'linewidth',2);
plot(advection,'linewidth',2);plot(pressureClin,'linewidth',2);plot(timestep,'linewidth',2);
errorbar(1:148,tot1,tot1-errL5,errU5-tot1,'--','linewidth',2)
plot(days(isClosed),err3(isClosed),'go')
legend('-d/dt','wind','vertical diffusion','dissipation','advection','pressure','timestep','total with errorbars','closed within error')

%% perror
forcetot=wind+dissipation+advection+timestep+pressureClin+diffusion;

perror=100*(tendency-forcetot)./abs(tendency); figure; plot(abs(perror))
%%
forcetot2=wind+dissipation+advection+timestep+pressure+diffusion;
figure; plot(forcetot-tendency); hold all; plot(pressureTrop); plot(forcetot2-tendency)















