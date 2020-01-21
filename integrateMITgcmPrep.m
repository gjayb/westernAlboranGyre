load('geometrySpinupSteady.mat');
%clear Urot Vrot
load('uvwDailyDepth1rotatedNF.mat');
u1=Urot; v1=Vrot; 
load('uvwDailyNativeNF.mat','W');
w1=squeeze(W(:,:,1,:)); clear W

%u1=squeeze(Urot);
%v1=squeeze(Vrot);

%u2=mean(Urot(:,:,122:152),3);
%v2=mean(Vrot(:,:,122:152),3);

[~,~,nt]=size(u1);
tvel=0:86400:(nt-1)*86400;

xmin=min(min(XC)); ymin=min(min(YC));
NX=150; NY=90;
    xmin2=-7.5;
    xmax2=1;
    ymin2=34.5;
    ymax2=37.5;

    xin=linspace(xmin2,xmax2,NX);
    yin=linspace(ymin2,ymax2,NY);

    [xin2, yin2]=meshgrid(xin,yin);
    %convert to meters from lat/lon; lat is y, lon is x
    xinM=(xin2-xmin*ones(size(xin2))).*111000.*cosd(yin2); 
    yinM=(yin2-ymin*ones(size(yin2))).*111000;
    
xvel2=(XC-xmin*ones(size(XC))).*111000.*cosd(YC); 
yvel2=(YC-ymin*ones(size(YC))).*111000;
%SHOULD have XU,YU for u and XV,YV for v, w uses *C 

xvel=min(min(xvel2)):1000:max(max(xvel2));
yvel=min(min(yvel2)):1000:max(max(yvel2));

sizexvel=size(xvel)
sizexvel2=size(xvel2)

[xvelg,yvelg]=meshgrid(xvel,yvel);
sizexvelg=size(xvelg)
for k=1:length(tvel)
    u3(:,:,k)=griddata(xvel2,yvel2,u1(:,:,k),xvelg,yvelg);
    v3(:,:,k)=griddata(xvel2,yvel2,v1(:,:,k),xvelg,yvelg);
   w3(:,:,k)=griddata(xvel2,yvel2,w1(:,:,k),xvelg,yvelg);
end
u=u3; v=v3; w=w3;
clear Urot Vrot W u1 v1 u2 v2 u3 v3 w1 w3 
%u=repmat(u3,[1 1 32]);
%v=repmat(v3,[1 1 32]);

z0=[xinM(:)';yinM(:)'];

%fn='uvHdepth1interpolated2.mat';
fn='uvwDailyDepth1interpolatedNF.mat';
save(fn,'-v7.3');
