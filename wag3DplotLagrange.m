%3D WAG plot with isopycnals' manifolds
%load
load('geometrySpinupSteady')
load('isoDepthsNF')
load('sigma2675ClosedNFall.mat', 'lonWagClosed')
load('sigma2675ClosedNFall.mat', 'latWagClosed')
np=find(latWagClosed(i,:)>0,1,'last');
lon2675=lonWagClosed(i,1:np);
lat2675=latWagClosed(i,1:np);
load('sigma275ClosedNF1to33.mat', 'lonWagClosed')
load('sigma275ClosedNF1to33.mat', 'latWagClosed')
np=find(latWagClosed(i,:)>0,1,'last');
lon275=lonWagClosed(i,1:np);
lat275=latWagClosed(i,1:np);
load('sigma263ClosedNF1to119.mat', 'lonWagClosed')
load('sigma263ClosedNF1to119.mat', 'latWagClosed')
np=find(latWagClosed(i,:)>0,1,'last');
lon263=lonWagClosed(i,1:np);
lat263=latWagClosed(i,1:np);
load('sigma27ClosedNF.mat', 'lonWagClosed')
load('sigma27ClosedNF.mat', 'latWagClosed')
np=find(latWagClosed(i,:)>0,1,'last');
lon27=lonWagClosed(i,1:np);
lat27=latWagClosed(i,1:np);
load('sigma265ClosedNF.mat')
np=find(latWagClosed(i,:)>0,1,'last');
lon265=lonWagClosed(i,1:np);
lat265=latWagClosed(i,1:np);
load('surfaceClosedNF.mat')
np=find(latWagClosed(i,:)>0,1,'last');
lonS=lonWagClosed(i,1:np);
latS=latWagClosed(i,1:np);
%% get z,d 
z263=griddata(XC,YC,-isoDepth263(:,:,20),lon263,lat263);
z265=griddata(XC,YC,-isoDepth265(:,:,20),lon265,lat265);
z2675=griddata(XC,YC,-isoDepth2675(:,:,20),lon2675,lat2675);
z27=griddata(XC,YC,-isoDepth27(:,:,20),lon27,lat27);
z275=griddata(XC,YC,-isoDepth275(:,:,20),lon275,lat275);

d275=griddata(XC,YC,d,lon275,lat275);
d27=griddata(XC,YC,d,lon27,lat27);
d2675=griddata(XC,YC,d,lon2675,lat2675);
d265=griddata(XC,YC,d,lon265,lat265);
d263=griddata(XC,YC,d,lon263,lat263);
dS=griddata(XC,YC,d,lonS,latS);

xS=lonS(dS>40);
yS=latS(dS>40);
x263=lon263(d263>40);
y263=lat263(d263>40);
z263b=z263(d263>40);
x265=lon265(d265>40);
y265=lat265(d265>40);
z265b=z265(d265>40);
x2675=lon2675(d2675>40);
y2675=lat2675(d2675>40);
z2675b=z2675(d2675>40);
x27=lon27(d27>40);
y27=lat27(d27>40);
z27b=z27(d27>40);
x275=lon275(d275>40);
y275=lat275(d275>40);
z275b=z275(d275>40);

%% prep to plot iso
hold275=-isoDepth275(:,:,20);
indices275=inpolygon(XC,YC,lon275,lat275);
zmesh275=zeros([700 200]);
zmesh275(:)=nan;
zmesh275(indices275)=hold275(indices275);
%% plot 3D with iso275
figure; hsurf=surf(XC,YC,zmesh275); shading 'flat'
alpha(0.5)
hold all; plot3(x275,y275,-z275b,'linewidth',2)
hold all; plot3(x27,y27,-z27b,'linewidth',2)
hold all; plot3(x2675,y2675,-z2675b,'linewidth',2)
hold all; plot3(x265,y265,-z265b,'linewidth',2)
hold all; plot3(x263,y263,-z263b,'linewidth',2)
hold all; plot(lonS,latS,'linewidth',2)
axis([-5.6 -3 35.1 36.4 -160 0])
set(gca,'fontsize',12)
xlabel('longitude'); ylabel('latitude'); zlabel('depth')
view(-111,26)
%% prep to plot curves
xunst=xS;
yunst=yS;
tunst=(1:1:length(xunst)); t1unst=linspace(1,length(xunst),200);
x0S=interp1(tunst,xunst,t1unst); 
y0S=interp1(tunst,yunst,t1unst);

xunst=x263;
yunst=y263;
zunst=z263b;
tunst=(1:1:length(xunst)); t1unst=linspace(1,length(xunst),200);
x0263=interp1(tunst,xunst,t1unst); 
y0263=interp1(tunst,yunst,t1unst);
z0263=interp1(tunst,zunst,t1unst);

xunst=x265;
yunst=y265;
zunst=z265b;
tunst=(1:1:length(xunst)); t1unst=linspace(1,length(xunst),200);
x0265=interp1(tunst,xunst,t1unst); 
y0265=interp1(tunst,yunst,t1unst);
z0265=interp1(tunst,zunst,t1unst);

xunst=x2675;
yunst=y2675;
zunst=z2675b;
tunst=(1:1:length(xunst)); t1unst=linspace(1,length(xunst),200);
x02675=interp1(tunst,xunst,t1unst); 
y02675=interp1(tunst,yunst,t1unst);
z02675=interp1(tunst,zunst,t1unst); 

xunst=x27;
yunst=y27;
zunst=z27b;
tunst=(1:1:length(xunst)); t1unst=linspace(1,length(xunst),200);
x027=interp1(tunst,xunst,t1unst); 
y027=interp1(tunst,yunst,t1unst);
z027=interp1(tunst,zunst,t1unst);

xunst=x275;
yunst=y275;
zunst=z275b;
tunst=(1:1:length(xunst)); t1unst=linspace(1,length(xunst),200);
x0275=interp1(tunst,xunst,t1unst); 
y0275=interp1(tunst,yunst,t1unst);
z0275=interp1(tunst,zunst,t1unst);

xall=[x0S; x0263; x0265; x02675; x027; x0275];
yall=[y0S; y0263; y0265; y02675; y027; y0275];
zall=[zeros(1,200); z0263; z0265; z02675; z027; z0275];

x1=[x0S; x0263];
y1=[y0S; y0263];
z1=[zeros(1,200); z0263];




