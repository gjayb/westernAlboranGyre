
%% load basics
load('geometrySpinupSteady')
load('distancesAreas','hFacW','raw','hFacS','ras')
rAw=reshape(raw,[700 200]);
rAs=reshape(ras,[700 200]);
%% make pressure field
load('momentumDiagnostics148dayNF1.mat', 'AdvU')
size(AdvU)
AdvU=AdvU(:,:,1:20,:);
load('momentumDiagnostics148dayNF1.mat', 'UDiss')
UDiss=UDiss(:,:,1:20,:);
Utot=AdvU+UDiss; clear AdvU UDiss
load('momentumDiagnostics148dayNF1.mat', 'Utend')
Utend=Utend(:,:,1:20,:)/86400;
Utot=Utot-Utend; clear Utend
load('momentumDiagnostics148dayNF2.mat', 'Uext')
Uext=Uext(:,:,1,:);
Utot(:,:,1,:)=Utot(:,:,1,:)+Uext; clear Uext
load('momentumDiagnostics148dayNF3.mat', 'AbU')
AbU=AbU(:,:,1:20,:);
Utot=Utot+AbU; clear Abu
dZ(1,1,1:20)=diff(dInterface(1:21));
cellVolU=repmat(rAw,[1 1 20]).*hFacW(:,:,1:20).*repmat(dZ,[700 200 1]);
load('momentumDiagnostics148dayNF2.mat', 'VisZU')
VisZU=VisZU(:,:,1:21,:);
UDif2a=(VisZU(:,:,2:end,:)-VisZU(:,:,1:end-1,:))./repmat(cellVolU,[1 1 1 148]); 
clear VisZU
Utot=Utot+UDif2a; clear UDif2a
pressU=-Utot; clear Utot
%%
load('momentumDiagnostics148dayNF1.mat', 'AdvV')
AdvV=AdvV(:,:,1:20,:);
load('momentumDiagnostics148dayNF1.mat', 'VDiss')
VDiss=VDiss(:,:,1:20,:);
Vtot=AdvV+VDiss; clear AdvV VDiss
load('momentumDiagnostics148dayNF1.mat', 'Vtend')
Vtend=Vtend(:,:,1:20,:)/86400;
Vtot=Vtot-Vtend; clear Vtend
load('momentumDiagnostics148dayNF2.mat', 'Vext')
Vext=Vext(:,:,1,:);
Vtot(:,:,1,:)=Vtot(:,:,1,:)+Vext; clear Vext
load('momentumDiagnostics148dayNF3.mat', 'AbV')
AbV=AbV(:,:,1:20,:);
Vtot=Vtot+AbV; clear AbV
dZ(1,1,1:20)=diff(dInterface(1:21));
cellVolV=repmat(rAs,[1 1 20]).*hFacS(:,:,1:20).*repmat(dZ,[700 200 1]);
disp('start diff')
load('momentumDiagnostics148dayNF2.mat', 'VisZV')
VisZV=VisZV(:,:,1:21,:);
VDif2a=(VisZV(:,:,2:end,:)-VisZV(:,:,1:end-1,:))./repmat(cellVolV,[1 1 1 148]); 
clear VisZV
Vtot=Vtot+VDif2a; clear VDif2a
pressV=-Vtot;
clear Vtot
disp('pressV done')
%%solve for native geostrophic velocities, clear pressure
fc=2*7.2921e-5.*sind(YC);
ugeo=pressV./repmat(fc,[1 1 20 148]);
vgeo=-pressU./repmat(fc,[1 1 20 148]);
disp('geoN exists')
%% east-west geostrophic, clear native
nt=148; geoE=ugeo.*repmat(AngleCS,[1 1 20 nt]) - vgeo.*repmat(AngleSN,[1 1 20 nt]);
geoN=ugeo.*repmat(AngleSN,[1 1 20 nt]) + vgeo.*repmat(AngleCS,[1 1 20 nt]);

geoNmean=mean(geoN,4); clear geoN
geoEmean=mean(geoE,4); clear geoE
%% movie of daily geostrophic velocities over mean ones' streamlines, layers 1,9
xmin=min(min(XC));
ymin=min(min(YC));
xcm=(XC-xmin).*111111.*cosd(YC);
ycm=(YC-ymin).*111111;

xs1=3e5:8500:6.4e5;
ys1=2.3e5:8500:4.3e5;
[xs,ys]=meshgrid(xs1,ys1);
%% regrid to mesh
xv=1000:1000:max(xcm(:));
yv=0:1000:max(ycm(:));
[xvg,yvg]=meshgrid(xv,yv);
for i=1:20
    i
geoNmg(:,:,i)=griddata(xcm,ycm,geoNmean(:,:,i),xvg,yvg);
geoEmg(:,:,i)=griddata(xcm,ycm,geoEmean(:,:,i),xvg,yvg);
end
%%
indices1=false([700 200]); indices1(200:10:500,1:20:end)=1;
v = VideoWriter('geoLevel9.avi');
v.FrameRate=3;
open(v)
f1=figure('renderer','zbuffer');

for i=1:148
    clf
    h3=streamline(stream2(xvg,yvg,geoEmg(:,:,9),geoNmg(:,:,9),xs,ys,[0.1 100]));
    set(h3,'Color','k')
    hold all
    holdE=geoE(:,:,9,i);
    holdN=geoN(:,:,9,i);
    quiver([3.5e5;xcm(indices1)],[4e5;ycm(indices1)],[0.5;holdE(indices1)],[0;holdN(indices1)],'r','linewidth',2)
    hold on; contour(xcm,ycm,inWag(:,:,9),[1 1],'m','linewidth',3)
    title(num2str(i))
    legend('mean geostrophic')
    axis tight
    writeVideo(v,getframe(gcf))
    writeVideo(v,getframe(gcf))
    writeVideo(v,getframe(gcf))
    writeVideo(v,getframe(gcf))
end
close(v)

%% calculate ageostrophic velocities: rotate total vel, remove geo
load('uvwDailyNativeNF.mat', 'U')
U=U(:,:,1:20,1:148);
uageo=U-ugeo; clear U

load('uvwDailyNativeNF.mat', 'V')
V=V(:,:,1:20,1:148);
vageo=V-vgeo; clear V

ageoE=uageo.*repmat(AngleCS,[1 1 20 nt]) - vageo.*repmat(AngleSN,[1 1 20 nt]);
ageoN=uageo.*repmat(AngleSN,[1 1 20 nt]) + vageo.*repmat(AngleCS,[1 1 20 nt]);

%% movie of daily ageo velocities over geo streamlines, layers 1,9
indices1=false([700 200]); indices1(240:10:470,1:20:end)=1;
v = VideoWriter('ageoSurf.avi');
v.FrameRate=3;
open(v)
f1=figure('renderer','zbuffer');

for i=1:148
    clf
    h3=streamline(stream2(xvg,yvg,geoEmg(:,:,1),geoNmg(:,:,1),xs,ys,[0.1 100]));
    set(h3,'Color','k')
    hold all
    %contour(xcm,ycm,inWag(:,:,1),[1 1],'m','linewidth',3)
    holdE=ageoE(:,:,1,i);
    holdN=ageoN(:,:,1,i);
    quiver([3.5e5;xcm(indices1)],[4e5;ycm(indices1)],[0.25/0.25;holdE(indices1)./0.25],[0;holdN(indices1)./0.25],'r','linewidth',2)
    title(num2str(i))
    legend('mean geostrophic')%,'salinity minimum')
    axis tight
    writeVideo(v,getframe(gcf))
    writeVideo(v,getframe(gcf))
    writeVideo(v,getframe(gcf))
    writeVideo(v,getframe(gcf))
end
close(v)
%vel scaled by 0.1m/s for level 9. 0.25m/s for level 1
%% movie of daily (a)geo velocities over daily mean salinity, layer 9
indices1=false([700 200]); indices1(240:10:450,10:20:end)=1;
v = VideoWriter('salinityGeoLevel9.avi');
v.FrameRate=3;
open(v)
f1=figure('renderer','zbuffer');

for i=1:148
    clf
    pcolor(XC,YC,PractSave(:,:,9,i)); shading 'flat'; caxis([36.2 37.5]); colorbar
    hold all
    %contour(XC,YC,inWag(:,:,9),[1 1],'r','linewidth',2)
    holdE=geoE(:,:,9,i);
    holdN=geoN(:,:,9,i);
    quiver([-5;XC(indices1)],[37;YC(indices1)],[0.25/0.25;holdE(indices1)./0.25],[0;holdN(indices1)./0.25],'m','linewidth',2)
    axis([-5.5 -2 35 37.5])
    title(num2str(i))
    legend('salinity','geostrophic velocities')
    %axis tight
    writeVideo(v,getframe(gcf))
    writeVideo(v,getframe(gcf))
    writeVideo(v,getframe(gcf))
    writeVideo(v,getframe(gcf))
end
close(v)
%% figures that are good examples of cross-mean geo by geo


%% figures that are good examples of cross-geo by ageo


