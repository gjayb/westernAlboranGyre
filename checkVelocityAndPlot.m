close all

addpath('/nobackup1/gbrett/mStuff')
cd('/nobackup1/gbrett/spinupsteady')

load('geometrySpinupSteady.mat');

Uave=rdmds('Uave',NaN);
[nx,ny,nz,nt]=size(Uave)
Vave=rdmds('Vave',NaN);
Wave=rdmds('Wave',NaN);

disp('Done loading')

areas=distX2.*distY2;
areas2=zeros([nx ny nz nt]);

ini=mod(nt,10);
days=ini:10:nt;

UArot=Uave.*repmat(AngleCS,[1 1 nz nt]) - Vave.*repmat(AngleSN,[1 1 nz nt]);  %U(:,:,:,i).*repmat(AngleCS,[1 1 46]) - V(:,:,:,i).*repmat(AngleSN,[1 1 46]);
VArot=Uave.*repmat(AngleSN,[1 1 nz nt]) + Vave.*repmat(AngleCS,[1 1 nz nt]); %U(:,:,:,i).*repmat(AngleSN,[1 1 46]) + V(:,:,:,i).*repmat(AngleCS,[1 1 46]);

for i=ini:10:nt
    disp(num2str(i))
    for j=1:nz
        UArot(:,:,j,i) = griddata(XU,YU,Uave(:,:,j,i),XC,YC);
        VArot(:,:,j,i) = griddata(XV,YV,Vave(:,:,j,i),XC,YC);
        areas2(:,:,j,i) = areas;
    end
end

figure; quiver(XC(1:10:end,1:10:end),YC(1:10:end,1:10:end),UArot(1:10:end,1:10:end,1,end),VArot(1:10:end,1:10:end,1,end)); hold all; plot(lonCoast,latCoast,'k')
save2pdf('UVsurfaceEnd')
close all

figure; quiver(XC(1:10:end,1:10:end),YC(1:10:end,1:10:end),UArot(1:10:end,1:10:end,17,end),VArot(1:10:end,1:10:end,17,end)); hold all; plot(lonCoast,latCoast,'k')
save2pdf('UV200mEnd')
close all

figure; quiver(XC(1:10:end,1:10:end),YC(1:10:end,1:10:end),UArot(1:10:end,1:10:end,22,end),VArot(1:10:end,1:10:end,22,end)); hold all; plot(lonCoast,latCoast,'k')
save2pdf('UV350mEnd')
close all

disp('KE calculation')

ke=areas2.*(UArot.^2 + VArot.^2);
keSeries=nansum(nansum(nansum(ke)));
keEnd=sum(ke(:,:,:,end),3);

%clear ke

figure; contourf(XC,YC,keEnd,[0 1e1 1e2 1e3 1e4 1e5 1e6 1e7 1e8 2e8 3e8]); colorbar; shading 'flat'; hold all; plot(lonCoast,latCoast,'k')
save2pdf('keEnd')
close all

if nt>10
keDiff=keEnd-sum(ke(:,:,:,end-10),3);
figure; contourf(XC,YC,keDiff); colorbar; shading 'flat'; hold all; plot(lonCoast,latCoast,'k')
save2pdf('keDiff')
close all
end

figure; plot(days,squeeze(keSeries(1,1,1,ini:10:nt)))
save2pdf('keSeries')
close all

fn=strcat('outputDay',num2str(nt));
save(fn,'-v7.3');

disp('finished')
