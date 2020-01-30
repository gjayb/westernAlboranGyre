%Reads in all SSHave, Save, Tave, Uave, Vave, Wave in the current folder
%uses geometrySpinupSteady
%makes and saves plots of latest T,S,UV
%plots timeseries of kinetic energy

close all

addpath('/nobackup1/gbrett/mStuff')
%cd('/nobackup1/gbrett/spinupsteady')

load('geometrySpinupSteady.mat');
dBin=0.5.*dInterface(1:end-1)+0.5.*dInterface(2:end);
Save=rdmds('Save',NaN);
[nx,ny,nz,nt]=size(Save);

figure; contourf(XC,YC,Save(:,:,1,end),33:0.25:40);  colorbar;
hold all; plot(lonCoast,latCoast,'k')
fn1=strcat('SaveDay',num2str(nt));
save2pdf(fn1)
close all

figure; contourf(YC(350,:),dBin,squeeze(Save(350,:,:,end)).',33:0.25:40); colorbar;
set(gca,'Ydir','reverse')
ylim([0 1000])
title('WAG S')
fn2=strcat(fn1,'sec2');
save2pdf(fn2)
close all

figure; contourf(YC(250,:),dBin,squeeze(Save(250,:,:,end)).',33:0.25:40); colorbar;
set(gca,'Ydir','reverse')
ylim([0 1000])
title('Strait S')
fn2=strcat(fn1,'sec1');
save2pdf(fn2)
close all

figure; contourf(YC(500,:),dBin,squeeze(Save(500,:,:,end)).',33:0.25:40); colorbar;
set(gca,'Ydir','reverse')
title('EAG S')
ylim([0 1000])
fn2=strcat(fn1,'sec3');
save2pdf(fn2)
close all

if nt>10
figure; contourf(XC,YC,Save(:,:,1,end)-Save(:,:,1,end-10)); shading 'flat'; colorbar; 
hold all; plot(lonCoast,latCoast,'k')
fn4='Save10DayChange';
save2pdf(fn4)
close all
end

figure; contourf(XC,YC,Save(:,:,21,end),33:0.25:40); shading 'flat'; colorbar; %caxis([34 40])
hold all; plot(lonCoast,latCoast,'k')
fn9=strcat('Save300mDay',num2str(nt));
save2pdf(fn9)
close all

figure; contourf(XC,YC,Save(:,:,28,end),33:0.25:40); shading 'flat'; colorbar; %caxis([34 40])
hold all; plot(lonCoast,latCoast,'k')
fn9=strcat('Save800mDay',num2str(nt));
save2pdf(fn9)
close all

if nt>10
figure; contourf(XC,YC,Save(:,:,21,end)-Save(:,:,21,end-10)); shading 'flat'; colorbar;
hold all; plot(lonCoast,latCoast,'k')
fn7=strcat('Sdiff300mDay',num2str(nt));
save2pdf(fn7)
close all

figure; contourf(XC,YC,Save(:,:,28,end)-Save(:,:,28,end-10)); shading 'flat'; colorbar;
hold all; plot(lonCoast,latCoast,'k')
fn9=strcat('Sdiff800mDay',num2str(nt));
save2pdf(fn9)
close all
end
%
Tave=rdmds('Tave',NaN);
figure; contourf(XC,YC,Tave(:,:,1,end),3:0.25:30); shading 'flat'; colorbar; caxis([0 30])
hold all; plot(lonCoast,latCoast,'k')
fn2=strcat('TaveDay',num2str(nt));
save2pdf(fn2)
close all

figure; contourf(YC(350,:),dBin,squeeze(Tave(350,:,:,end)).',3:0.25:26); colorbar;
set(gca,'Ydir','reverse')
title('WAG S')
ylim([0 1000])
fn3=strcat(fn2,'sec2');
save2pdf(fn3)
close all

figure; contourf(YC(250,:),dBin,squeeze(Tave(250,:,:,end)).',3:0.25:26); colorbar;
set(gca,'Ydir','reverse')
title('Strait S')
ylim([0 1000])
fn3=strcat(fn2,'sec1');
save2pdf(fn3)
close all

figure; contourf(YC(500,:),dBin,squeeze(Tave(500,:,:,end)).',3:0.25:26); colorbar;
set(gca,'Ydir','reverse')
title('EAG S')
ylim([0 1000])
fn3=strcat(fn2,'sec3');
save2pdf(fn3)
close all

if nt>10
figure; contourf(XC,YC,Tave(:,:,1,end)-Tave(:,:,1,end-10)); shading 'flat'; colorbar; 
hold all; plot(lonCoast,latCoast,'k')
fn5='Tave10DayChange';
save2pdf(fn5)
close all
end

figure; contourf(XC,YC,Tave(:,:,21,end),0:1:30); shading 'flat'; colorbar; caxis([0 30])
hold all; plot(lonCoast,latCoast,'k')
fn8=strcat('Tave300mDay',num2str(nt));
save2pdf(fn8)
close all



figure; contourf(XC,YC,Tave(:,:,28,end),0:1:30); shading 'flat'; colorbar; caxis([0 30])
hold all; plot(lonCoast,latCoast,'k')
fn10=strcat('Tave800mDay',num2str(nt));
save2pdf(fn10)
close all


if nt>10
figure; contourf(XC,YC,Tave(:,:,21,end)-Tave(:,:,21,end-10)); shading 'flat'; colorbar;
hold all; plot(lonCoast,latCoast,'k')
fn8=strcat('Tdiff300mDay',num2str(nt));
save2pdf(fn8)
close all



figure; contourf(XC,YC,Tave(:,:,28,end)-Tave(:,:,28,end-10)); shading 'flat'; colorbar;
hold all; plot(lonCoast,latCoast,'k')
fn10=strcat('Tdiff800mDay',num2str(nt));
save2pdf(fn10)
close all
end
%
SSHave=rdmds('SSHave',NaN);
%disp(strcat('So far there are ',num2str(nt),' days of data.'))
figure; contourf(XC,YC,SSHave(:,:,end),[-1:0.1:-0.1 -0.075 -0.05 -0.025 -0.01 0 0.01]); shading 'flat'; colorbar;
hold all; plot(lonCoast,latCoast,'k')
fn3=strcat('SSHaveDay',num2str(nt));
save2pdf(fn3)
close all

if nt>10
figure; contourf(XC,YC,SSHave(:,:,end)-SSHave(:,:,end-10)); shading 'flat'; colorbar;
hold all; plot(lonCoast,latCoast,'k')
fn6='SSHave10DayChange';
save2pdf(fn6)
close all
end


%
Uave=rdmds('Uave',NaN);
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
        UArot(:,:,j,i) = griddata(XU,YU,UArot(:,:,j,i),XC,YC);
        VArot(:,:,j,i) = griddata(XV,YV,VArot(:,:,j,i),XC,YC);
        areas2(:,:,j,i) = areas;
    end
end

figure; quiver(XC(1:10:end,1:10:end),YC(1:10:end,1:10:end),UArot(1:10:end,1:10:end,1,end),VArot(1:10:end,1:10:end,1,end)); hold all; plot(lonCoast,latCoast,'k')
fn3=strcat('UVsurfaceDay',num2str(nt));
save2pdf(fn3)
close all

figure; quiver(XC(1:10:end,1:10:end),YC(1:10:end,1:10:end),UArot(1:10:end,1:10:end,17,end),VArot(1:10:end,1:10:end,17,end)); hold all; plot(lonCoast,latCoast,'k')
fn3=strcat('UV200mDay',num2str(nt));
save2pdf(fn3)
close all

figure; quiver(XC(1:10:end,1:10:end),YC(1:10:end,1:10:end),UArot(1:10:end,1:10:end,22,end),VArot(1:10:end,1:10:end,22,end)); hold all; plot(lonCoast,latCoast,'k')
fn3=strcat('UV350mDay',num2str(nt));
save2pdf(fn3)
close all

disp('KE calculation')

ke=areas2.*(UArot.^2 + VArot.^2);
keSeries=nansum(nansum(nansum(ke)));
keEnd=sum(ke(:,:,:,end),3);
if nt>10
    keDiff=keEnd-sum(ke(:,:,:,end-10),3);
end
clear ke

figure; contourf(XC,YC,keEnd,[0 1e1 1e2 1e3 1e4 1e5 1e6 1e7 1e8 2e8 3e8]); colorbar; shading 'flat'; hold all; plot(lonCoast,latCoast,'k')
save2pdf('keEnd')
close all

if nt>10
figure; contourf(XC,YC,keDiff); colorbar; shading 'flat'; hold all; plot(lonCoast,latCoast,'k')
save2pdf('keDiff')
close all
end

figure; plot(days,squeeze(keSeries(1,1,1,ini:10:nt)))
save2pdf('keSeries')
close all

% fn=strcat('outputDay',num2str(nt));
% save(fn,'-v7.3');

disp('finished')
