%reads in all S and T, plots last 2 snapshots

close all

addpath('/nobackup1/gbrett/mStuff')
cd('/nobackup1/gbrett/spinupsteady')

load('geometrySpinupSteady.mat');
S=rdmds('S',NaN);
[nx,ny,nz,nt]=size(S);

figure; contourf(XC,YC,S(:,:,1,end),30:0.5:40); shading 'flat'; colorbar; %caxis([34 40])
hold all; plot(lonCoast,latCoast,'k')
fn1=strcat('SHour',num2str(nt));
save2pdf(fn1)
close all

figure; contourf(XC,YC,S(:,:,1,end-1),30:0.5:40); shading 'flat'; colorbar; %caxis([34 40])
hold all; plot(lonCoast,latCoast,'k')
fn1=strcat('SHour',num2str(nt-1));
save2pdf(fn1)
close all

figure; contourf(XC,YC,S(:,:,21,end),30:0.5:40); shading 'flat'; colorbar; %caxis([34 40])
hold all; plot(lonCoast,latCoast,'k')
fn7=strcat('S300mHour',num2str(nt));
save2pdf(fn7)
close all

figure; contourf(XC,YC,S(:,:,21,end-1),30:0.5:40); shading 'flat'; colorbar; %caxis([34 40])
hold all; plot(lonCoast,latCoast,'k')
fn7=strcat('S300mHour',num2str(nt-1));
save2pdf(fn7)
close all

clear S

load('geometrySpinupSteady.mat');
T=rdmds('T',NaN);

figure; contourf(XC,YC,T(:,:,1,end),30); shading 'flat'; colorbar; %caxis([34 40])
hold all; plot(lonCoast,latCoast,'k')
fn1=strcat('THour',num2str(nt));
save2pdf(fn1)
close all

figure; contourf(XC,YC,T(:,:,1,end-1),30); shading 'flat'; colorbar; %caxis([34 40])
hold all; plot(lonCoast,latCoast,'k')
fn1=strcat('THour',num2str(nt-1));
save2pdf(fn1)
close all

figure; contourf(XC,YC,T(:,:,21,end),30); shading 'flat'; colorbar; %caxis([34 40])
hold all; plot(lonCoast,latCoast,'k')
fn7=strcat('T300mHour',num2str(nt));
save2pdf(fn7)
close all

figure; contourf(XC,YC,T(:,:,21,end-1),30); shading 'flat'; colorbar; %caxis([34 40])
hold all; plot(lonCoast,latCoast,'k')
fn7=strcat('T300mHour',num2str(nt-1));
save2pdf(fn7)
close all

clear T