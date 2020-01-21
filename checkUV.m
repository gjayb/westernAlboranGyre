addpath ../mStuff
load('uvwDailyNativeGrid.mat')

[nx,ny,nz,nt]=size(u)

u1=rdmds('Uave',NaN);
[nx1,ny1,nz1,nta]=size(u1)
umaxdiff=max(max(max(max(abs(u-u1(1:nx,1:ny,1:nz,1:nt))))))
clear u

v1=rdmds('Vave',NaN);
vmaxdiff=max(max(max(max(abs(v-v1(1:nx,1:ny,1:nz,1:nt))))))

clear
