

addpath('/nobackup1/gbrett/mStuff')
load('geometrySpinupSteady.mat');
dBin=0.5.*dInterface(1:end-1)+0.5.*dInterface(2:end);

%load u,v, select the depth, clear the rest
U=rdmds('Uave',NaN);
V=rdmds('Vave',NaN);
[nx,ny,nz,nt]=size(U)
j=1
ui=U(:,:,j,:);
vi=V(:,:,j,:);
clear U V

%rotate to east-north
Urot=ui.*repmat(AngleCS,[1 1 1 nt]) - vi.*repmat(AngleSN,[1 1 1 nt]);  
Vrot=ui.*repmat(AngleSN,[1 1 1 nt]) + vi.*repmat(AngleCS,[1 1 1 nt]); 

areas=distX2.*distY2;
areas2=zeros([nx ny nz nt]);
for i=1:nt
    disp(num2str(i))
        Urot(:,:,j,i) = griddata(XU,YU,Urot(:,:,j,i),XC,YC);
        Vrot(:,:,j,i) = griddata(XV,YV,Vrot(:,:,j,i),XC,YC);
        areas2(:,:,j,i) = areas;
end

%save
fn=strcat('uvAdepth',num2str(i),'.mat');
save(fn,'-v7.3');