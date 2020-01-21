

addpath('/nobackup1/gbrett/mStuff')
load('geometrySpinupSteady.mat');
dBin=0.5.*dInterface(1:end-1)+0.5.*dInterface(2:end);
%load('steady148all.mat')

%load u,v, select the depth, clear the rest
%U=rdmds('Uave',NaN);
%V=rdmds('Vave',NaN);
%[nx,ny,nz,nt]=size(U)
%disp('load 1')
%ui=mean(U(:,:,:,1:148),4);
%vi=mean(V(:,:,:,1:148),4);
%clear U V

%rotate to east-north
%Urot=ui.*repmat(AngleCS,[1 1  nz ]) - vi.*repmat(AngleSN,[1 1 nz]);  
%Vrot=ui.*repmat(AngleSN,[1 1  nz ]) + vi.*repmat(AngleCS,[1 1 nz]);
%clear ui vi 
%disp('rotated')
%areas=distX2.*distY2;
%areas2=zeros([nx ny nt]);
%for i=1:nt
%    disp(num2str(i))
%        Urot(:,:,i) = griddata(XU,YU,Urot(:,:,i),XC,YC);
%        Vrot(:,:,i) = griddata(XV,YV,Vrot(:,:,i),XC,YC);
       % areas2(:,:,i) = areas;
%end
%disp('gridded')
sal1=rdmds('Save',NaN);
temp1=rdmds('Tave',NaN);
%w1=rdmds('Wave',NaN);
disp('loads done')
%Sal=mean(sal1(:,:,:,1:148),4); clear sal1
%Temp=mean(temp1(:,:,:,1:148),4); clear temp1
%W=mean(w1(:,:,:,1:148),4); clear w1
%disp('means taken')
%clear u* v* U* V* *1

S=sal1(:,:,1:16,1:148); clear sal1
T=temp1(:,:,1:16,1:148); clear temp1

fn=('varying148ts16levels.mat');%strcat('steady148all.mat');
disp('saving')
save(fn,'-v7.3');
disp('all done')
