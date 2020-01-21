

addpath('/nobackup1/gbrett/mStuff')
load('geometrySpinupSteady.mat');
%dBin=0.5.*dInterface(1:end-1)+0.5.*dInterface(2:end);

%load u,v, select the depth, clear the rest
%tEnds=[1 180:360:8640*150];
times=8640:8640:(162*8640);
U=rdmds('Uave',times);%tEnds(239:721));
%U=U(:,:,1:15,1:148);
%j=1
%ui=squeeze(U(:,:,j,:)); clear U

V=rdmds('Vave',times);%tEnds(239:721));
W=rdmds('Wave',times);
disp('save 1')
save('uvwDailyNativeNF.mat','-v7.3')
[nx,ny,nz,nt]=size(U)
j=1
ui=squeeze(U(:,:,j,:));
vi=squeeze(V(:,:,j,:));
%W=squeeze(W(:,:,j,:));
clear U V
%Wave=squeeze(mean(W(:,:,:,1:148),4));
%clear W

%rotate to east-north
Urot=ui.*repmat(AngleCS,[1 1  nt]) - vi.*repmat(AngleSN,[1 1 nt]);  
Vrot=ui.*repmat(AngleSN,[1 1  nt]) + vi.*repmat(AngleCS,[1 1 nt]); 
clear ui vi

%areas=distX2.*distY2;
%areas2=zeros([nx ny nt]);
for i=1:nt
    disp(num2str(i))
        Urot(:,:,i) = griddata(XU,YU,Urot(:,:,i),XC,YC);
        Vrot(:,:,i) = griddata(XV,YV,Vrot(:,:,i),XC,YC);
%%        areas2(:,:,i) = areas;
end

disp('save 2')
save('uvwDailyDepth1rotatedNF.mat','-v7.3')

%U=U(:,:,1:15,1:148);
%V=V(:,:,1:15,1:148);
%W=W(:,:,1:16,1:148);
%save
%fn=strcat('uvHdays10to30depth',num2str(j),'.mat');
%fn='uvwAve148.mat';
%fn='uva148levels15native.mat'
%SSHave=rdmds('SSHave',NaN);
%ssh148=SSHave(:,:,1:148);
%clear SSHave
%fn='ssh148.mat';
%save(fn,'-v7.3');
