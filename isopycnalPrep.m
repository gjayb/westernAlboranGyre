%loads potential temp and salinity daily records, calculates density for a
%subsection, interpolates the density to a fine grid, saves

addpath('/nobackup1/gbrett/mStuff')
load('geometrySpinupSteady.mat');

dBin=0.5.*dInterface(1:end-1)+0.5.*dInterface(2:end);

times=8640:8640:8640*176;%175 days, so could run 27 day manifold for 148th day
xDeg=XC(230:500,1:130);
yDeg=YC(230:500,1:130);
xmin=min(min(xDeg));
ymin=min(min(yDeg));
xM=111000*cosd(yDeg).*(xDeg-xmin*ones(size(xDeg)));
yM=111000*(yDeg-ymin.*ones(size(ymin)));

% %% load
% PT=zeros([271 130 37 length(times)]);
% S=PT;
% Rho=PT;
% 
% for i=1:length(times)
%     i
%     t1=rdmds('Tave',times(i));
%     s1=rdmds('Save',times(i));
%     PT(:,:,:,i)=t1(230:500,1:130,1:37);
%     S(:,:,:,i)=s1(230:500,1:130,1:37);
% end
% %% calculate ct
% CT=gsw_CT_from_pt(S,PT); clear PT
% %% rho
% 
% pbin=dBin./10;
% 
% for j=1:37
%     j
%     Rho(:,:,j,:)=gsw_rho(S(:,:,j,:),CT(:,:,j,:),pbin(j));
% end
%     clear S CT
% %%  save
% disp('saving rhoVary175.mat')
% save('rhoVary175.mat','-v7.3')

%% velocities
%clear Rho

u=zeros([700 200 46 length(times)]);
v=u;

for i=1:length(times)
    i
    U1=rdmds('Uave',times(i));
    V1=rdmds('Vave',times(i));
    for j=1:46
        U(:,:,j) = griddata(XU,YU,U1(:,:,j),XC,YC);
        V(:,:,j) = griddata(XV,YV,V1(:,:,j),XC,YC);
    end
    u(:,:,:,i)=U.*repmat(AngleCS,[1 1 46]) - V.*repmat(AngleSN,[1 1 46]);  
    v(:,:,:,i)=U.*repmat(AngleSN,[1 1 46]) + V.*repmat(AngleCS,[1 1 46]); 

end
clear U V U1 V1

w=rdmds('Wave',times);
%% better grid

xmin=min(min(XC));
ymin=min(min(YC));
xinM=111000*cosd(YC).*(XC-xmin*ones(size(XC)));
yinM=111000*(YC-ymin.*ones(size(YC)));
xvel=0:1000:max(max(xinM));
yvel=0:1000:max(max(yinM));

[xvelg,yvelg]=meshgrid(xvel,yvel);
[ng,mg]=size(xvelg);
U=zeros([ng mg 46 length(times)]);
V=U; W=U;
for k=1:length(times)
    k
    for j=1:46
        U(:,:,j,k)=griddata(xinM,yinM,u(:,:,j,k),xvelg,yvelg);
        V(:,:,j,k)=griddata(xinM,yinM,v(:,:,j,k),xvelg,yvelg);
        W(:,:,j,k)=griddata(xinM,yinM,w(:,:,j,k),xvelg,yvelg);
    end
end
clear u v w
%% save
disp('saving')
save('uvwDaily175.mat','-v7.3')
disp('done')