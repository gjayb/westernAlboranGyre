%see what happens to nitrate levels using NP model and mean horizontal
%surface circulation; NP along 28day trajectories, binned, movie

%trajectories
load('uvwAve148.mat','*rot');
load('uvwAve148.mat','*C');
u=Urot(:,:,1);
v=Vrot(:,:,1);
load('nitrateMEDAR.mat')
load('coastAnd20110915.mat', 'latCoast')
load('coastAnd20110915.mat', 'lonCoast')
lonIn=long;%(45:77,51:111);
latIn=latg;%(45:77,51:111);

xmin=min(min(XC));
ymin=min(min(YC));
xcm=(XC-xmin).*111000.*cosd(YC);
ycm=(YC-ymin).*111000;
x1=(double(lon)-xmin).*111000.*cosd(36);
y1=(double(lat)-ymin).*111000;
[xg,yg]=meshgrid(x1,y1);
xin=(lonIn-xmin).*111000.*cosd(latIn);
yin=(latIn-ymin).*111000;
z0=[xin(:)';yin(:)'];

%%
u0=griddata(xcm,ycm,u,xg,yg);
v0=griddata(xcm,ycm,v,xg,yg);

timesWanted=0:1:28;
options=odeset('RelTol',10^(-5),'AbsTol',10^(-8));
disp('entering integration')
[ttr,positions]=ode45(@HamEqSolver_BiLin_Irina,timesWanted*86400,z0(:),options,repmat(u0,[1 1 2]),repmat(v0,[1 1 2]),x1,y1,86400*[-1 timesWanted(end)+1]);

xtr=positions(:,1:2:end-1);
ytr=positions(:,2:2:end);
disp('integration done')
lattr1=ones(size(ytr)).*ymin+ytr./111000;
lontr1=ones(size(xtr)).*xmin+xtr./(111000.*cosd(lattr1));
%% biology, NP, Irina model
mu1=1; eps1=0.1; sigma1=0.1; %per day; mu 0.4 from Irina, 1 from Dennis
n0=0.5;%mmol nitrate/m^3
p0=0.01;%mmol Chl/m^3
nIn=nitrateW(21,:,:);

npin=[nIn(:).';p0*ones(size(nIn(:))).'];
npin=npin(:);

disp('entering integration')
options=odeset('RelTol',10^(-13),'AbsTol',10^(-17));
[ttr,np]=ode45(@dnpdt,timesWanted,npin,options,mu1,n0,eps1,sigma1);
disp('integration done')
  Ntr=np(:,1:2:end-1); 
  Ptr=np(:,2:2:end);
Ptr(Ptr<0)=0;
Ntr(Ntr<0)=0;
  %% binning
  
  xbin=xg(1:3:end,1:3:end);
  ybin=yg(1:3:end,1:3:end);
  [ny,nx]=size(xbin);
  
  
ntraj=zeros([ny-1 nx-1 29]);
binN=ntraj;
binP=ntraj;
for k=1:29
    k
  for i=1:ny-1
      for j=1:nx-1
          in1=find(xtr(k,:)>=xbin(i,j) & xbin(i,j+1)>xtr(k,:) & ytr(k,:)>=ybin(i,j) & ybin(i+1,j)>ytr(k,:));
          ntraj(i,j,k)=length(in1);
          if length(in1)>0
              binN(i,j,k)=nanmean(Ntr(k,in1));
              binP(i,j,k)=nanmean(Ptr(k,in1));
          end
      end
  end
end
  
%% plot
  daysi=[1 8 15 22; 2 9 16 23;...
  3 10 17 24;4 11 18 25;5 12 19 26;6 13 20 27; 7 14 21 28];
meanN=zeros([ny-1 nx-1 7]);
meanP=meanN;
meanNT=meanN;
for k=1:7
  meanN(:,:,k)=nansum(binN(:,:,daysi(k,:)).*ntraj(:,:,daysi(k,:)),3)./nansum(ntraj(:,:,daysi(k,:)),3);
  meanP(:,:,k)=nansum(binP(:,:,daysi(k,:)).*ntraj(:,:,daysi(k,:)),3)./nansum(ntraj(:,:,daysi(k,:)),3);
  meanNT(:,:,k)=nansum(ntraj(:,:,daysi(k,:)),3);
end
  
  
  figure; pcolor(xbin(1:end-1,1:end-1),ybin(1:end-1,1:end-1),mean(ntraj,3)); shading 'flat'; colorbar
  title('mean number of trajectories per bin')
  
  for k=1:7
  figure; subplot(1,2,1); pcolor(xbin(1:end-1,1:end-1),ybin(1:end-1,1:end-1),meanN(:,:,k)); shading 'flat'; colorbar
  caxis([0 1])
  title('mean N')
subplot(1,2,2); pcolor(xbin(1:end-1,1:end-1),ybin(1:end-1,1:end-1),meanP(:,:,k)); shading 'flat'; colorbar
  title('mean P')
  caxis([0 1])
  end
  
%     for k=1:4:29
%   figure; subplot(1,2,1); pcolor(xbin(1:end-1,1:end-1),ybin(1:end-1,1:end-1),binN(:,:,k)); shading 'flat'; colorbar
%   caxis([0 1])
%   title('bin N')
% subplot(1,2,2); pcolor(xbin(1:end-1,1:end-1),ybin(1:end-1,1:end-1),binP(:,:,k)); shading 'flat'; colorbar
%   title('bin P')
%   caxis([0 1])
%     end
  
  %% plots releasing every 4 days
  
    daysi2=[1 5 9 13 17 21 25; 2 6 10 14 18 22 26;...
            3 7 11 15 19 23 27; 4 8 12 16 20 24 28];
meanN2=zeros([ny-1 nx-1 4]);
meanP2=meanN2;
meanNT2=meanN2;
for k=1:4
  meanN2(:,:,k)=nansum(binN(:,:,daysi2(k,:)).*ntraj(:,:,daysi2(k,:)),3)./nansum(ntraj(:,:,daysi2(k,:)),3);
  meanP2(:,:,k)=nansum(binP(:,:,daysi2(k,:)).*ntraj(:,:,daysi2(k,:)),3)./nansum(ntraj(:,:,daysi2(k,:)),3);
  meanNT2(:,:,k)=nansum(ntraj(:,:,daysi2(k,:)),3);
end
  
  
  for k=1:4
  figure; subplot(1,2,1); pcolor(xbin(1:end-1,1:end-1),ybin(1:end-1,1:end-1),meanN2(:,:,k)); shading 'flat'; colorbar
  caxis([0 1])
  title('mean N2')
subplot(1,2,2); pcolor(xbin(1:end-1,1:end-1),ybin(1:end-1,1:end-1),meanP2(:,:,k)); shading 'flat'; colorbar
  title('mean P2')
  caxis([0 1])
  end
 %% every day version
 meanN3=nansum(binN.*ntraj,3)./nansum(ntraj,3);
  meanP3=nansum(binP.*ntraj,3)./nansum(ntraj,3);
  meanNT3=nansum(ntraj,3);
  
  figure; pcolor(xbin(1:end-1,1:end-1),ybin(1:end-1,1:end-1),meanN3); shading 'flat'; colorbar
  caxis([0 1])
  title('mean N3')
figure; pcolor(xbin(1:end-1,1:end-1),ybin(1:end-1,1:end-1),meanP3); shading 'flat'; colorbar
  title('mean P3')
  caxis([0 1])
  
  figure; pcolor(xbin(1:end-1,1:end-1),ybin(1:end-1,1:end-1),meanNT3); shading 'flat'; colorbar
  title('points per bin 3')
  caxis([1 20])
  
  %%
Ntr0=repmat(nIn(:).',[29 1]);
Ptr0=zeros(size(Ntr0));
binN0=zeros(size(ntraj));
for k=1:29
    k
  for i=1:ny-1
      for j=1:nx-1
          in1=find(xtr(k,:)>=xbin(i,j) & xbin(i,j+1)>xtr(k,:) & ytr(k,:)>=ybin(i,j) & ybin(i+1,j)>ytr(k,:));
          if length(in1)>0
              binN0(i,j,k)=nanmean(Ntr0(k,in1));
          end
      end
  end
end

meanN03=nansum(binN0.*ntraj,3)./nansum(ntraj,3);

figure; pcolor(xbin(1:end-1,1:end-1),ybin(1:end-1,1:end-1),meanN03); shading 'flat'; colorbar
caxis([0.1 3]); title('N, advection-only 28 daily release mean circulation')
figure; pcolor(xbin(1:end-1,1:end-1),ybin(1:end-1,1:end-1),binN0(:,:,1)); shading 'flat'; colorbar
caxis([0.1 3]); title('Initial N, binned')
figure; pcolor(xbin(1:end-1,1:end-1),ybin(1:end-1,1:end-1),meanN03-binN0(:,:,1)); shading 'flat'; colorbar
title('change in N, advection only')
%%
figure; contourf(xbin(1:end-1,1:end-1),ybin(1:end-1,1:end-1),meanN3-binN0(:,:,1),-1:0.1:1);  colorbar
  caxis([-0.9 0.9])
  title('change in N, advection+biology')
%% decide if steady-state
cumN03=zeros(size(binN0));
cumN03(:,:,1)=binN0(:,:,1);
for k=2:29
    cumN03(:,:,k)=nansum(binN0(:,:,1:k).*ntraj(:,:,1:k),3)./nansum(ntraj(:,:,1:k),3);
end

for k=[2 8 15 22 29]
    figure; pcolor(xbin(1:end-1,1:end-1),ybin(1:end-1,1:end-1),cumN03(:,:,k)); shading 'flat'; colorbar
end