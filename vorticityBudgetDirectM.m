%% load for circulation
load('edgesWAGeuler2017NF.mat', 'inWag','XC','YC')
%load('uvwSSHDailyDepth1rotated148F.mat','Urot','Vrot','SSHa')

load('distancesAreas')
load('geometrySpinupSteady')
xmin=min(XC(:)); ymin=min(YC(:));
xm=111111*cosd(YC).*(XC-xmin*ones(size(XC)));
ym=111111*(YC-ymin*ones(size(YC)));
xg=111111*cosd(YG).*(XG-xmin*ones(size(XG)));
yg=111111*(YG-ymin*ones(size(YG)));
xcoast=111111*cosd(latCoast).*(lonCoast-xmin*ones(size(lonCoast)));
ycoast=111111*(latCoast-ymin*ones(size(latCoast)));
xum=111111*cosd(YU).*(XU-xmin*ones(size(XU)));
yum=111111*(YU-ymin*ones(size(YU)));
xvm=111111*cosd(YV).*(XV-xmin*ones(size(XV)));
yvm=111111*(YV-ymin*ones(size(YV)));
dZ(1,1,:)=diff(dInterface);
DXC=reshape(dxc,[700 200]); DYC=reshape(dyc,[700 200]);
rAw=reshape(raw,[700 200]); rAs=reshape(ras,[700 200]);
cellVolU=repmat(rAw,[1 1 46]).*hFacW.*repmat(dZ,[700 200 1]);
cellVolV=repmat(rAs,[1 1 46]).*hFacS.*repmat(dZ,[700 200 1]);
%% update inWag with new layers
%inWag3=inWag(:,:,1).*double(d>200);
%figure; [c,h]=contour(xm,ym,inWag3(:,:,1),[1 1]);
%edgeWagX=c(1,2:444);
%edgeWagY=c(2,2:444);
%inWag5=inWag(:,:,1).*double(d>50);
%figure; [c,h]=contour(xm,ym,inWag5(:,:,1),[1 1]);
%edgeWagX=c(1,2:489);
%edgeWagY=c(2,2:489);
inWag2=inWag.*hFacW.*hFacS;
%figure(1); hold all;
%for i=1:16
%%inWag2(:,:,i)=inWag(:,:,i).*hFacW(:,:,i).*hFacS(:,:,i);
%[c,h]=contour(xm,ym,inWag2(:,:,i),[1 1]);
%np=c(2,1);
%edgeWagX(i,1:np)=c(1,2:np+1);
%edgeWagY(i,1:np)=c(2,2:np+1);
%inside(:,:,i)=inpolygon(xg,yg,edgeWagX(i,1:np),edgeWagY(i,1:np));%&(hFacW(:,:,i)==1)&(hFacS(:,:,i)==1);
%%hole(:,:,i)=inpolygon(xg,yg,edgeWagX(i,1:np),edgeWagY(i,1:np))&(~inside(:,:,i));
%end
%close(1)
inside=logical(inWag2);
insideSize=size(inside)
%%

%% square edge
% xi1=375;%300
% xi2=425;%350
% yi1=100;%60
% yi2=150;%100
% x1=xm(xi1,yi1:yi2);%left bottom to top, !corners!
% y1=ym(xi1,yi1:yi2);
% x2=xm(xi1:xi2,yi2);%top left to right
% y2=ym(xi1:xi2,yi2);
% x3=xm(xi2,yi1:yi2);%right bottom to top
% y3=ym(xi2,yi1:yi2);
% x4=xm(xi1:xi2,yi1);%bottom left to right
% y4=ym(xi1:xi2,yi1);
% x2=x2(:); y2=y2(:); x1=x1(:); y1=y1(:); 
% x4=x4(:); y4=y4(:); x3=x3(:); y3=y3(:); 
% xSquare=[x4(1:end-1); x3(1:end-1);x2(end:-1:2);x1(end:-1:1)];
% ySquare=[y4(1:end-1); y3(1:end-1);y2(end:-1:2);y1(end:-1:1)];
% ulogic2=false([700 200]); ulogic4=ulogic2; ulogic4(xi1+1:xi2,yi1)=true; ulogic2(xi1+1:xi2,yi2)=true;
% vlogic1=false([700 200]); vlogic3=vlogic1; vlogic1(xi1,yi1+1:yi2)=true; vlogic3(xi2,yi1+1:yi2)=true;
% load('distancesAreas')
% DYU=reshape(dyu,[700 200]);
% DXV=reshape(dxv,[700 200]);
% DYC=reshape(dyc,[700 200]);
% DXC=reshape(dxc,[700 200]);
% %earthR=6371.0*1e3;
% ds1=DYC(xi1,yi1+1:yi2);%sqrt(diff([x4(1); x1(:)]).^2 +diff([y4(1); y1(:)]).^2);%
% ds2=DXC(xi1+1:xi2,yi2);%sqrt(diff([x1(end); x2(:)]).^2 +diff([y1(end);y2(:)]).^2);%
% ds3=DYC(xi2,yi1+1:yi2);%sqrt(diff([x3(:); x2(end)]).^2 +diff([y3(:); y2(end)]).^2);%
% ds4=DXC(xi1+1:xi2,yi1);%sqrt(diff([x4(:); x3(1)]).^2 +diff([y4(:); y3(1)]).^2);%
% ds1=ds1(:);ds2=ds2(:);ds3=ds3(:);ds4=ds4(:);
% %[~,~,nt]=size(Uext);
% dx=diff([xSquare(:); xSquare(1)]); dy=diff([ySquare(:); ySquare(1)]);
% vortLogic=false([700 200]); vortLogic(xi1+1:xi2,yi1+1:yi2)=true;
%% any edge
%inside=false([700 200]); %inside(xi1+1:xi2,yi1+1:yi2)=true;%this works! %inside=logical(inWag(:,:,1));%this doesn't work- probably coast
openN=false([700 200 16]); openE=openN; openS=openN; openW=openN;
for i=1:16
openN(:,1:end-1,i)=inside(:,1:end-1,i)&~inside(:,2:end,i);
openS(:,1:end-1,i)=inside(:,2:end,i)&~inside(:,1:end-1,i);
openE(1:end-1,:,i)=inside(1:end-1,:,i)&~inside(2:end,:,i);
openW(1:end-1,:,i)=inside(2:end,:,i)&~inside(1:end-1,:,i);
end
sizeON=size(openN)%
islog=islogical(openN)
%%
signs=ones([700 200 16]);
% signs(hole)=-1;
% for i=7:14
%    signs(:,1:end-1,i)=min(signs(:,1:end-1,i),signs(:,2:end,i));
%    signs(1:end-1,:,i)=min(signs(1:end-1,:,i),signs(2:end,:,i)); 
%    signs(:,2:end,i)=min(signs(:,2:end,i),signs(:,1:end-1,i));
%    signs(2:end,:,i)=min(signs(2:end,:,i),signs(1:end-1,:,i));
% end

%% load momentum terms
%% diffusion prep
DYU=reshape(dyu,[700 200]);
DXV=reshape(dxv,[700 200]);
DYC=reshape(dyc,[700 200]);
DXC=reshape(dxc,[700 200]);
rAz=reshape(raz,[700 200]);
vortLogic=inside;
dZ(1,1,:)=diff(dInterface);
rAw=reshape(raw,[700 200]);
rAs=reshape(ras,[700 200]);
cellVolU=repmat(rAw,[1 1 46]).*hFacW.*repmat(dZ,[700 200 1]);
cellVolV=repmat(rAs,[1 1 46]).*hFacS.*repmat(dZ,[700 200 1]);
load('momentumDiagnostics148dayNF2.mat','VisZ*')
UDif2a=(VisZU(:,:,2:end,:)-VisZU(:,:,1:end-1,:))./cellVolU(:,:,1:end-1,:);
VDif2a=(VisZV(:,:,2:end,:)-VisZV(:,:,1:end-1,:))./cellVolV(:,:,1:end-1,:);
UDif2a(isnan(UDif2a))=0;
VDif2a(isnan(VDif2a))=0;
clear VisZ* cellVol*
load('momentumDiagnostics148dayNF1.mat','*tend')

%% vorticity budget components, by layer

for k=1:16
    k
    hfw=hFacW(:,:,k);
    hfs=hFacS(:,:,k);
    hfc=hFacC(:,:,k);
    signs1=signs(:,:,k);
for i=1:148%nt
     hold1=Utend(:,:,k,i)/86400;
    hold2=Vtend(:,:,k,i)/86400;
    tendency(k,i)=sum(-hold1(openN(:,:,k)).*DXC(openN(:,:,k)).*hfw(openN(:,:,k)))+sum(hold1(openS(:,:,k)).*DXC(openS(:,:,k)).*hfw(openS(:,:,k)))...
        +sum(-hold2(openW(:,:,k)).*DYC(openW(:,:,k)).*hfs(openW(:,:,k)))+sum(hold2(openE(:,:,k)).*DYC(openE(:,:,k)).*hfs(openE(:,:,k))); 
hold1=UDif2a(:,:,k,i);%VisZU(:,:,k,i)./cellVolU(:,:,k); 
    hold2=VDif2a(:,:,k,i);%VisZV(:,:,k,i)./cellVolV(:,:,k);
      hold1(isnan(hold1))=0; hold2(isnan(hold2))=0;
      difftop(k,i)=sum(-hold1(openN(:,:,k)).*DXC(openN(:,:,k)).*hfw(openN(:,:,k)))+sum(hold1(openS(:,:,k)).*DXC(openS(:,:,k)).*hfw(openS(:,:,k)))...
        +sum(-hold2(openW(:,:,k)).*DYC(openW(:,:,k)).*hfs(openW(:,:,k)))+sum(hold2(openE(:,:,k)).*DYC(openE(:,:,k)).*hfs(openE(:,:,k)));
 hold1(isnan(hold1))=0; hold2(isnan(hold2))=0;
      diffbot(k,i)=sum(-hold1(openN(:,:,k)).*DXC(openN(:,:,k)).*hfw(openN(:,:,k)))+sum(hold1(openS(:,:,k)).*DXC(openS(:,:,k)).*hfw(openS(:,:,k)))...
        +sum(-hold2(openW(:,:,k)).*DYC(openW(:,:,k)).*hfs(openW(:,:,k)))+sum(hold2(openE(:,:,k)).*DYC(openE(:,:,k)).*hfs(openE(:,:,k)));
end
end
clear Utend Vtend UDif2a VDif2a

load('momentumDiagnostics148dayNF1.mat','*Diss','Adv*')

for k=1:16
    k
    hfw=hFacW(:,:,k);
    hfs=hFacS(:,:,k);
    hfc=hFacC(:,:,k);
    signs1=signs(:,:,k);
    
   for i=1:148%nt
    
    hold1=AdvU(:,:,k,i);
    hold2=AdvV(:,:,k,i);
    advection(k,i)=sum(-hold1(openN(:,:,k)).*DXC(openN(:,:,k)).*hfw(openN(:,:,k)))+sum(hold1(openS(:,:,k)).*DXC(openS(:,:,k)).*hfw(openS(:,:,k)))...
        +sum(-hold2(openW(:,:,k)).*DYC(openW(:,:,k)).*hfs(openW(:,:,k)))+sum(hold2(openE(:,:,k)).*DYC(openE(:,:,k)).*hfs(openE(:,:,k)));
    hold1=UDiss(:,:,k,i);
    hold2=VDiss(:,:,k,i);
    dissipation(k,i)=sum(-hold1(openN(:,:,k)).*DXC(openN(:,:,k)).*hfw(openN(:,:,k)))+sum(hold1(openS(:,:,k)).*DXC(openS(:,:,k)).*hfw(openS(:,:,k)))...
        +sum(-hold2(openW(:,:,k)).*DYC(openW(:,:,k)).*hfs(openW(:,:,k)))+sum(hold2(openE(:,:,k)).*DYC(openE(:,:,k)).*hfs(openE(:,:,k)));%
   end
end
clear AdvU AdvV UDiss VDiss

load('momentumDiagnostics148dayNF2.mat')
clear VisZU VisZV
load('momentumCori.mat')
for k=1:16
    k
    hfw=hFacW(:,:,k);
    hfs=hFacS(:,:,k);
    hfc=hFacC(:,:,k);
    signs1=signs(:,:,k);
    
   for i=1:148%nt
              hold1=Uext(:,:,k,i); hold2=Vext(:,:,k,i);
              windstress(k,i)=sum(-hold1(openN(:,:,k)).*DXC(openN(:,:,k)).*hfw(openN(:,:,k)))+sum(hold1(openS(:,:,k)).*DXC(openS(:,:,k)).*hfw(openS(:,:,k)))...
        +sum(-hold2(openW(:,:,k)).*DYC(openW(:,:,k)).*hfs(openW(:,:,k)))+sum(hold2(openE(:,:,k)).*DYC(openE(:,:,k)).*hfs(openE(:,:,k))); 
        hold1=VisZU(:,:,k+1,i)./cellVolU(:,:,k); hold2=VisZV(:,:,k+1,i)./cellVolV(:,:,k);
         hold1=UdPdx(:,:,k,i);
    hold2=VdPdy(:,:,k,i);
    pressClin(k,i)=sum(-hold1(openN(:,:,k)).*DXC(openN(:,:,k)).*hfw(openN(:,:,k)))+sum(hold1(openS(:,:,k)).*DXC(openS(:,:,k)).*hfw(openS(:,:,k)))...
        +sum(-hold2(openW(:,:,k)).*DYC(openW(:,:,k)).*hfs(openW(:,:,k)))+sum(hold2(openE(:,:,k)).*DYC(openE(:,:,k)).*hfs(openE(:,:,k)));
hold1=UCori(:,:,k,i); hold2=VCori(:,:,k,i);
coriolis(k,i)=sum(-hold1(openN(:,:,k)).*DXC(openN(:,:,k)).*hfw(openN(:,:,k)))+sum(hold1(openS(:,:,k)).*DXC(openS(:,:,k)).*hfw(openS(:,:,k)))...
        +sum(-hold2(openW(:,:,k)).*DYC(openW(:,:,k)).*hfs(openW(:,:,k)))+sum(hold2(openE(:,:,k)).*DYC(openE(:,:,k)).*hfs(openE(:,:,k)));

   end
end
clear Uext Vext UdPdx VdPdy UCori VCori
load('momentumDiagnostics148dayNF3.mat')
clear SSH
for k=1:16
    k
    hfw=hFacW(:,:,k);
    hfs=hFacS(:,:,k);
    hfc=hFacC(:,:,k);
    signs1=signs(:,:,k);
    
   for i=1:148%nt
           hold1=AbU(:,:,k,i);
    hold2=AbV(:,:,k,i);
    timestep(k,i)=sum(-hold1(openN(:,:,k)).*DXC(openN(:,:,k)).*hfw(openN(:,:,k)))+sum(hold1(openS(:,:,k)).*DXC(openS(:,:,k)).*hfw(openS(:,:,k)))...
        +sum(-hold2(openW(:,:,k)).*DYC(openW(:,:,k)).*hfs(openW(:,:,k)))+sum(hold2(openE(:,:,k)).*DYC(openE(:,:,k)).*hfs(openE(:,:,k)));
 hold1=Usd(:,:,k,i);
    hold2=Vsd(:,:,k,i);
    sidedrag(k,i)=sum(-hold1(openN(:,:,k)).*DXC(openN(:,:,k)).*hfw(openN(:,:,k)))+sum(hold1(openS(:,:,k)).*DXC(openS(:,:,k)).*hfw(openS(:,:,k)))...
        +sum(-hold2(openW(:,:,k)).*DYC(openW(:,:,k)).*hfs(openW(:,:,k)))+sum(hold2(openE(:,:,k)).*DYC(openE(:,:,k)).*hfs(openE(:,:,k)));
hold1=Ubd(:,:,k,i);
    hold2=Vbd(:,:,k,i);
    botdrag(k,i)=sum(-hold1(openN(:,:,k)).*DXC(openN(:,:,k)).*hfw(openN(:,:,k)))+sum(hold1(openS(:,:,k)).*DXC(openS(:,:,k)).*hfw(openS(:,:,k)))...
        +sum(-hold2(openW(:,:,k)).*DYC(openW(:,:,k)).*hfs(openW(:,:,k)))+sum(hold2(openE(:,:,k)).*DYC(openE(:,:,k)).*hfs(openE(:,:,k)));
end

end
disp('saving')
save('vorticityBudget2D148.mat')
%% difftop(1,:)=-windstress;
%%
runthis=false
if runthis
for i=[2 9]%6:15
% figure; plot(-tendency(i,:),'k','linewidth',2); hold all; plot(-difftop(i,:),'linewidth',2);plot(diffbot(i,:),'linewidth',2);
% plot(advection(i,:),'linewidth',2);plot(dissipation(i,:),'linewidth',2);plot(pressClin(i,:),'linewidth',2);plot(timestep(i,:),'linewidth',2);
% plot(-tendency(i,:)-difftop(i,:)+dissipation(i,:)+advection(i,:)+timestep(i,:)+pressClin(i,:)+diffbot(i,:),'g--','linewidth',2)
% legend('-d/dt','top stress diff','bot diffusion','advection','dissipation','baroclinic pressure','timestep','total')
% title(num2str(i))

% figure; plot(-tendency(i,:)-difftop(i,:)+diffbot(i,:)+advection(i,:)+dissipation(i,:),'linewidth',2);
% hold all; plot(pressClin(i,:),'linewidth',2);plot(timestep(i,:),'linewidth',2);
% plot(-tendency(i,:)-difftop(i,:)+dissipation(i,:)+advection(i,:)+timestep(i,:)+pressClin(i,:)+diffbot(i,:),'g--','linewidth',2)
% legend('-d/dt+diff+adv+drag','baroclinic pressure','timestep','total')
% title(num2str(i))
% 
figure; plot(-tendency(i,:),'k','linewidth',2); hold all; plot(-difftop(i,:),'linewidth',2);plot(diffbot(i,:),'linewidth',2);
plot(advection(i,:),'linewidth',2);plot(dissipation(i,:),'linewidth',2);
plot(-tendency(i,:)-difftop(i,:)+dissipation(i,:)+advection(i,:)+diffbot(i,:),'g--','linewidth',2)
legend('-d/dt','top stress diff','bot diffusion','advection','dissipation','total')
title(num2str(i))
end

% figure; plot(-tendency,'k','linewidth',2); hold all; plot(wind,'linewidth',2);plot(diffusion,'linewidth',2);plot(dissipation,'linewidth',2);
% plot(advection,'linewidth',2);plot(pressure,'linewidth',2);plot(timestep,'linewidth',2);
% plot(-tendency+wind+dissipation+advection+timestep+pressure+diffusion,'g--','linewidth',2)
% legend('-d/dt','wind','vertical diffusion','dissipation','advection','full pressure','timestep','total')
%%
figure; plot(-sum(tendency,1),'k','linewidth',2); hold all; plot(windstress,'linewidth',2); plot(sum(diffbot,1)-sum(difftop(2:6,:),1),'linewidth',2)
plot(sum(advection,1),'linewidth',2); plot(sum(dissipation,1),'linewidth',2); %plot(sum(pressClin,1)); plot(sum(timestep,1))
plot(-sum(tendency,1)+windstress+sum(diffbot,1)-sum(difftop(2:6,:),1)+sum(advection,1)+sum(dissipation,1)+sum(pressClin,1)+sum(timestep,1),'g--','linewidth',2)
legend('-d/dt','windstress','diffusion','advection','dissipation','total')
title('Vertical sum of vorticity budget, 185m, Euler WAG','fontsize',14)
xlabel('days','fontsize',12)
ylabel('m^2/s, negative speeds up gyre','fontsize',12)
set(gca,'fontsize',12)
axis tight
%% depth-integrated time series
weights=repmat(squeeze(dZ(1:16)),[1 148]);
figure; plot(-sum(weights.*tendency,1),'linewidth',1.5); hold all; %plot(5.*windstress,'linewidth',1.5); 
%plot(sum(weights.*diffbot,1)-sum(weights(:,:).*difftop(:,:),1),'linewidth',1.5)
plot(sum(weights.*diffbot,1)+sum(windstress(:,:).*difftop(:,:),1),'linewidth',1.5)
plot(sum(weights.*advection,1),'linewidth',1.5); plot(sum(weights.*dissipation,1),'linewidth',1.5); %plot(sum(pressClin,1)); plot(sum(timestep,1))
plot(-sum(weights.*tendency,1)+sum(weights.*windstress,1)+sum(weights.*diffbot,1)+sum(weights.*advection,1)+sum(weights.*dissipation,1)+sum(weights.*pressClin,1)+sum(weights.*timestep,1),'k','linewidth',1.5)
legend('-d/dt','windstress and diffusion','advection','dissipation','total')
title('Depth-integrated vorticity budget, 185m, Euler WAG','fontsize',14)
xlabel('days','fontsize',12)
ylabel('m^3/s^2, negative speeds up gyre','fontsize',12)
set(gca,'fontsize',12)
axis tight
%%
figure; plot(-sum(weights.*tendency,1),'k','linewidth',2); hold all;  plot(sum(weights.*diffbot,1)-sum(weights.*difftop,1),'linewidth',2)
plot(sum(weights.*advection,1),'linewidth',2); plot(sum(weights.*dissipation,1),'linewidth',2); %plot(sum(pressClin,1)); plot(sum(timestep,1))
plot(-sum(weights.*tendency,1)+sum(weights.*windstress,1)+sum(weights.*diffbot,1)-sum(weights(2:16,:).*difftop(2:16,:),1)+sum(weights.*advection,1)+sum(weights.*dissipation,1)+sum(weights.*pressClin,1)+sum(weights.*timestep,1),'g--','linewidth',2)
legend('-d\Gamma/dt','wind + v. diffusion','advection','drag + h. diffusion','total')
title('Depth-integrated vorticity budget, 185m, Euler WAG','fontsize',14)
xlabel('days','fontsize',12)
ylabel('m^3/s^2, negative speeds up gyre','fontsize',12)
set(gca,'fontsize',12)
axis tight
%% error
err1=-sum(weights.*tendency,1)+sum(weights.*diffbot,1)-sum(weights.*difftop,1)+sum(weights.*advection,1)+sum(weights.*dissipation,1)+sum(weights.*pressClin,1)+sum(weights.*timestep,1);
perr1=100*err1./sum(weights.*tendency,1);
figure; plot(perr1)
%% depth-integrated time series, 1 week moving mean
L = filter(ones(7,1)/7,1,[-sum(weights.*tendency,1).'; zeros(6,1)]);
out1 = L(7:end);
figure; plot(out1,'k','linewidth',2); hold all; 
L = filter(ones(7,1)/7,1,[5.*windstress.'; zeros(6,1)]);
out2 = L(7:end);
plot(out2,'linewidth',2); 
L = filter(ones(7,1)/7,1,[(sum(weights.*diffbot,1)-sum(weights(2:6,:).*difftop(2:6,:))).'; zeros(6,1)]);
out3 = L(7:end);
plot(out3,'linewidth',2)
L = filter(ones(7,1)/7,1,[-sum(weights.*advection,1).'; zeros(6,1)]);
out4 = L(7:end);
plot(out4,'linewidth',2); 
L = filter(ones(7,1)/7,1,[-sum(weights.*dissipation,1).'; zeros(6,1)]);
out5 = L(7:end);
plot(out5,'linewidth',2); %plot(sum(pressClin,1)); plot(sum(timestep,1))
plot(out1+out2+out3+out4+out5,'g--','linewidth',2)
legend('-d/dt','windstress','diffusion','advection','dissipation','total')
title('Depth-integrated vorticity budget, 1-week moving average','fontsize',14)
xlabel('days','fontsize',12)
ylabel('m^3/s, negative speeds up gyre','fontsize',12)
set(gca,'fontsize',12)
axis tight
%%
allterms(:,:,1)=tendency; allterms(:,:,2)=diffbot-difftop; allterms(:,:,3)=advection; allterms(:,:,4)=dissipation;
figure; bar(1:16, squeeze(mean(allterms,2)))
axis tight
legend('d/dt','vertical diffusion','advection','drag')
title('Mean Vorticity Terms, 16 layers','fontsize',14)
set(gca,'xtick',1:16)
set(gca,'xticklabel',{'0-5m','5-11m','11-17m','17-24m','24-31m','31-39m','39-48m','48-58m','58-69m','69-81m','81-94m','94-109m','109-125m','125-143m','143-163m','163-185m'})
set(gca,'fontsize',12)
ylabel('m^2/s, negative speeds up gyre','fontsize',12)

figure; bar(1:16, squeeze(mean(abs(allterms),2)))
axis tight
legend('d/dt','vertical diffusion','advection','drag')
title('Mean Daily Magnitudes of Vorticity Terms, 16 layers','fontsize',14)
set(gca,'xtick',1:16)
set(gca,'xticklabel',{'0-5m','5-11m','11-17m','17-24m','24-31m','31-39m','39-48m','48-58m','58-69m','69-81m','81-94m','94-109m','109-125m','125-143m','143-163m','163-185m'})
set(gca,'fontsize',12)
ylabel('m^2/s','fontsize',12)
%%
allterms2=squeeze(sum(allterms,1));
[rhos,ps]=corr(allterms2)
weights2=repmat(squeeze(dZ(1:16)),[1 148 4]);
allterms3=squeeze(sum(weights2.*allterms,1));
[rhoi,pi]=corr(allterms3)
%advection and d/dt correlate; diffusion and dissipation too
%advection+d/dt are ~ negative diffusion+dissipation
figure; plot(-sum(weights.*tendency,1),'k','linewidth',2); hold all; plot(sum(weights.*advection,1),'linewidth',2);
 plot(sum(weights.*diffbot,1)-sum(weights.*difftop,1),'linewidth',2)
 plot(sum(weights.*dissipation,1),'linewidth',2); %plot(sum(pressClin,1)); plot(sum(timestep,1))
plot(-sum(weights.*tendency,1)+5.*windstress+sum(weights.*diffbot,1)-sum(weights(2:16,:).*difftop(2:16,:),1)+sum(weights.*advection,1)+sum(weights.*dissipation,1)+sum(weights.*pressClin,1)+sum(weights.*timestep,1),'g--','linewidth',2)
legend('-d/dt','windstress','diffusion','advection','dissipation','total')

%% monthly means
for i=1:4
    figure; bar(1:16, squeeze(mean(allterms(:,1+(i-1)*30:i*30,:),2)))
axis tight
legend('d/dt','vertical diffusion','advection','drag')
t1=strcat('Mean Vorticity Terms, 16 layers, month ',num2str(i));
title(t1,'fontsize',14)
set(gca,'xtick',1:16)
set(gca,'xticklabel',{'0-5m','5-11m','11-17m','17-24m','24-31m','31-39m','39-48m','48-58m','58-69m','69-81m','81-94m','94-109m','109-125m','125-143m','143-163m','163-185m'})
set(gca,'fontsize',12)
ylabel('m^2/s, negative speeds up gyre','fontsize',12)

figure; bar(1:16, squeeze(mean(abs(allterms(:,1+(i-1)*30:i*30,:)),2)))
axis tight
legend('d/dt','vertical diffusion','advection','drag')
t1=strcat('Mean Magnitudes of Vorticity Terms, 16 layers, month ',num2str(i));
title(t1,'fontsize',14)
set(gca,'xtick',1:16)
set(gca,'xticklabel',{'0-5m','5-11m','11-17m','17-24m','24-31m','31-39m','39-48m','48-58m','58-69m','69-81m','81-94m','94-109m','109-125m','125-143m','143-163m','163-185m'})
set(gca,'fontsize',12)
ylabel('m^2/s','fontsize',12)
end

i=5;

  figure; bar(1:16, squeeze(mean(allterms(:,1+(i-1)*30:end,:),2)))
axis tight
legend('d/dt','vertical diffusion','advection','drag')
t1=strcat('Mean Vorticity Terms, 16 layers, month ',num2str(i));
title(t1,'fontsize',14)
set(gca,'xtick',1:16)
set(gca,'xticklabel',{'0-5m','5-11m','11-17m','17-24m','24-31m','31-39m','39-48m','48-58m','58-69m','69-81m','81-94m','94-109m','109-125m','125-143m','143-163m','163-185m'})
set(gca,'fontsize',12)
ylabel('m^2/s, negative speeds up gyre','fontsize',12)

figure; bar(1:16, squeeze(mean(abs(allterms(:,1+(i-1)*30:end,:)),2)))
axis tight
legend('d/dt','vertical diffusion','advection','drag')
t1=strcat('Mean Magnitudes of Vorticity Terms, 16 layers, month ',num2str(i));
title(t1,'fontsize',14)
set(gca,'xtick',1:16)
set(gca,'xticklabel',{'0-5m','5-11m','11-17m','17-24m','24-31m','31-39m','39-48m','48-58m','58-69m','69-81m','81-94m','94-109m','109-125m','125-143m','143-163m','163-185m'})
set(gca,'fontsize',12)
ylabel('m^2/s','fontsize',12)

%%
for i=1:4
    figure; bar(1:4, squeeze(mean(sum(weights2(:,1+(i-1)*30:i*30,:).*allterms(:,1+(i-1)*30:i*30,:),1),2)))
    set(gca,'xtick',1:4)
set(gca,'xticklabel',{'d/dt','diff','adv','diss'})
t1=strcat('mean depth-int terms, month ',num2str(i));
title(t1)

    figure; bar(1:4, squeeze(mean(abs(sum(weights2(:,1+(i-1)*30:i*30,:).*allterms(:,1+(i-1)*30:i*30,:),1)),2)))
    set(gca,'xtick',1:4)
set(gca,'xticklabel',{'d/dt','diff','adv','diss'})
t1=strcat('mean mag depth-int terms, month ',num2str(i));
title(t1)
end
i=5;
    figure; bar(1:4, squeeze(mean(sum(weights2(:,1+(i-1)*30:end,:).*allterms(:,1+(i-1)*30:end,:),1),2)))
    set(gca,'xtick',1:4)
set(gca,'xticklabel',{'d/dt','diff','adv','diss'})
t1=strcat('mean depth-int terms, month ',num2str(i));
title(t1)

    figure; bar(1:4, squeeze(mean(abs(sum(weights2(:,1+(i-1)*30:end,:).*allterms(:,1+(i-1)*30:end,:),1)),2)))
    set(gca,'xtick',1:4)
set(gca,'xticklabel',{'d/dt','diff','adv','diss'})
t1=strcat('mean mag depth-int terms, month ',num2str(i));
title(t1)

%%
weights2=repmat(squeeze(dZ(1:16)),[1 148 4]);
figure; bar(1:16, squeeze(mean(weights2.*allterms,2)))
axis tight
legend('d/dt','vertical diffusion','advection','drag')
title('Mean Vorticity Terms*depth, 16 layers','fontsize',14)
set(gca,'xtick',1:16)
set(gca,'xticklabel',{'0-5m','5-11m','11-17m','17-24m','24-31m','31-39m','39-48m','48-58m','58-69m','69-81m','81-94m','94-109m','109-125m','125-143m','143-163m','163-185m'})
set(gca,'fontsize',12)
ylabel('m^3/s, negative speeds up gyre','fontsize',12)

figure; bar(1:16, squeeze(mean(weights2.*abs(allterms),2)))
axis tight
legend('d/dt','vertical diffusion','advection','drag')
title('Mean Daily Magnitudes*depth of Vorticity Terms, 16 layers','fontsize',14)
set(gca,'xtick',1:16)
set(gca,'xticklabel',{'0-5m','5-11m','11-17m','17-24m','24-31m','31-39m','39-48m','48-58m','58-69m','69-81m','81-94m','94-109m','109-125m','125-143m','143-163m','163-185m'})
set(gca,'fontsize',12)
ylabel('m^3/s','fontsize',12)

%% plot mean vorticity, surface velocity, wind velocity
fc=2*7.2921e-5.*sind(YC);
load('windForEkman.mat', 'tauMean*','l*Cg')
lonCg=lonCg.';
latCg=latCg.';
%load('vorticitySurface.mat', 'vorta')
%load('uvwtsSigma148meanNF.mat', 'UmeanN','VmeanN')
uE=UmeanN(:,:,1).*AngleCS - VmeanN(:,:,1).*AngleSN;  
vN=UmeanN(:,:,1).*AngleSN+ VmeanN(:,:,1).*AngleCS;
xw=XC(1:10:end,1:10:end);
yw=YC(1:10:end,1:10:end);
uE=uE(1:10:end,1:10:end);
vN=vN(1:10:end,1:10:end);
scaleW=0.5;%m/s
scaleA=0.05;%N/m^2, Pascals, Pa
figure; pcolor(XC,YC,mean(vorta(:,:,1:148),3)./fc); shading 'flat'; colorbar
caxis([-1 1])
colormap(b2r(-1,1))
hold all
quiver([-6;lonCg(:)],[37;latCg(:)],[1;tauMean148E(:)/scaleA],[0;tauMean148N(:)/scaleA],'k')
quiver([-6;xw(:)],[37;yw(:)],[0;uE(:)/scaleW],[1;vN(:)/scaleW],'m')
axis([-7 0 35 37.5])
contour(XC,YC,inWag2(:,:,1),[1 1],'linewidth',2,'Color',[0 0.5 0])
legend('\zeta/f','wind','water','WAG')
set(gca,'fontsize',12)
%% plot mean magnitude vs depth
zBin=-0.5*(dInterface(1:16)+dInterface(2:17));

figure; plot(squeeze(mean(abs(allterms),2)),zBin,'linewidth',2);
legend('d/dt','diffusion','advection','dissipation')
axis([0 0.3 -185 0])
set(gca,'fontsize',12)
ylabel('Cell center depth, m')
xlabel('Vorticity budget term magnitudes, m^2/s^2','fontsize',12)
title('Depth structure of vorticity budget, mean magnitudes','fontsize',14)

%% plot some days' terms vs depth 
for i=1:10:101
figure; plot(squeeze(allterms(:,i,:)),zBin,'linewidth',2);
legend('d/dt','diffusion','advection','dissipation')
ylim([-200 0])
set(gca,'fontsize',12)
ylabel('Cell center depth, m')
xlabel('Vorticity budget term magnitudes, m^2/s^2','fontsize',12)
title(num2str(i),'fontsize',14)
end
end%if runthis
