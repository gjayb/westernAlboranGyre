%find cross-manifold flow, the difference between 86400s* daily mean vel
%at today's manifold and the location of the next day's manifold
%similar to crossIsoVel

%% load for surface, fill in days
load('wagAreaAndFluxGateS2filled.mat','*losed','xcm','ycm')
load('geometrySpinupSteady.mat','*C*')
xmin=min(XC(:)); ymin=min(YC(:));
lonWagClosedS=lonWagClosed;
latWagClosedS=latWagClosed;
load('wagAreaAndFluxGate263v2filled.mat','*losed')
for i=[38:42 51 52 80 81 88 109 110 145]
        np=find(latWagClosed(i,:)>0,1,'last');
        latWagClosedS(i,1:np)=latWagClosed(i,1:np);
        lonWagClosedS(i,1:np)=lonWagClosed(i,1:np);
end
load('wagAreaAndFluxGate265v2filled.mat','*losed')
for i=[82:85 87 106:108]
        np=find(latWagClosed(i,:)>0,1,'last');
        latWagClosedS(i,1:np)=latWagClosed(i,1:np);
        lonWagClosedS(i,1:np)=lonWagClosed(i,1:np);
end
latWagClosed=latWagClosedS; clear latWagClosedS
lonWagClosed=lonWagClosedS; clear lonWagClosedS
latWagClosed(latWagClosed==0)=nan;
lonWagClosed(latWagClosed==0)=nan;
 xm=111000*cosd(latWagClosed).*(lonWagClosed-xmin*ones(size(lonWagClosed)));
 ym=111000*(latWagClosed-ymin*ones(size(latWagClosed)));
 xcoast=111000*cosd(latCoast).*(lonCoast-xmin*ones(size(lonCoast)));
 ycoast=111000*(latCoast-ymin*ones(size(latCoast)));
load('uvwSSHDailyDepth1rotated148F.mat', 'Urot','Vrot') %these are ueast, vnorth, at XC,YC

%% surface manifolds
%add first point to end point of xm,ym
[nt,~]=size(xm);
for i=1:nt
    np=find(ym(i,:)>0,1,'last');
    xm(i,np+1)=xm(i,1);
    ym(i,np+1)=ym(i,1);
end
dx=diff(xm,1,2);
dy=diff(ym,1,2);
ds=sqrt(dx.^2+dy.^2);
xc=0.5*xm(:,1:end-1)+0.5*xm(:,2:end);
yc=0.5*ym(:,1:end-1)+0.5*ym(:,2:end);
%is it coast
flagCoast1=ismember(xm,xcoast);
flagCoast2=ismember(ym,ycoast);
flagCoast=flagCoast1&flagCoast2;
coastXYC=flagCoast(:,1:end-1)&flagCoast(:,2:end);%segment is on coast
%surface coast points
for i=1:145
    np=find(xc(i,:)>0,1,'last');
    for k=1:np
        dxhold=xcoast-xc(i,k);
        dyhold=ycoast-yc(i,k);
        [coastDist(i,k),~]=min(sqrt(dxhold.^2+dyhold.^2));

    end
end
figure; plot(xc(20,:),yc(20,:)); hold on
scatter(xc(20,:),yc(20,:),36,coastDist(20,:),'filled')
colorbar
%% u perpendicular
for i=2:nt+1
    uc(i-1,:)=griddata(xcm,ycm,Urot(:,:,i),xc(i-1,:),yc(i-1,:));%indexing issues
    vc(i-1,:)=griddata(xcm,ycm,Vrot(:,:,i),xc(i-1,:),yc(i-1,:)); 
end
vecU=[uc(:).'; vc(:).'];
vecPerp=[(dy(:)./ds(:)).'; (-dx(:)./ds(:)).'];
perpU=reshape(dot(vecU,vecPerp),size(xc));%this is giving negative of what I want
perpU=-perpU;
portionPerp=perpU./sqrt(uc.^2+vc.^2);
%coast removal attempt
perpU2=perpU; perpU2(coastXYC)=nan;
%manifold vel
xcPredict=xc+86400*perpU.*dy./ds;
ycPredict=yc-86400*perpU.*dx./ds;
% lineX1=xc-3*86400*abs(perpU).*dy./ds;
% lineX2=xc+3*86400*abs(perpU).*dy./ds;
% lineY1=yc+3*86400*abs(perpU).*dx./ds;
% lineY2=yc-3*86400*abs(perpU).*dx./ds;
lineX1=xc-20000.*dy./ds;
lineX2=xc+20000.*dy./ds;
lineY1=yc+20000.*dx./ds;
lineY2=yc-20000.*dx./ds;
%% loop through surface for manifold velocity
manDist=zeros(size(xc(2:end,:)));
multipleManDist=zeros(size(xc(2:end,:)));
for i=1:145
    np=find(xc(i,:)>0,1,'last');
    np2=find(xc(i+1,:)>0,1,'last');
    for k=1:np
        [xhold,yhold]=polyxpoly([lineX1(i,k) lineX2(i,k)],[lineY1(i,k) lineY2(i,k)],xc(i+1,1:np2),yc(i+1,1:np2));
        if length(xhold)>0
        dxhold=xhold-xc(i,k);
        dyhold=yhold-yc(i,k);
        [manDist(i,k),index1]=min(sqrt(dxhold.^2+dyhold.^2));
       hold1=cross([dx(i,k) dy(i,k) 0],[dxhold(index1) dyhold(index1) 0]);
        manSign(i,k)=sign(hold1(3));
        xcNext(i,k)=xhold(index1);
        ycNext(i,k)=yhold(index1);
            if length(xhold)>1
                multipleManDist(i,k)=1;
            end
        else
            manDist(i,k)=nan; %to try to remove errors
        end
       
       
    end
    
end
xcNext(xcNext==0)=nan;
ycNext(ycNext==0)=nan;
manDistNoCoast=manDist; manDistNoCoast(coastXYC)=nan;
manDistOld=manDist; manDistOld(isnan(manDistOld))=0;
%% cross-manifold velocity and volume
dz=5;
velMan=manDistNoCoast.*manSign./86400;
velManOld=manDistOld.*manSign./86400;
velErr=perpU(1:145,:)-velManOld;
velErr2=perpU2(1:145,:)-velMan; velErr2(coastXYC)=nan;
%velErr3=-perpU2(1:145,:)-velMan; velErr3(coastXYC)=nan;%fixed already
velCoast=perpU(1:145,:).*double(coastXYC(1:145,:));%-velManOld.*double(coastXYC(1:145,:));coast doesn't move!!
velCoast2=velCoast; velCoast2(~coastXYC(1:145,:))=nan;
coastLength=nansum(ds(1:nt-1,:).*double(coastXYC(1:145,:)),2);
dz=5;%bin depth for now, this may need SSH
volErr=ds(1:nt-1,:).*velErr.*dz;
volErr2=ds(1:nt-1,:).*velErr2.*dz;%no coast
volErrTot=nansum(volErr,2);
volErrTot2=nansum(volErr2,2);
volErrCoast=nansum(ds(1:nt-1,:).*velCoast.*dz,2);
perpVol=ds(1:nt,:).*perpU.*dz;
perpVol2=ds(1:nt,:).*perpU2.*dz;
gateEst=nanmax(perpVol,[],2);
gateEst2=nanmax(perpVol2,[],2);
gateFlag=(xc<4.4e5)&(xc>4.35e5);
gateTot=nansum(volErr2.*double(gateFlag(1:nt-1,:)).*double(~coastXYC(1:nt-1,:)),2);
eastFlag=(xc>5.75e5)&(~coastXYC);
westFlag=(xc<3.54e5)&(~coastXYC);
swFlag=(xc>3.54e5)&(xc<4.35e5)&(yc<3.2e5)&(~coastXYC);
eastTot=nansum(volErr2.*double(eastFlag(1:nt-1,:)),2);
westTot=nansum(volErr2.*double(westFlag(1:nt-1,:)),2);
swTot=nansum(volErr2.*double(swFlag(1:nt-1,:)),2);
multipleIntTot=nansum(volErr2.*multipleManDist(1:nt-1,:).*double(~(swFlag(1:nt-1,:)|eastFlag(1:nt-1,:)|westFlag(1:nt-1,:))),2);



figure; plot(volErrTot-gateEst(1:145)); hold all; plot(volErrTot2-gateEst2(1:145)); plot(volErrTot-gateEst(1:145)-volErrTot2+gateEst2(1:145))
plot(volErrCoast); legend('old estimate','new estimate','difference','coast')
figure;  plot(volErrTot2-gateEst2(1:145)); hold all; plot(volErrTot2-gateTot)
plot(volErrCoast); legend('old gate estimate','new gate estimate','coast')

figure; plot(volErrTot2-gateTot); hold all; plot(eastTot); plot(westTot); plot(swTot); plot(multipleIntTot); 
plot(volErrTot2-gateTot-eastTot-westTot-swTot-multipleIntTot)
legend('err','from far east','far west','sw','multiple ints','remainder')
figure; plot(volErrTot2-gateTot); hold all; 
plot(volErrTot2-gateTot-eastTot-westTot-swTot-multipleIntTot)
%figure; plot(abs(volErrTot-gateEst(1:145))); hold all; plot(abs(volErrTot2-gateEst2(1:145))); legend('old estimate','man vel defined and not coast estimate')
%%
for i=21:40
    figure; plot(xc(i,:),yc(i,:)); hold on
    scatter(xc(i,:),yc(i,:),36,velErr2(i,:),'filled'); 
    caxis([-0.25 0.25]); colorbar; hold all
    plot(xc(i,multipleManDist(i,:)==1),yc(i,multipleManDist(i,:)==1),'m*')
    %plot(xc(i,isnan(manDistNoCoast(i,:))),yc(i,isnan(manDistNoCoast(i,:))),'m*')
    plot(xc(i+1,:),yc(i+1,:));
    title(num2str(i))
end
%%

%%
notCoast1=coastDist>5000;
notCoast2=coastDist>15000;
xLim=(xc<6e5).*(xc>3e5);
volErr3=nansum(volErr.*double(notCoast1),2);
volErr4=nansum(volErr.*double(notCoast2),2);
volErr5=nansum(volErr.*double(notCoast2).*double(xLim(1:145,:)),2);
figure; plot(volErr2); hold all; plot(volErr3); plot(volErr4); plot(volErr5)
filldaysS=[38:42 51 52 80:85 87 88 106:109 110 145];
plot(filldaysS,volErr4(filldaysS),'k*')
%volErr5=volErr4; volErr5(filldaysS)=nan; volErr5(filldaysS-1)=nan;
%plot(volErr5)
%% plots surface
figure; plot(volErr2-gateEst(1:145),'linewidth',2); hold all; plot(nansum(volErr.*double(coastDist<15000),2)-gateEst(1:145),'linewidth',2); 
plot(nansum(volErr.*double(xc(1:145,:)>6e5),2)-gateEst(1:145),'linewidth',2); plot(nansum(volErr.*double(xc(1:145,:)<3e5),2)-gateEst(1:145),'linewidth',2)
plot(volErr5-gateEst(1:145),'k','linewidth',2)
legend('cross-manifold volume flux','manifold near coast','manifold in SoG','manifold far east','remainder')
axis tight
set(gca,'fontsize',12)
xlabel('simulation day')
ylabel('volume flux in top 5m bin, m^3/s')
title('Surface cross-manifold flux')
%%
load('gateAdvectionS263indexed.mat', 'gateFluxBin1')
figure; plot(gateFluxBin1); hold all; plot(max(volErr,[],2))
%%
figure; plot(abs(abs(volErr2-gateFluxBin1(1:145).')./gateFluxBin1(2:146).'),'linewidth',2); hold all; plot(abs(abs(volErr5-gateFluxBin1(1:145).')./gateFluxBin1(2:146).'),'--','linewidth',2); 
axis tight
set(gca,'fontsize',12)
set(gca,'yscale','log')
xlabel('simulation day')
ylabel('volume flux proportion in top 5m bin, unitless')
title('Surface cross-manifold flux (not including gate)/gate flux')
legend('full cross-manifold flux','cross-manifold flux away from coast, in correct WAG region')
%% load for \sigma=26.3



%% 26.3 manifolds


%%

%% surface euler via crossCurveFlux
depth=5+SSHa;
field=ones(size(depth));
[ integral1,vec1 ] = crossCurveFlux( xm,ym,Urot,Vrot,edgeWagX,edgeWagY,field,depth );

for i=1:length(edgeWagX)
[coastDist(i),~]=min(sqrt((edgeWagX(i)-xcoast).^2+(edgeWagY(i)-ycoast).^2));
end

coastEuler=(coastDist(1:end-1)<5000)&(coastDist(2:end)<5000);%segment is 'on coast'
for i=1:148
eulerCrossCoast(i)=nansum(vec1(coastEuler,i));
end
%% randn errors along length of manifold
[na,nb]=size(velErr);
velErrRand=0.02.*(randn(na,nb));
volErrRand=abs(ds(1:nt-1,:)).*velErrRand.*dz;
velErrSize=0.02.*sign(randi(2,na,nb)-1.5);
volErrSize=ds(1:nt-1,:).*velErrSize.*dz;
volErrMax=nansum(0.01.*dz.*abs(ds(1:nt-1,:)),2);
volErrMin=nansum(-0.01.*dz.*abs(ds(1:nt-1,:)),2);
% figure; plot(vec1(repmat(coastEuler',[1 148])))
% hold all; plot(volErrRand(coastXYC(1:nt-1,:)))

volRandCoast=nansum(volErrRand.*double(coastXYC(1:nt-1,:)),2);
volErrSizeC=nansum(volErrSize.*double(coastXYC(1:nt-1,:)),2);
volErrSizeTot=nansum(volErrSize,2);
volErrRandTot=nansum(volErrRand,2);
figure; plot(volErrCoast); hold all; plot(volRandCoast); plot(volErrSizeC);
figure; plot(volErrTot2-gateTot); hold all; plot(volErrSizeTot); plot(volErrRandTot); plot(volErrMax); plot(volErrMin)
velErrMean=(volErrTot2)./nansum(dz.*abs(ds(1:nt-1,:)),2);
velErrCoastMean=(volErrCoast)./nansum(dz.*abs(ds(1:nt-1,:).*coastXYC(1:nt-1,:)),2);
velErrMeanRem=(volErrTot2-gateTot)./nansum(dz.*abs(ds(1:nt-1,:)),2);
velErrMeanRem2=(volErrTot2-(gateTot+multipleIntTot+swTot+eastTot+westTot))./nansum(dz.*abs(ds(1:nt-1,:)),2);
figure; plot(velErrMean); hold all; plot(velErrMeanRem); plot(velErrMeanRem2)
 %%