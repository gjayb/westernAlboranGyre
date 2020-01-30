%finds surface WAG area for each day from manifolds

%load('manifoldsSurfaceInt14NFanalyzed.mat')%often doesn't close
%or load('analyzedSurfaceMan8day.mat') ?
%load('analyzedSurfaceMan14day.mat','*C*','lon*','lat*')
load('geometrySpinupSteady')
[nt,~]=size(lontrFall)

lontrFall(lontrFall==0)=nan;
lontrBall(lontrBall==0)=nan;
lattrFall(lattrFall==0)=nan;
lattrBall(lattrBall==0)=nan;

xmin=min(min(XC)); ymin=min(min(YC));

clear Sigma CT Angle* density* open* wagR


%% find 'initial' points of manifolds
for i=1:nt
    i
    [~,~,ii]=polyxpoly(lonCoast(575:700),latCoast(575:700),lontrFall(i,:),lattrFall(i,:));
    index=find(ii(:,1)==max(ii(:,1)),1,'last');
    if index>0
    coastIndex1(i)=574+ii(index,1);
    np=find(lattrFall(i,:)>0,1,'last');
    lonFall(i,1:np)=[lontrFall(i,ii(index,2):np) lontrFall(i,1:ii(index,2)-1)];
    latFall(i,1:np)=[lattrFall(i,ii(index,2):np) lattrFall(i,1:ii(index,2)-1)];
    end
    
    [xi,yi,ii]=polyxpoly(lonCoast(700:1000),latCoast(700:1000),lontrBall(i,:),lattrBall(i,:));
    index=find(ii(:,1)==min(ii(:,1)),1,'first');
    if index>0
    coastIndex2(i)=700+ii(index,1);
    np=find(lattrBall(i,:)>0,1,'last');
    lonBall(i,1:np)=[lontrBall(i,ii(index,2):-1:1) lontrBall(i,np:-1:ii(index,2)+1)];
    latBall(i,1:np)=[lattrBall(i,ii(index,2):-1:1) lattrBall(i,np:-1:ii(index,2)+1)];
    end
    
    if mod(i,15)==0
        lonFall(latFall==0)=nan;
        latFall(latFall==0)=nan;
        lonBall(latBall==0)=nan;
        latBall(latBall==0)=nan;
        figure; plot(lonCoast(coastIndex1(i):coastIndex2(i)),latCoast(coastIndex1(i):coastIndex2(i)),'k')
        hold on
        plot(lonFall(i,:),latFall(i,:),'r'); plot(lonBall(i,:),latBall(i,:),'b')
        title(num2str(i))
    end
end
disp('section 1 done')

% %% pick gate longitude
% lonOpt=-5.4:0.1:-3.1;
% for j=1:length(lonOpt)
% for i=1:nt
%     index1=find(lonFall(i,:)>lonOpt(j));
%     index2=find(lonBall(i,:)<lonOpt(j));
% 
%     if length(index1)>0 && length(index2)>0
%         gateOpt(j,i)=1;
%     end
% end
% end
% figure; plot(lonOpt,sum(gateOpt,2))
% %% find intersections with the gate

lonWagB=zeros([nt 1]);
lonWagF=lonWagB;
latWagB=lonWagB;
latWagF=lonWagB;
areaWagGate=lonWagB;
lonWagClosed=lonWagB;
latWagClosed=lonWagB;
latGate=zeros([nt 2]);

for i=1:nt
    i
    [~,~,ii]=polyxpoly([-4.4 -4.4],[36.6 35.8],lonFall(i,:),latFall(i,:));
    if length(ii)>1
    index=min(ii(:,2));
    lonWagF(i,1:index+1)=[lonFall(i,1:index) -4.4];
    latWagF(i,1:index+1)=[latFall(i,1:index) latFall(i,index)];
    latGate(i,1)=latFall(i,index);
    end%if intersection forward
    
    [~,~,ii]=polyxpoly([-4.4 -4.4],[36.6 35.8],lonBall(i,:),latBall(i,:));%used to have [36.6 35.3], caused unreasonable results
    if length(ii)>1
    index=min(ii(:,2));
    lonWagB(i,1:index+1)=[lonBall(i,1:index) -4.4];
    latWagB(i,1:index+1)=[latBall(i,1:index) latBall(i,index)];
    latGate(i,2)=latBall(i,index);
    end%if intersection backward
    
    if latWagB(i,1)>0 & latWagF(i,1)>0
        npF=find(latWagF(i,:)>0,1,'last');
        npB=find(latWagB(i,:)>0,1,'last');
        np=length(coastIndex1(i):coastIndex2(i))+npB+npF;

        lonWagClosed(i,1:np)=[lonCoast(coastIndex1(i):coastIndex2(i)).' lonWagB(i,1:npB) lonWagF(i,npF:-1:1)];
        latWagClosed(i,1:np)=[latCoast(coastIndex1(i):coastIndex2(i)).' latWagB(i,1:npB) latWagF(i,npF:-1:1)];
        xm=111000*cosd(latWagClosed(i,1:np)).*(lonWagClosed(i,1:np)-xmin*ones(size(lonWagClosed(i,1:np))));
        ym=111000*(latWagClosed(i,1:np)-ymin*ones(size(latWagClosed(i,1:np))));
        ym(ym==0)=nan;
        xm(ym==0)=nan;
        areaWagGate(i)=polyarea(xm,ym);
        
%         if mod(i,10)==0
%         figure; plot(lonWagClosed(i,1:np),latWagClosed(i,1:np))
%         title(num2str(i))
%         end%if mod(i,10)==0
        
    end %if intersections found for both
end%for i


%%
%from other file, need to update
%% find cells inside WAG each day

load('geometrySpinupSteady.mat')

for i=8:148
    np=find(latWagClosed(i,:)>0,1,'last');
    if length(np)>0
    inside1=inpolygon(XC,YC,lonWagClosed(i,1:np),latWagClosed(i,1:np));
    inWAG(:,:,i)=inside1;
    end
end

% wag volume v1
load('distancesAreas.mat','RAC')
load('isoDepthsNF.mat')
wagVol263=zeros(148,1);
%wagVol27=zeros(148,1);
for i=8:148
        d=max(isoDepth263(:,:,i),5);
        %d2=max(iso27depth(:,:,i),5);
        inthis=inWAG(:,:,i);
        
        for k=1:14
    wagVol263(i)=wagVol263(i)+nansum(nansum((d>dInterface(k)).*(dInterface(k+1)-dInterface(k)).*RAC.*inthis));
    %wagVol27(i)=wagVol27(i)+nansum(nansum((d2>dInterface(k)).*(dInterface(k+1)-dInterface(k)).*RAC.*inthis));

        end
end

%% gate flux: find gate for each day
ymin2=34.5; xmin2=-7.5;
for i=1:nt
    np=find(latWagClosed(i,:)>0,1,'last');
    if length(np)>0 & min(latGate(i,:))>0
        yg1=111000*(min(latGate(i,:))-ymin);
        yg2=111000*(max(latGate(i,:))-ymin);
        yG=[yg1:1000:yg2 yg2]; np2=length(yG);
        yGate(i,1:np2)=yG;
        gateFluxSign(i)=sign(latGate(i,2)-latGate(i,1));
    end
    
end
lonGate=-4.4*ones(size(latGate));
xGate=111000*cosd(36.3).*(-4.4-xmin)*ones(size(yGate));
xGate(yGate==0)=nan; yGate(yGate==0)=nan;
%% gate flux: load u, interpolate
load('uvwSSHDailyDepth1rotated148F.mat', 'Urot','XC','YC')
u=Urot; clear Urot
xu=111000*cosd(YC).*(XC-xmin);
yu=111000*(YC-ymin);

%% gate flux: integrate, just top 5m
for i=8:nt
    if yGate(i,1)>0
        np=find(yGate(i,:)>0,1,'last');
        uGate=griddata(xu,yu,u(:,:,i),xGate(i,1:np),yGate(i,1:np));
        uGate2=0.5.*uGate(2:end)+0.5.*uGate(1:end-1);
        gateFluxS(i)=5*gateFluxSign(i)*sum(abs(diff(yGate(i,1:np))).*uGate2);
    end
end

figure; plot(gateFluxS,'o');

%% check flux vs area changes, top 5m
holdvar=areaWagGate;
iNan=find(holdvar==0|isnan(holdvar));
holdvar(iNan)=nan;

wagVolGate=zeros(size(wagVol263));
for i=8:148
    if ismember(i,iNan)
        wagVolGate(i)=nan;
    elseif ismember(i-1,iNan)
        wagVolGate(i)=5*holdvar(i);
    else
        wagVolGate(i)=wagVolGate(i-1)+86400*gateFluxS(i);
    end
end

figure; plot(1:148,5*holdvar,'bo')
hold on
plot(1:148,5*squeeze(sum(sum(inWAG.*repmat(RAC,[1 1 148])))),'ro')
plot(1:148,wagVolGate,'c')

%% gate flux: integrate surface to sigma=26.3

load('uvwSSHDailyDepth1rotated148F.mat','SSHa')
ssh148=SSHa; clear SSHa
load('isoDepthsNF.mat', 'isoDepth263')
xcm=111000*cosd(YC).*(XC-xmin*ones(size(XC)));
ycm=111000*(YC-ymin*ones(size(YC)));
load('gateLocations.mat')

for i=8:nt
    if yGate(i,1)>0
        np=find(yGate(i,:)>0,1,'last');
    end
    if ismember(i,[38:42 51 52 80 81 88 109 110 145])
        np=find(yGate263(i,:)>0,1,'last');
        xGate(i,1:np)=xGate263(i,1:np);
        yGate(i,1:np)=yGate263(i,1:np);
        gateFluxSign(i)=gateFluxSign263(i);
    elseif ismember(i,[82:85 87 106:108])
        np=find(yGate265(i,:)>0,1,'last');
        xGate(i,1:np)=xGate265(i,1:np);
        yGate(i,1:np)=yGate265(i,1:np);
        gateFluxSign(i)=gateFluxSign265(i);
    end
        uGate=griddata(xu,yu,u(:,:,i),xGate(i,1:np),yGate(i,1:np));
        dGate=griddata(xcm,ycm,isoDepth263(:,:,i),xGate(i,1:np),yGate(i,1:np));
        hGate=griddata(xcm,ycm,ssh148(:,:,i),xGate(i,1:np),yGate(i,1:np));
        uGate2=0.5.*uGate(2:end)+0.5.*uGate(1:end-1);
        meanUS(i)=mean(uGate)*gateFluxSign(i);
        dGate2=0.5.*dGate(2:end)+0.5.*dGate(1:end-1);
        hGate2=0.5.*hGate(2:end)+0.5.*hGate(1:end-1);
        meanDS(i)=mean(dGate+hGate);
        gateFluxS26(i)=gateFluxSign(i)*sum(abs(diff(yGate(i,1:np))).*uGate2.*(hGate2+dGate2));

end

figure; plot(gateFluxS26,'o');
%%
gatelon=-4.4;
clear u*
save('wagAreaAndFluxGateS2filled.mat')
%% sigma=26 manifolds, gate, and gate flux sign
%here! 11/13/2017
clear; %%had problems before w/o this clear!!
load('manifoldsIso263Int14NFanalyzed.mat')
[nt,~]=size(lontrFall)

lontrFall(lontrFall==0)=nan;
lontrBall(lontrBall==0)=nan;
lattrFall(lattrFall==0)=nan;
lattrBall(lattrBall==0)=nan;
load('geometrySpinupSteady.mat')
xmin=min(min(XC)); ymin=min(min(YC));
lonWagB=zeros([nt 1]);
lonWagF=lonWagB;
latWagB=lonWagB;
latWagF=lonWagB;
areaWagGate=lonWagB;
lonWagClosed=lonWagB;
latWagClosed=lonWagB;
latGate=zeros([nt 2]);

for i=14:nt
    i
    [~,~,ii]=polyxpoly(lonCoast(575:700),latCoast(575:700),lontrFall(i,:),lattrFall(i,:));
    index=find(ii(:,1)==max(ii(:,1)),1,'last');
    coastIndex1(i)=574+ii(index,1);
    np=find(lattrFall(i,:)>0,1,'last');
    lonFall(i,1:np)=[lontrFall(i,ii(index,2):np) lontrFall(i,1:ii(index,2)-1)];
    latFall(i,1:np)=[lattrFall(i,ii(index,2):np) lattrFall(i,1:ii(index,2)-1)];
    
    [xi,yi,ii]=polyxpoly(lonCoast(700:1000),latCoast(700:1000),lontrBall(i,:),lattrBall(i,:));
    index=find(ii(:,1)==min(ii(:,1)),1,'first');
    coastIndex2(i)=700+ii(index,1);
    np=find(lattrBall(i,:)>0,1,'last');
    lonBall(i,1:np)=[lontrBall(i,ii(index,2):-1:1) lontrBall(i,np:-1:ii(index,2)+1)];
    latBall(i,1:np)=[lattrBall(i,ii(index,2):-1:1) lattrBall(i,np:-1:ii(index,2)+1)];

    [~,~,ii]=polyxpoly([-4.4 -4.4],[36.6 35.8],lonFall(i,:),latFall(i,:));
    if length(ii)>1
    index=min(ii(:,2));
    lonWagF(i,1:index+1)=[lonFall(i,1:index) -4.4];
    latWagF(i,1:index+1)=[latFall(i,1:index) latFall(i,index)];
    latGate(i,1)=latFall(i,index);
    end%if intersection forward
    
    [~,~,ii]=polyxpoly([-4.4 -4.4],[36.6 35.8],lonBall(i,:),latBall(i,:));
    if length(ii)>1
    index=min(ii(:,2));
    lonWagB(i,1:index+1)=[lonBall(i,1:index) -4.4];
    latWagB(i,1:index+1)=[latBall(i,1:index) latBall(i,index)];
    latGate(i,2)=latBall(i,index);
    end%if intersection backward
    
    if latWagB(i,1)>0 & latWagF(i,1)>0
        npF=find(latWagF(i,:)>0,1,'last');
        npB=find(latWagB(i,:)>0,1,'last');
        np=length(coastIndex1(i):coastIndex2(i))+npB+npF;

        lonWagClosed(i,1:np)=[lonCoast(coastIndex1(i):coastIndex2(i)).' lonWagB(i,1:npB) lonWagF(i,npF:-1:1)];
        latWagClosed(i,1:np)=[latCoast(coastIndex1(i):coastIndex2(i)).' latWagB(i,1:npB) latWagF(i,npF:-1:1)];
        xm=111000*cosd(latWagClosed(i,1:np)).*(lonWagClosed(i,1:np)-xmin*ones(size(lonWagClosed(i,1:np))));
        ym=111000*(latWagClosed(i,1:np)-ymin*ones(size(latWagClosed(i,1:np))));
        ym(ym==0)=nan;
        xm(ym==0)=nan;
        areaWagGate(i)=polyarea(xm,ym);
                
        
        
        np=find(latWagClosed(i,:)>0,1,'last');
    if length(np)>0 & min(latGate(i,:))>0
        yg1=111000*(min(latGate(i,:))-ymin);
        yg2=111000*(max(latGate(i,:))-ymin);
        yG=[yg1:1000:yg2 yg2]; np2=length(yG);
        yGate(i,1:np2)=yG;
        gateFluxSign(i)=sign(latGate(i,2)-latGate(i,1));
    end

    end %if intersections found for both
    
end

lonGate=-4.4*ones(size(latGate));
xGate=111000*cosd(36.3).*(-4.4-xmin)*ones(size(yGate));
xGate(yGate==0)=nan; yGate(yGate==0)=nan;
%% gate flux sigma 26.3 to 26.5
 
load('uvwNativeGridIsoDepth263.mat','uIso')
load('isoDepthsNF.mat', 'isoDepth263')
load('isoDepthsNF.mat', 'isoDepth265')
xcm=111000*cosd(YC).*(XC-xmin*ones(size(XC)));
ycm=111000*(YC-ymin*ones(size(YC)));
for i=14:nt
    if yGate(i,1)>0
        np=find(yGate(i,:)>0,1,'last');
    end
    if ismember(i,[60:63 71])
        np=find(yGateS(i,:)>0,1,'last');
        xGate(i,1:np)=xGateS(i,1:np);
        yGate(i,1:np)=yGateS(i,1:np);
        gateFluxSign(i)=gateFluxSignS(i);
    elseif ismember(i,[82:87 106:108])
        np=find(yGate265(i,:)>0,1,'last');
        xGate(i,1:np)=xGate265(i,1:np);
        yGate(i,1:np)=yGate265(i,1:np);
        gateFluxSign(i)=gateFluxSign265(i);
    end
        uGate=griddata(xcm,ycm,uIso(:,:,i),xGate(i,1:np),yGate(i,1:np));
        uGate2=0.5.*uGate(2:end)+0.5.*uGate(1:end-1);
         meanUg26(i)=gateFluxSign(i)*mean(uGate);
        
        dGate26=griddata(xcm,ycm,isoDepth263(:,:,i),xGate(i,1:np),yGate(i,1:np));
        dGate262=0.5.*dGate26(2:end)+0.5.*dGate26(1:end-1);
        dGate265=griddata(xcm,ycm,isoDepth265(:,:,i),xGate(i,1:np),yGate(i,1:np));
        dGate2652=0.5.*dGate265(2:end)+0.5.*dGate265(1:end-1);
        meanD26265(i)=mean(dGate265-dGate26);
        
        gateFlux26265(i)=gateFluxSign(i)*sum(abs(diff(yGate(i,1:np))).*uGate2.*(dGate2652-dGate262));
    
end

figure; plot(gateFlux26265,'o--'); title('Gate Flux, \sigma=26.3 to \sigma=26.5')
%%
clear u*
save('wagAreaAndFluxGate263v2filled.mat')
%% sigma=26.5 manifolds, gate, and gate flux sign
%here 11/26
clear;
load('manifoldsIso265Int14NFanalyzed.mat')
[nt,~]=size(lontrFall)

lontrFall(lontrFall==0)=nan;
lontrBall(lontrBall==0)=nan;
lattrFall(lattrFall==0)=nan;
lattrBall(lattrBall==0)=nan;
load('geometrySpinupSteady.mat')
xmin=min(min(XC)); ymin=min(min(YC));
lonWagB=zeros([nt 1]);
lonWagF=lonWagB;
latWagB=lonWagB;
latWagF=lonWagB;
areaWagGate=lonWagB;
lonWagClosed=lonWagB;
latWagClosed=lonWagB;
latGate=zeros([nt 2]);

for i=14:nt
    i
    [~,~,ii]=polyxpoly(lonCoast(575:700),latCoast(575:700),lontrFall(i,:),lattrFall(i,:));
    index=find(ii(:,1)==max(ii(:,1)),1,'last');
    coastIndex1(i)=574+ii(index,1);
    np=find(lattrFall(i,:)>0,1,'last');
    lonFall(i,1:np)=[lontrFall(i,ii(index,2):np) lontrFall(i,1:ii(index,2)-1)];
    latFall(i,1:np)=[lattrFall(i,ii(index,2):np) lattrFall(i,1:ii(index,2)-1)];
    
    [xi,yi,ii]=polyxpoly(lonCoast(700:1000),latCoast(700:1000),lontrBall(i,:),lattrBall(i,:));
    index=find(ii(:,1)==min(ii(:,1)),1,'first');
    coastIndex2(i)=700+ii(index,1);
    np=find(lattrBall(i,:)>0,1,'last');
    lonBall(i,1:np)=[lontrBall(i,ii(index,2):-1:1) lontrBall(i,np:-1:ii(index,2)+1)];
    latBall(i,1:np)=[lattrBall(i,ii(index,2):-1:1) lattrBall(i,np:-1:ii(index,2)+1)];

    [~,~,ii]=polyxpoly([-4.4 -4.4],[36.6 35.8],lonFall(i,:),latFall(i,:));
    if length(ii)>1
    index=min(ii(:,2));
    lonWagF(i,1:index+1)=[lonFall(i,1:index) -4.4];
    latWagF(i,1:index+1)=[latFall(i,1:index) latFall(i,index)];
    latGate(i,1)=latFall(i,index);
    end%if intersection forward
    
    [~,~,ii]=polyxpoly([-4.4 -4.4],[36.6 35.8],lonBall(i,:),latBall(i,:));
    if length(ii)>1
    index=min(ii(:,2));
    lonWagB(i,1:index+1)=[lonBall(i,1:index) -4.4];
    latWagB(i,1:index+1)=[latBall(i,1:index) latBall(i,index)];
    latGate(i,2)=latBall(i,index);
    end%if intersection backward
    
    if latWagB(i,1)>0 & latWagF(i,1)>0
        npF=find(latWagF(i,:)>0,1,'last');
        npB=find(latWagB(i,:)>0,1,'last');
        np=length(coastIndex1(i):coastIndex2(i))+npB+npF;

        lonWagClosed(i,1:np)=[lonCoast(coastIndex1(i):coastIndex2(i)).' lonWagB(i,1:npB) lonWagF(i,npF:-1:1)];
        latWagClosed(i,1:np)=[latCoast(coastIndex1(i):coastIndex2(i)).' latWagB(i,1:npB) latWagF(i,npF:-1:1)];
        xm=111000*cosd(latWagClosed(i,1:np)).*(lonWagClosed(i,1:np)-xmin*ones(size(lonWagClosed(i,1:np))));
        ym=111000*(latWagClosed(i,1:np)-ymin*ones(size(latWagClosed(i,1:np))));
        ym(ym==0)=nan;
        xm(ym==0)=nan;
        areaWagGate(i)=polyarea(xm,ym);
                
        
        
        np=find(latWagClosed(i,:)>0,1,'last');
    if length(np)>0 & min(latGate(i,:))>0
        yg1=111000*(min(latGate(i,:))-ymin);
        yg2=111000*(max(latGate(i,:))-ymin);
        yG=[yg1:1000:yg2 yg2]; np2=length(yG);
        yGate(i,1:np2)=yG;
        gateFluxSign(i)=sign(latGate(i,2)-latGate(i,1));
    end

    end %if intersections found for both
    
end

lonGate=-4.4*ones(size(latGate));
xGate=111000*cosd(36.3).*(-4.4-xmin)*ones(size(yGate));
xGate(yGate==0)=nan; yGate(yGate==0)=nan;
%% gate flux sigma 26.5 to 26.75
clear u*
load('uvwNativeGridIsoDepth265.mat','uIso')
load('isoDepthsNF.mat', 'isoDepth265')
load('isoDepthsNF.mat', 'isoDepth2675')
xcm=111000*cosd(YC).*(XC-xmin*ones(size(XC)));
ycm=111000*(YC-ymin*ones(size(YC)));
for i=14:nt
    if yGate(i,1)>0
        np=find(yGate(i,:)>0,1,'last');
    end
    if ismember(i,[60:63 71])
        np=find(yGateS(i,:)>0,1,'last');
        xGate(i,1:np)=xGateS(i,1:np);
        yGate(i,1:np)=yGateS(i,1:np);
        gateFluxSign(i)=gateFluxSignS(i);
    elseif ismember(i,[64 66:70 103])
        np=find(yGate263(i,:)>0,1,'last');
        xGate(i,1:np)=xGate263(i,1:np);
        yGate(i,1:np)=yGate263(i,1:np);
        gateFluxSign(i)=gateFluxSign263(i);
    end 
        uGate=griddata(xcm,ycm,uIso(:,:,i),xGate(i,1:np),yGate(i,1:np));
        uGate2=0.5.*uGate(2:end)+0.5.*uGate(1:end-1);
         meanUg265(i)=gateFluxSign(i)*mean(uGate);
        
        dGate265=griddata(xcm,ycm,isoDepth265(:,:,i),xGate(i,1:np),yGate(i,1:np));
        dGate2652=0.5.*dGate265(2:end)+0.5.*dGate265(1:end-1);
        dGate27=griddata(xcm,ycm,isoDepth2675(:,:,i),xGate(i,1:np),yGate(i,1:np));
        dGate272=0.5.*dGate27(2:end)+0.5.*dGate27(1:end-1);
        meanD2652675(i)=mean(dGate27-dGate265);
        
        gateFlux2652675(i)=gateFluxSign(i)*sum(abs(diff(yGate(i,1:np))).*uGate2.*(dGate272-dGate2652));
    
end

figure; plot(gateFlux2652675,'o'); title('Gate Flux, \sigma=26.5 to \sigma=26.75')
%%
clear u*
save('wagAreaAndFluxGate265v2filled.mat')

%% sigma=26.75 manifolds, gate, and gate flux sign
%here 11/27/2017
clear;
load('manifoldsIso2675Int14NFanalyzed.mat')
[nt,~]=size(lontrFall)

lontrFall(lontrFall==0)=nan;
lontrBall(lontrBall==0)=nan;
lattrFall(lattrFall==0)=nan;
lattrBall(lattrBall==0)=nan;
load('geometrySpinupSteady.mat')
xmin=min(min(XC)); ymin=min(min(YC));
lonWagB=zeros([nt 1]);
lonWagF=lonWagB;
latWagB=lonWagB;
latWagF=lonWagB;
areaWagGate=lonWagB;
lonWagClosed=lonWagB;
latWagClosed=lonWagB;
latGate=zeros([nt 2]);

for i=14:nt
    i
    [~,~,ii]=polyxpoly(lonCoast(575:700),latCoast(575:700),lontrFall(i,:),lattrFall(i,:));
    index=find(ii(:,1)==max(ii(:,1)),1,'last');
    coastIndex1(i)=574+ii(index,1);
    np=find(lattrFall(i,:)>0,1,'last');
    lonFall(i,1:np)=[lontrFall(i,ii(index,2):np) lontrFall(i,1:ii(index,2)-1)];
    latFall(i,1:np)=[lattrFall(i,ii(index,2):np) lattrFall(i,1:ii(index,2)-1)];
    
    [xi,yi,ii]=polyxpoly(lonCoast(700:1000),latCoast(700:1000),lontrBall(i,:),lattrBall(i,:));
    index=find(ii(:,1)==min(ii(:,1)),1,'first');
    coastIndex2(i)=700+ii(index,1);
    np=find(lattrBall(i,:)>0,1,'last');
    lonBall(i,1:np)=[lontrBall(i,ii(index,2):-1:1) lontrBall(i,np:-1:ii(index,2)+1)];
    latBall(i,1:np)=[lattrBall(i,ii(index,2):-1:1) lattrBall(i,np:-1:ii(index,2)+1)];

    [~,~,ii]=polyxpoly([-4.4 -4.4],[36.6 35.8],lonFall(i,:),latFall(i,:));
    if length(ii)>1
    index=min(ii(:,2));
    lonWagF(i,1:index+1)=[lonFall(i,1:index) -4.4];
    latWagF(i,1:index+1)=[latFall(i,1:index) latFall(i,index)];
    latGate(i,1)=latFall(i,index);
    end%if intersection forward
    
    [~,~,ii]=polyxpoly([-4.4 -4.4],[36.6 35.8],lonBall(i,:),latBall(i,:));
    if length(ii)>1
    index=min(ii(:,2));
    lonWagB(i,1:index+1)=[lonBall(i,1:index) -4.4];
    latWagB(i,1:index+1)=[latBall(i,1:index) latBall(i,index)];
    latGate(i,2)=latBall(i,index);
    end%if intersection backward
    
    if latWagB(i,1)>0 & latWagF(i,1)>0
        npF=find(latWagF(i,:)>0,1,'last');
        npB=find(latWagB(i,:)>0,1,'last');
        np=length(coastIndex1(i):coastIndex2(i))+npB+npF;

        lonWagClosed(i,1:np)=[lonCoast(coastIndex1(i):coastIndex2(i)).' lonWagB(i,1:npB) lonWagF(i,npF:-1:1)];
        latWagClosed(i,1:np)=[latCoast(coastIndex1(i):coastIndex2(i)).' latWagB(i,1:npB) latWagF(i,npF:-1:1)];
        xm=111000*cosd(latWagClosed(i,1:np)).*(lonWagClosed(i,1:np)-xmin*ones(size(lonWagClosed(i,1:np))));
        ym=111000*(latWagClosed(i,1:np)-ymin*ones(size(latWagClosed(i,1:np))));
        ym(ym==0)=nan;
        xm(ym==0)=nan;
        areaWagGate(i)=polyarea(xm,ym);
                
        
        
        np=find(latWagClosed(i,:)>0,1,'last');
    if length(np)>0 & min(latGate(i,:))>0
        yg1=111000*(min(latGate(i,:))-ymin);
        yg2=111000*(max(latGate(i,:))-ymin);
        yG=[yg1:1000:yg2 yg2]; np2=length(yG);
        yGate(i,1:np2)=yG;
        gateFluxSign(i)=sign(latGate(i,2)-latGate(i,1));
    end

    end %if intersections found for both
    
end

lonGate=-4.4*ones(size(latGate));
xGate=111000*cosd(36.3).*(-4.4-xmin)*ones(size(yGate));
xGate(yGate==0)=nan; yGate(yGate==0)=nan;
%% gate flux sigma 26.75 to 27
clear u*
load('uvwNativeGridIsoDepth2675.mat','uIso')%native grid
load('isoDepthsNF.mat', 'isoDepth2675')
load('isoDepthsNF.mat', 'isoDepth27')
xcm=111000*cosd(YC).*(XC-xmin*ones(size(XC)));
ycm=111000*(YC-ymin*ones(size(YC)));
for i=14:nt
    if yGate(i,1)>0
        np=find(yGate(i,:)>0,1,'last');
    end
    if ismember(i,[62 63])
        np=find(yGateS(i,:)>0,1,'last');
        xGate(i,1:np)=xGateS(i,1:np);
        yGate(i,1:np)=yGateS(i,1:np);
        gateFluxSign(i)=gateFluxSignS(i);
    elseif ismember(i,67:70)
        np=find(yGate263(i,:)>0,1,'last');
        xGate(i,1:np)=xGate263(i,1:np);
        yGate(i,1:np)=yGate263(i,1:np);
        gateFluxSign(i)=gateFluxSign263(i);
    elseif ismember(i,[15 43 48 51:56 72 112])
        np=find(yGate265(i,:)>0,1,'last');
        xGate(i,1:np)=xGate265(i,1:np);
        yGate(i,1:np)=yGate265(i,1:np);
        gateFluxSign(i)=gateFluxSign265(i);
    elseif ismember(i,[61 64])
        np=find(yGate27(i,:)>0,1,'last');
        xGate(i,1:np)=xGate27(i,1:np);
        yGate(i,1:np)=yGate27(i,1:np);
        gateFluxSign(i)=gateFluxSign27(i);
    end 
        uGate=griddata(xcm,ycm,uIso(:,:,i),xGate(i,1:np),yGate(i,1:np));
        uGate2=0.5.*uGate(2:end)+0.5.*uGate(1:end-1);
        meanUg2675(i)=gateFluxSign(i)*mean(uGate);
        
        dGate27=griddata(xcm,ycm,isoDepth27(:,:,i),xGate(i,1:np),yGate(i,1:np));
        dGate2752=0.5.*dGate27(2:end)+0.5.*dGate27(1:end-1);
        dGate2675=griddata(xcm,ycm,isoDepth2675(:,:,i),xGate(i,1:np),yGate(i,1:np));
        dGate272=0.5.*dGate2675(2:end)+0.5.*dGate2675(1:end-1);
        meanD267527(i)=mean(dGate27-dGate2675);
        
        gateFlux267527(i)=gateFluxSign(i)*sum(abs(diff(yGate(i,1:np))).*uGate2.*(dGate2752-dGate272));
end

%SO- this flux is in m^3/s, calculated each day. the total VOLUME that goes
%is 86400s/day times this thing. DO I WANT THAT?????

figure; plot(gateFlux267527,'o--'); title('Gate Flux, \sigma=26.75 to \sigma=27')
%%
clear uIso
save('wagAreaAndFluxGate2675v2filled.mat')
%% sigma 27 to 27.5
clear;
load('manifoldsIso27Int14NFanalyzed.mat')
[nt,~]=size(lontrFall)

lontrFall(lontrFall==0)=nan;
lontrBall(lontrBall==0)=nan;
lattrFall(lattrFall==0)=nan;
lattrBall(lattrBall==0)=nan;
load('geometrySpinupSteady.mat')
xmin=min(min(XC)); ymin=min(min(YC));
lonWagB=zeros([nt 1]);
lonWagF=lonWagB;
latWagB=lonWagB;
latWagF=lonWagB;
areaWagGate=lonWagB;
lonWagClosed=lonWagB;
latWagClosed=lonWagB;
latGate=zeros([nt 2]);

for i=14:nt
    i
    [~,~,ii]=polyxpoly(lonCoast(575:700),latCoast(575:700),lontrFall(i,:),lattrFall(i,:));
    index=find(ii(:,1)==max(ii(:,1)),1,'last');
    coastIndex1(i)=574+ii(index,1);
    np=find(lattrFall(i,:)>0,1,'last');
    lonFall(i,1:np)=[lontrFall(i,ii(index,2):np) lontrFall(i,1:ii(index,2)-1)];
    latFall(i,1:np)=[lattrFall(i,ii(index,2):np) lattrFall(i,1:ii(index,2)-1)];
    
    [xi,yi,ii]=polyxpoly(lonCoast(700:1000),latCoast(700:1000),lontrBall(i,:),lattrBall(i,:));
    index=find(ii(:,1)==min(ii(:,1)),1,'first');
    coastIndex2(i)=700+ii(index,1);
    np=find(lattrBall(i,:)>0,1,'last');
    lonBall(i,1:np)=[lontrBall(i,ii(index,2):-1:1) lontrBall(i,np:-1:ii(index,2)+1)];
    latBall(i,1:np)=[lattrBall(i,ii(index,2):-1:1) lattrBall(i,np:-1:ii(index,2)+1)];

    [~,~,ii]=polyxpoly([-4.4 -4.4],[36.6 35.8],lonFall(i,:),latFall(i,:));
    if length(ii)>1
    index=min(ii(:,2));
    lonWagF(i,1:index+1)=[lonFall(i,1:index) -4.4];
    latWagF(i,1:index+1)=[latFall(i,1:index) latFall(i,index)];
    latGate(i,1)=latFall(i,index);
    end%if intersection forward
    
    [~,~,ii]=polyxpoly([-4.4 -4.4],[36.6 35.8],lonBall(i,:),latBall(i,:));
    if length(ii)>1
    index=min(ii(:,2));
    lonWagB(i,1:index+1)=[lonBall(i,1:index) -4.4];
    latWagB(i,1:index+1)=[latBall(i,1:index) latBall(i,index)];
    latGate(i,2)=latBall(i,index);
    end%if intersection backward
    
    if latWagB(i,1)>0 & latWagF(i,1)>0
        npF=find(latWagF(i,:)>0,1,'last');
        npB=find(latWagB(i,:)>0,1,'last');
        np=length(coastIndex1(i):coastIndex2(i))+npB+npF;

        lonWagClosed(i,1:np)=[lonCoast(coastIndex1(i):coastIndex2(i)).' lonWagB(i,1:npB) lonWagF(i,npF:-1:1)];
        latWagClosed(i,1:np)=[latCoast(coastIndex1(i):coastIndex2(i)).' latWagB(i,1:npB) latWagF(i,npF:-1:1)];
        xm=111000*cosd(latWagClosed(i,1:np)).*(lonWagClosed(i,1:np)-xmin*ones(size(lonWagClosed(i,1:np))));
        ym=111000*(latWagClosed(i,1:np)-ymin*ones(size(latWagClosed(i,1:np))));
        ym(ym==0)=nan;
        xm(ym==0)=nan;
        areaWagGate(i)=polyarea(xm,ym);
                
        
        
        np=find(latWagClosed(i,:)>0,1,'last');
    if length(np)>0 & min(latGate(i,:))>0
        yg1=111000*(min(latGate(i,:))-ymin);
        yg2=111000*(max(latGate(i,:))-ymin);
        yG=[yg1:1000:yg2 yg2]; np2=length(yG);
        yGate(i,1:np2)=yG;
        gateFluxSign(i)=sign(latGate(i,2)-latGate(i,1));
    end

    end %if intersections found for both
    
end

lonGate=-4.4*ones(size(latGate));
xGate=111000*cosd(36.3).*(-4.4-xmin)*ones(size(yGate));
xGate(yGate==0)=nan; yGate(yGate==0)=nan;
%% gate flux sigma 27 to 27.5
clear u*
load('uvwNativeGridIsoDepth27.mat','uIso')%native grid
load('isoDepthsNF.mat', 'isoDepth275')
load('isoDepthsNF.mat', 'isoDepth27')
xcm=111000*cosd(YC).*(XC-xmin*ones(size(XC)));
ycm=111000*(YC-ymin*ones(size(YC)));
for i=14:nt
    if yGate(i,1)>0
        np=find(yGate(i,:)>0,1,'last');
    end
    if ismember(i,[62 63])
        np=find(yGateS(i,:)>0,1,'last');
        xGate(i,1:np)=xGateS(i,1:np);
        yGate(i,1:np)=yGateS(i,1:np);
        gateFluxSign(i)=gateFluxSignS(i);
    elseif ismember(i,67:70)
        np=find(yGate263(i,:)>0,1,'last');
        xGate(i,1:np)=xGate263(i,1:np);
        yGate(i,1:np)=yGate263(i,1:np);
        gateFluxSign(i)=gateFluxSign263(i);
    elseif ismember(i,[71 117])
        np=find(yGate2675(i,:)>0,1,'last');
        xGate(i,1:np)=xGate2675(i,1:np);
        yGate(i,1:np)=yGate2675(i,1:np);
        gateFluxSign(i)=gateFluxSign2675(i);
    end 
        uGate=griddata(xcm,ycm,uIso(:,:,i),xGate(i,1:np),yGate(i,1:np));
        uGate2=0.5.*uGate(2:end)+0.5.*uGate(1:end-1);
        meanUg27(i)=gateFluxSign(i)*mean(uGate);
        
        dGate275=griddata(xcm,ycm,isoDepth275(:,:,i),xGate(i,1:np),yGate(i,1:np));
        dGate2752=0.5.*dGate275(2:end)+0.5.*dGate275(1:end-1);
        dGate27=griddata(xcm,ycm,isoDepth27(:,:,i),xGate(i,1:np),yGate(i,1:np));
        dGate272=0.5.*dGate27(2:end)+0.5.*dGate27(1:end-1);
        meanD27275(i)=mean(dGate275-dGate27);
        
        gateFlux27275(i)=gateFluxSign(i)*sum(abs(diff(yGate(i,1:np))).*uGate2.*(dGate2752-dGate272));
end

%SO- this flux is in m^3/s, calculated each day. the total VOLUME that goes
%is 86400s/day times this thing. DO I WANT THAT?????

figure; plot(gateFlux27275,'o--'); title('Gate Flux, \sigma=27 to \sigma=27.5')
%%
clear uIso
save('wagAreaAndFluxGate27v2filled.mat')
%% sigma 27.5
clear;
load('manifoldsIso275Int14NFanalyzed.mat')
[nt,~]=size(lontrFall)

lontrFall(lontrFall==0)=nan;
lontrBall(lontrBall==0)=nan;
lattrFall(lattrFall==0)=nan;
lattrBall(lattrBall==0)=nan;
load('geometrySpinupSteady.mat')
xmin=min(min(XC)); ymin=min(min(YC));
lonWagB=zeros([nt 1]);
lonWagF=lonWagB;
latWagB=lonWagB;
latWagF=lonWagB;
areaWagGate=lonWagB;
lonWagClosed=lonWagB;
latWagClosed=lonWagB;
latGate=zeros([nt 2]);

for i=14:nt
    i
    [~,~,ii]=polyxpoly(lonCoast(575:700),latCoast(575:700),lontrFall(i,:),lattrFall(i,:));
    index=find(ii(:,1)==max(ii(:,1)),1,'last');
    coastIndex1(i)=574+ii(index,1);
    np=find(lattrFall(i,:)>0,1,'last');
    lonFall(i,1:np)=[lontrFall(i,ii(index,2):np) lontrFall(i,1:ii(index,2)-1)];
    latFall(i,1:np)=[lattrFall(i,ii(index,2):np) lattrFall(i,1:ii(index,2)-1)];
    
    [xi,yi,ii]=polyxpoly(lonCoast(700:1000),latCoast(700:1000),lontrBall(i,:),lattrBall(i,:));
    index=find(ii(:,1)==min(ii(:,1)),1,'first');
    coastIndex2(i)=700+ii(index,1);
    np=find(lattrBall(i,:)>0,1,'last');
    lonBall(i,1:np)=[lontrBall(i,ii(index,2):-1:1) lontrBall(i,np:-1:ii(index,2)+1)];
    latBall(i,1:np)=[lattrBall(i,ii(index,2):-1:1) lattrBall(i,np:-1:ii(index,2)+1)];

    [~,~,ii]=polyxpoly([-4.4 -4.4],[36.6 35.8],lonFall(i,:),latFall(i,:));
    if length(ii)>1
    index=min(ii(:,2));
    lonWagF(i,1:index+1)=[lonFall(i,1:index) -4.4];
    latWagF(i,1:index+1)=[latFall(i,1:index) latFall(i,index)];
    latGate(i,1)=latFall(i,index);
    end%if intersection forward
    
    [~,~,ii]=polyxpoly([-4.4 -4.4],[36.6 35.8],lonBall(i,:),latBall(i,:));
    if length(ii)>1
    index=min(ii(:,2));
    lonWagB(i,1:index+1)=[lonBall(i,1:index) -4.4];
    latWagB(i,1:index+1)=[latBall(i,1:index) latBall(i,index)];
    latGate(i,2)=latBall(i,index);
    end%if intersection backward
    
    if latWagB(i,1)>0 & latWagF(i,1)>0
        npF=find(latWagF(i,:)>0,1,'last');
        npB=find(latWagB(i,:)>0,1,'last');
        np=length(coastIndex1(i):coastIndex2(i))+npB+npF;

        lonWagClosed(i,1:np)=[lonCoast(coastIndex1(i):coastIndex2(i)).' lonWagB(i,1:npB) lonWagF(i,npF:-1:1)];
        latWagClosed(i,1:np)=[latCoast(coastIndex1(i):coastIndex2(i)).' latWagB(i,1:npB) latWagF(i,npF:-1:1)];
        xm=111000*cosd(latWagClosed(i,1:np)).*(lonWagClosed(i,1:np)-xmin*ones(size(lonWagClosed(i,1:np))));
        ym=111000*(latWagClosed(i,1:np)-ymin*ones(size(latWagClosed(i,1:np))));
        ym(ym==0)=nan;
        xm(ym==0)=nan;
        areaWagGate(i)=polyarea(xm,ym);
                
        
        
        np=find(latWagClosed(i,:)>0,1,'last');
    if length(np)>0 & min(latGate(i,:))>0
        yg1=111000*(min(latGate(i,:))-ymin);
        yg2=111000*(max(latGate(i,:))-ymin);
        yG=[yg1:1000:yg2 yg2]; np2=length(yG);
        yGate(i,1:np2)=yG;
        gateFluxSign(i)=sign(latGate(i,2)-latGate(i,1));
    end

    end %if intersections found for both
    
end

lonGate=-4.4*ones(size(latGate));
xGate=111000*cosd(36.3).*(-4.4-xmin)*ones(size(yGate));
xGate(yGate==0)=nan; yGate(yGate==0)=nan;
%% gate flux sigma 27.5 to 28
clear u*
load('uvwNativeGridIsoDepth275.mat','uIso')%native grid
load('isoDepthsNF.mat', 'isoDepth275')
load('isoDepthsNF.mat', 'isoDepth28')
xcm=111000*cosd(YC).*(XC-xmin*ones(size(XC)));
ycm=111000*(YC-ymin*ones(size(YC)));
for i=14:nt
    if yGate(i,1)>0
        np=find(yGate(i,:)>0,1,'last');
        uGate=griddata(xcm,ycm,uIso(:,:,i),xGate(i,1:np),yGate(i,1:np));
        uGate2=0.5.*uGate(2:end)+0.5.*uGate(1:end-1);
        meanUg27(i)=gateFluxSign(i)*mean(uGate);
        
        dGate28=griddata(xcm,ycm,isoDepth28(:,:,i),xGate(i,1:np),yGate(i,1:np));
        dGate282=0.5.*dGate28(2:end)+0.5.*dGate28(1:end-1);
        dGate275=griddata(xcm,ycm,isoDepth275(:,:,i),xGate(i,1:np),yGate(i,1:np));
        dGate2752=0.5.*dGate275(2:end)+0.5.*dGate275(1:end-1);
        meanD27528(i)=mean(dGate28-dGate275);
        
        gateFlux27528(i)=gateFluxSign(i)*sum(abs(diff(yGate(i,1:np))).*uGate2.*(dGate282-dGate2752));
    end
end

%SO- this flux is in m^3/s, calculated each day. the total VOLUME that goes
%is 86400s/day times this thing. DO I WANT THAT?????

figure; plot(gateFlux27528,'o'); title('Gate Flux, \sigma=27 to \sigma=27.5')
%%
clear uIso
save('wagAreaAndFluxGate275v2.mat')

%% figure of total gate flux and its component layers
load('wagAreaAndFluxGateS2.mat','gateFlux*')
load('wagAreaAndFluxGate275v2.mat','gateFlux27528')
load('wagAreaAndFluxGate27v2.mat','gateFlux27275')
load('wagAreaAndFluxGate2675v2.mat','gateFlux267527')
load('wagAreaAndFluxGate265v2.mat','gateFlux2652675')
load('wagAreaAndFluxGate263v2.mat','gateFlux26265')

gateFluxTot=gateFluxS26(1:146)+gateFlux26265(1:146)+gateFlux2652675+gateFlux267527+gateFlux27275+gateFlux27528;
gateFluxTot2=gateFlux26265(1:146)+gateFlux2652675+gateFlux267527+gateFlux27275+gateFluxS26(1:146);

figure; plot(gateFluxS26); hold all
plot(gateFlux26265)
plot(gateFlux2652675)
plot(gateFlux267527)
plot(gateFlux27275)
%plot(gateFlux27528)
plot(gateFluxTot2,'k--','linewidth',2)
%plot(gateFluxTot,'k--','linewidth',2)
%% layer volume budgets- see lagrangeWAGbudget.m

%% precip
load('wagAreaAndFluxGateS2.mat','inWAG')
load('precip.mat')
for i=1:148
   precip(:,:,i)=RAC.*double(inWAG(:,:,i)).*fwflux(:,:,i)/1000;%fwflux is kg/m^2/s, I want m^3/s, suppose 1000kg/m^3
end
precipTot=squeeze(nansum(nansum(precip)));
%% SSH storage term
%load('uvwSSHDailyDepth1rotated148F.mat', 'SSHa')
for i=1:148
   volSSH(:,:,i)=SSHa(:,:,i).*RAC.*double(inWAG(:,:,i)); %m^3
end
storageFlux=diff(squeeze(nansum(nansum(volSSH))))/86400;
%%
load('volumeWagLagrangeFilled.mat')

gateFTot=gateFluxS26;
gateFTot(1:length(gateFlux26265))=gateFTot(1:length(gateFlux26265))+gateFlux26265;
gateFTot(1:length(gateFlux2652675))=gateFTot(1:length(gateFlux2652675))+gateFlux2652675;
gateFTot(1:length(gateFlux267527))=gateFTot(1:length(gateFlux267527))+gateFlux267527;

figure; plot(gateFTot,'k','linewidth',2); 
hold on; plot(15:130,diff(volumeWagLagrangeC(15:131))/86400,'m--','linewidth',2)
hold all; plot(gateFluxS26)
plot(gateFlux26265)
plot(gateFlux2652675)
plot(gateFlux267527)
title('WAG Gate Flux')
xlabel('simulation day')
ylabel('m^3/s')
legend('total','surface to \sigma=26','\sigma=26 to \sigma=26.5','\sigma=26.5 to \sigma=27','\sigma=27 to \sigma=27.5')

figure; plot(gateFluxS26); hold all; plot(diff(volumeWagS26)/86400)
figure; plot(gateFlux26265); hold all; plot(diff(volumeWag26265)/86400)
figure; plot(gateFlux2652675); hold all; plot(diff(volumeWag26527)/86400)
figure; plot(gateFlux267527); hold all; plot(diff(volumeWag27275)/86400)

gateFE=abs(gateFluxS26);
gateFE(1:length(gateFlux26265))=gateFE(1:length(gateFlux26265))+abs(gateFlux26265);
gateFE(1:length(gateFlux2652675))=gateFE(1:length(gateFlux2652675))+abs(gateFlux2652675);
gateFE(1:length(gateFlux267527))=gateFE(1:length(gateFlux267527))+abs(gateFlux267527);

figure; plot(gateFE./gateFE); hold all; plot(abs(gateFluxS26)./gateFE)
plot(abs(gateFlux26265)./gateFE(1:length(gateFlux26265)))
plot(abs(gateFlux2652675)./gateFE(1:length(gateFlux2652675)))
plot(abs(gateFlux267527)./gateFE(1:length(gateFlux267527)))
title('Fraction of WAG Gate Flux Exchange')
xlabel('simulation day')
legend('total','surface to \sigma=26','\sigma=26 to \sigma=26.5','\sigma=26.5 to \sigma=27','\sigma=27 to \sigma=27.5')
axis([1 148 0 1.05])
%% figure of gate widths
load('wagAreaAndFluxGateS.mat','yGate')
yS=yGate;
load('wagAreaAndFluxGate27.mat','yGate')
y27=yGate;
load('wagAreaAndFluxGate265.mat','yGate')
y265=yGate;
load('wagAreaAndFluxGate263.mat','yGate')
y263=yGate;
load('wagAreaAndFluxGate2675.mat','yGate')
y2675=yGate;
load('wagAreaAndFluxGate275.mat','yGate')
y275=yGate;

figure; plot(max(yS,[],2)-min(yS,[],2),'-o')
hold all; plot(max(y263,[],2)-min(y263,[],2),'-o')
plot(max(y265,[],2)-min(y265,[],2),'-o')
plot(max(y2675,[],2)-min(y2675,[],2),'-o')
plot(max(y27,[],2)-min(y27,[],2),'-o')
plot(max(y275,[],2)-min(y275,[],2),'-o')
xlabel('simulation day')
legend('surface','\sigma=26.3','\sigma=26.5','\sigma=26.75','\sigma=27','\sigma=27.5')
title('Gate Width')
ylabel('width in meters')

%% figure of gate depths and velocities
load('wagAreaAndFluxGateS.mat','mean*')
load('wagAreaAndFluxGate275.mat','mean*')
load('wagAreaAndFluxGate27.mat','mean*')
load('wagAreaAndFluxGate2675.mat','mean*')
load('wagAreaAndFluxGate265.mat','mean*')
load('wagAreaAndFluxGate263.mat','mean*')

figure; plot(meanDS,'-o'); hold all
plot(meanD26265,'-o'); plot(meanD2652675,'-o'); plot(meanD267527,'-o'); plot(meanD27275,'-o'); plot(meanD27528,'-o')
xlabel('simulation day')
legend('surface','\sigma=26.3','\sigma=26.5','\sigma=26.75','27','27.5')
title('Gate Depth')
ylabel('depth in meters')

figure; plot(meanUS,'-o'); hold all
plot(meanUg26,'-o'); plot(meanUg265,'-o'); plot(meanUg27,'-o')
xlabel('simulation day')
legend('surface','\sigma=26','\sigma=26.5','\sigma=27')
title('Gate Velocity')
ylabel('m/s')

figure; plot(sign(meanUS),'-o'); hold all
plot(1.1*sign(meanUg26),'-o'); plot(1.2*sign(meanUg265),'-o'); plot(1.3*sign(meanUg27),'-o')
xlabel('simulation day')
legend('surface','\sigma=26','\sigma=26.5','\sigma=27')
title('Flux direction')
ylabel('magnitudes offset for readability')

figure; plot(meanUS.*meanDS.*(max(yS,[],2).'-min(yS,[],2).'));
hold all; plot(gateFluxS)
%% WAG volume over time, estimate

load('diffFluxHorizS26.mat','inWAG2')
inS=inWAG2;
load('diffFluxHoriz27275.mat','inWAG2')
inB=inWAG2;
load('isopycnalDepths.mat')

for i=8:131
    sshVol(i)=nansum(nansum(RAC(inS(:,:,i))));
    bVol(i)=nansum(nansum(RAC(inB(:,:,i)).*iso275depth(inB(:,:,i))));%-iso27depth(inB(:,:,i))));
end
    
figure; plot(diff(sshVol)/(86400)); hold all; plot(diff(bVol)/86400)

%% WAG volume over time, properly
load('geometrySpinupSteady.mat','XC','YC','dInterface','*Coast')
load('distancesAreas.mat','RAC','hFacC')
load('wagAreaAndFluxGate.mat','*WagClosed')
load('isopycnalDepths2.mat','iso*')
load('ssh148.mat','ssh148')
dZ(1,1,:)=diff(dInterface);
cellVol=repmat(dZ,[700 200 1]).*repmat(RAC,[1 1 46]).*hFacC;

volumeWagS26=zeros(148,1);
volumeWag26265=volumeWagS26;
volumeWag26527=volumeWagS26;
volumeWag27275=volumeWagS26;
volumeWagS26c=volumeWagS26;
volumeWag26265c=volumeWagS26;
volumeWag26527c=volumeWagS26;
volumeWag27275c=volumeWagS26;

figure;
for i=9:148
    np=find(latWagClosed(i-8,:)>0,1,'last');
    if np>1
    inNow(:,:,i)=inpolygon(XC,YC,lonWagClosed(i-8,1:np),latWagClosed(i-8,1:np));
    volumeWagS26(i)=squeeze(nansum(nansum((iso26depth(:,:,i)+ssh148(:,:,i)).*inNow(:,:,i).*RAC)));
    [c,h]=contour(XC,YC,double(inNow(:,:,i)),[1 1]);
    if length(c)>0
        npointsi=find(c(2,:)==floor(c(2,:)));
        varhold=max(c(2,npointsi));
        i1=1+find(c(2,:)==varhold);
        iend=varhold+i1-1;
        inWAG2(:,:,i)=inpolygon(XC,YC,c(1,i1:iend),c(2,i1:iend));
        volumeWagS26c(i)=squeeze(nansum(nansum((iso26depth(:,:,i)+ssh148(:,:,i)).*inWAG2(:,:,i).*RAC)));
    end
    end
end

inWagS26=inNow; clear inNow
inWagS26c=inWAG2; clear inWAG2
%%
wagAreaS1=squeeze(nansum(nansum(inWagS26.*repmat(RAC,[1 1 148]))));
wagAreaS2=squeeze(nansum(nansum(inWagS26c.*repmat(RAC,[1 1 148]))));

figure; plot(wagAreaS1);
hold all; plot(wagAreaS2);

%%
load('wagAreaAndFluxGate26.mat','*WagClosed')
for i=1:138
    np=find(latWagClosed(i,:)>0,1,'last');
    if np>1
    inNow(:,:,i)=inpolygon(XC,YC,lonWagClosed(i,1:np),latWagClosed(i,1:np));
    volumeWag26265(i)=squeeze(nansum(nansum((iso265depth(:,:,i)-iso26depth(:,:,i)).*inNow(:,:,i).*RAC)));
    [c,h]=contour(XC,YC,double(inNow(:,:,i)),[1 1]);
    if length(c)>0
        npointsi=find(c(2,:)==floor(c(2,:)));
        varhold=max(c(2,npointsi));
        i1=1+find(c(2,:)==varhold);
        iend=varhold+i1-1;
        inWAG2(:,:,i)=inpolygon(XC,YC,c(1,i1:iend),c(2,i1:iend));
        volumeWag26265c(i)=squeeze(nansum(nansum((-iso26depth(:,:,i)+iso265depth(:,:,i)).*inWAG2(:,:,i).*RAC)));
    end
    end
end
inWag26265=inNow; clear inNow
inWag26265c=inWAG2; clear inWAG2


load('wagAreaAndFluxGate265.mat','*WagClosed')
for i=1:131
    np=find(latWagClosed(i,:)>0,1,'last');
    if np>1
    inNow=inpolygon(XC,YC,lonWagClosed(i,1:np),latWagClosed(i,1:np));
    volumeWag26527(i)=squeeze(nansum(nansum((iso27depth(:,:,i)-iso265depth(:,:,i)).*inNow.*RAC)));
    [c,h]=contour(XC,YC,double(inNow),[1 1]);
    if length(c)>0
        npointsi=find(c(2,:)==floor(c(2,:)));
        varhold=max(c(2,npointsi));
        i1=1+find(c(2,:)==varhold);
        iend=varhold+i1-1;
        inWAG2=inpolygon(XC,YC,c(1,i1:iend),c(2,i1:iend));
        volumeWag26527c(i)=squeeze(nansum(nansum((iso27depth(:,:,i)-iso265depth(:,:,i)).*inWAG2.*RAC)));
    end
    end
end

%%
load('wagAreaAndFluxGate27.mat','*WagClosed')
for i=1:131
    np=find(latWagClosed(i,:)>0,1,'last');
    if np>1
    inNow(:,:,i)=inpolygon(XC,YC,lonWagClosed(i,1:np),latWagClosed(i,1:np));
    volumeWag27275(i)=squeeze(nansum(nansum((iso275depth(:,:,i)-iso27depth(:,:,i)).*inNow(:,:,i).*RAC)));
    [c,h]=contour(XC,YC,double(inNow(:,:,i)),[1 1]);
    if length(c)>0
        npointsi=find(c(2,:)==floor(c(2,:)));
        varhold=max(c(2,npointsi));
        i1=1+find(c(2,:)==varhold);
        iend=varhold+i1-1;
        inWAG2(:,:,i)=inpolygon(XC,YC,c(1,i1:iend),c(2,i1:iend));
        volumeWag27275c(i)=squeeze(nansum(nansum((iso275depth(:,:,i)-iso27depth(:,:,i)).*inWAG2(:,:,i).*RAC)));
    end
    end
end
inWag27275=inNow; clear inNow
inWag27275c=inWAG2; clear inWAG2

for i=50:58
    figure; pcolor(XC,YC,inWag27275(:,:,i)+inWag27275c(:,:,i)); shading 'flat'
    title(num2str(i))
end

%%
volumeWagLagrange=volumeWagS26+volumeWag26265+volumeWag26527+volumeWag27275;
volumeWagLagrangeC=volumeWagS26c+volumeWag26265c+volumeWag26527c+volumeWag27275c;
save('volumeWagLagrange.mat','volume*')
%%
figure; plot(volumeWagLagrange,'k','linewidth',2)
hold all
plot(volumeWagS26); plot(volumeWag26265); plot(volumeWag26527); plot(volumeWag27275)
legend('surface to \sigma=26','\sigma=26 to \sigma=26.5','\sigma=26.5 to \sigma=27','\sigma=27 to \sigma=27.5')
title('Lagrangian WAG Layer Volumes')
xlabel('simulation day')
ylabel('m^3')

figure; plot(volumeWagLagrangeC,'k','linewidth',2)
hold all
plot(volumeWagS26c); plot(volumeWag26265c); plot(volumeWag26527c); plot(volumeWag27275c)
legend('surface to \sigma=26','\sigma=26 to \sigma=26.5','\sigma=26.5 to \sigma=27','\sigma=27 to \sigma=27.5')
title('Lagrangian WAG Layer Volumes via Contours')
xlabel('simulation day')
ylabel('m^3')

figure; plot(volumeWagLagrangeC-volumeWagLagrange,'k','linewidth',2)
hold all
plot(volumeWagS26c-volumeWagS26); plot(volumeWag26265c-volumeWag26265); plot(volumeWag26527c-volumeWag26527); plot(volumeWag27275c-volumeWag27275)
legend('surface to \sigma=26','\sigma=26 to \sigma=26.5','\sigma=26.5 to \sigma=27','\sigma=27 to \sigma=27.5')
title('Lagrangian WAG Layer Volumes, Contours-Originals')
xlabel('simulation day')
ylabel('m^3')

%% fill in missing sections
for i=9:148
   if volumeWagS26(i)==0 & volumeWagS26(i-1)>0
       volumeWagS26(i)=volumeWagS26(i-1);
   end
   if volumeWagS26c(i)==0 & volumeWagS26c(i-1)>0
       volumeWagS26c(i)=volumeWagS26c(i-1);
   end
end
for i=9:138
    if volumeWag26265(i)==0 & volumeWag26265(i-1)>0
       volumeWag26265(i)=volumeWag26265(i-1);
    end
   if volumeWag26265c(i)==0 & volumeWag26265c(i-1)>0
       volumeWag26265c(i)=volumeWag26265c(i-1);
   end
end

for i=15:131
    if volumeWag26527(i)==0 & volumeWag26527(i-1)>0
       volumeWag26527(i)=volumeWag26527(i-1);
    end
   if volumeWag27275(i)==0 & volumeWag27275(i-1)>0
   volumeWag27275(i)=volumeWag27275(i-1);
   end
   
   if volumeWag26527c(i)==0 & volumeWag26527c(i-1)>0
       volumeWag26527c(i)=volumeWag26527c(i-1);
    end
   if volumeWag27275c(i)==0 & volumeWag27275c(i-1)>0
   volumeWag27275c(i)=volumeWag27275c(i-1);
   end
   
end

volumeWagLagrange=volumeWagS26+volumeWag26265+volumeWag26527+volumeWag27275;
volumeWagLagrangeC=volumeWagS26c+volumeWag26265c+volumeWag26527c+volumeWag27275c;
save('volumeWagLagrangeFilled.mat','volume*')
%%
figure; plot(diff(volumeWagLagrange),'k','linewidth',2)
hold all
plot(diff(volumeWagS26)); plot(diff(volumeWag26265)); plot(diff(volumeWag26527)); plot(diff(volumeWag27275))
legend('total','surface to \sigma=26','\sigma=26 to \sigma=26.5','\sigma=26.5 to \sigma=27','\sigma=27 to \sigma=27.5')
title('Changes in Lagrangian WAG Layer Volumes')
xlabel('simulation day')
ylabel('m^3/day')

figure; plot(diff(volumeWagLagrangeC),'k','linewidth',2)
hold all
plot(diff(volumeWagS26c)); plot(diff(volumeWag26265c)); plot(diff(volumeWag26527c)); plot(diff(volumeWag27275c))
legend('total','surface to \sigma=26','\sigma=26 to \sigma=26.5','\sigma=26.5 to \sigma=27','\sigma=27 to \sigma=27.5')
title('Changes in Lagrangian WAG Contoured Layer Volumes')
xlabel('simulation day')
ylabel('m^3/day')

figure; plot(diff(volumeWagLagrangeC)-diff(volumeWagLagrange),'k','linewidth',2)
hold all
plot(diff(volumeWagS26c)-diff(volumeWagS26)); plot(diff(volumeWag26265c)-diff(volumeWag26265)); plot(diff(volumeWag26527c)-diff(volumeWag26527)); plot(diff(volumeWag27275c)-diff(volumeWag27275))
legend('total','surface to \sigma=26','\sigma=26 to \sigma=26.5','\sigma=26.5 to \sigma=27','\sigma=27 to \sigma=27.5')
title('Differences in Changes in Lagrangian WAG Contoured vs Not Layer Volumes')
xlabel('simulation day')
ylabel('m^3/day')

%% 'fill' volume errors with NaNs
load('volumeWagLagrange.mat')

volumeWagS26(volumeWagS26==0)=NaN;
volumeWagS26c(volumeWagS26c==0)=NaN;
volumeWag26265(volumeWag26265==0)=NaN;
volumeWag26265c(volumeWag26265c==0)=NaN;
volumeWag26527(volumeWag26527==0)=NaN;
volumeWag26527c(volumeWag26527c==0)=NaN;
volumeWag27275(volumeWag27275==0)=NaN;
volumeWag27275c(volumeWag27275c==0)=NaN;

volumeWagLagrange=volumeWagS26+volumeWag26265+volumeWag26527+volumeWag27275;
volumeWagLagrangeC=volumeWagS26c+volumeWag26265c+volumeWag26527c+volumeWag27275c;

save('volumeWagLagrangeNan.mat')

load('wagAreaAndFluxGate27.mat','gateFlux27275')
load('wagAreaAndFluxGate265.mat','gateFlux26527')
load('wagAreaAndFluxGate26.mat','gateFlux26265')
load('wagAreaAndFluxGate.mat','gateFluxS26')

figure; plot(gateFluxS26); hold all; plot(diff(volumeWagS26)/86400)
figure; plot(gateFlux26265); hold all; plot(diff(volumeWag26265)/86400)
figure; plot(gateFlux2652675); hold all; plot(diff(volumeWag26527)/86400)
figure; plot(gateFlux267527); hold all; plot(diff(volumeWag27275)/86400)

%% make dVolume/dT
load('volumeWagLagrangeFilled.mat')%Filled or Nan

volumeFluxS26=diff(volumeWagS26)/86400;
volumeFlux26265=diff(volumeWag26265)/86400;
volumeFlux26527=diff(volumeWag26527)/86400;
volumeFlux27275=diff(volumeWag27275)/86400;

volumeFTot=diff(volumeWagS26+volumeWag26265+volumeWag26527+volumeWag27275)/86400;

load('wagAreaAndFluxGate27.mat','gateFlux27275')
load('wagAreaAndFluxGate265.mat','gateFlux26527')
load('wagAreaAndFluxGate26.mat','gateFlux26265')
load('wagAreaAndFluxGate.mat','gateFluxS26')
load('volumesfluxesWAGeuler2017b.mat','vfluxs')
vfluxs2=[0.5*vfluxs(2:2:end)+0.5*vfluxs(1:2:end-1), vfluxs(end)].';
vfluxs2=vfluxs2(1:131);

%gateFluxS26(gateFluxS26==0)=NaN;
%gateFlux26265(gateFlux26265==0)=NaN;
%gateFlux26527(gateFlux26527==0)=NaN;
%gateFlux27275(gateFlux27275==0)=NaN;

gateFtot=gateFluxS26(1:131)+gateFlux26265(1:131)+gateFlux2652675(1:131)+gateFlux267527(1:131);

%0=-dV/dt +gateFlux +(precip-evap) flux +cross-isopycnal flux 
% = -volumeFlux +gateFlux +vfluxs2
%crossisopycnal flux 26= dV/dt -gateFlux -(precip-evap) flux
%                    = volumeFlux -gateFlux -vfluxs2
%0=-dV/dt +gateFlux +(precip-evap) flux +crossisopycnal flux bottom -crossisopycnal flux top
%crossisopycnal flux bottom 265= dV/dt -gateFlux -(precip-evap) flux +crossisopycnal flux top
cross275a=-vfluxs2+volumeFTot(1:131)-gateFtot(1:131).';
figure; plot(gateFtot,'k','linewidth',2)
hold on; plot(-volumeFTot,'m','linewidth',2)
plot(vfluxs2,'b','linewidth',2)%need to redo this for moving boundary
plot(cross275a,'c')
title('total')

cross26a=-vfluxs2+volumeFluxS26(1:131)-gateFluxS26(1:131).';
figure; plot(gateFluxS26,'k','linewidth',2)
hold on; plot(-volumeFluxS26,'m','linewidth',2)
plot(vfluxs2,'b','linewidth',2)
plot(cross26a,'c')
title('S26')

cross26265=volumeFlux26265(1:131)-gateFlux26265(1:131).';
figure; plot(gateFlux26265,'k','linewidth',2)
hold on; plot(-volumeFlux26265,'m','linewidth',2)
plot(cross26265,'c')
title('26265')

cross26527=volumeFlux26527(1:131)-gateFlux2652675(1:131).';
figure; plot(gateFlux2652675,'k','linewidth',2)
hold on; plot(-volumeFlux26527,'m','linewidth',2)
plot(cross26265,'c')
title('26527')

cross27275=volumeFlux27275(1:131)-gateFlux267527(1:131).';
figure; plot(gateFlux267527,'k','linewidth',2)
hold on; plot(-volumeFlux27275,'m','linewidth',2)
plot(cross27275,'c')
title('27275')

%cross26265=cross265-cross26a, so cross265=cross26265+cross26a
cross265=cross26265+cross26a;
cross27=cross26527+cross265;
cross275=cross27275+cross27;

figure; plot(vfluxs2,'b','linewidth',2); hold all
plot(cross26a); plot(cross265); plot(cross27); plot(cross275)
plot(cross275a,'k','linewidth',2)
legend('precip','surface layer cross','265','27','275','275 from total volume')
%% plot sample boundaries
% load('wagAreaAndFluxGate26.mat')
% for i=85:88
% figure; pcolor(XC,YC,inWag3(:,:,2,i)+0.5*edgeCell(:,:,2,i)); shading 'flat'; hold on
% plot(lonCoast,latCoast,'k','linewidth',2);
% plot(lontrFall(i,:),lattrFall(i,:),'r')
% plot(lontrBall(i,:),lattrBall(i,:),'b')
% np=find(latWagClosed(i,:)>0,1,'last');
% plot(lonWagClosed(i,1:np),latWagClosed(i,1:np),'g--','linewidth',2)
% title(strcat('\sigma=26 Manifolds and WAG edge, day ',num2str(i)))
% xlabel('longitude'); ylabel('latitude')
% legend('cell type','coast','unstable manifold','stable manifold','WAG edge')
% %colormap(cbrewer('div','PRGn',3))
% axis([-5.5 -1 35 37])
% end

%% surface plot boundaries sample
% load('wagAreaAndFluxGate.mat')
% for i=85:88
% figure; pcolor(XC,YC,inWag3(:,:,1,i)+0.5*edgeCell(:,:,1,i)); shading 'flat'; hold on
% plot(lonCoast,latCoast,'k','linewidth',2);
% plot(lontrFall(i-8,:),lattrFall(i-8,:),'r')
% plot(lontrBall(i-8,:),lattrBall(i-8,:),'b')
% np=find(latWagClosed(i-8,:)>0,1,'last');
% plot(lonWagClosed(i-8,1:np),latWagClosed(i-8,1:np),'g--','linewidth',2)
% title(strcat('Surface Manifolds and WAG edge, day ',num2str(i)))
% xlabel('longitude'); ylabel('latitude')
% legend('cell type','coast','unstable manifold','stable manifold','WAG edge')
% %colormap(cbrewer('div','PRGn',3))
% axis([-5.5 -1 35 37])
% end

%% possible error bars on volume in Lagrange WAG
% load('lagrangeWAGboundary.mat','inW*')
% load('volumeWagLagrangeFilled.mat')
% k=1;
% inWag3(:,:,k,:)=inWAGS(:,:,1:131);
% 
% edgeCell=false([700 200 16 131]);
% 
% edgeCell(2:end-1,2:end-1,2:end-1,:)=(inWag3(2:end-1,2:end-1,2:end-1,:)==1)&((inWag3(1:end-2,2:end-1,2:15,:)==0)|(inWag3(3:end,2:end-1,2:15,:)==0)...
%     |(inWag3(2:end-1,1:end-2,2:15,:)==0)|(inWag3(2:end-1,3:end,2:15,:)==0)|(inWag3(2:end-1,2:end-1,3:16,:)==0)|(inWag3(2:end-1,2:end-1,1:14,:)==0));
% 
% edgeCell(2:end-1,2:end-1,1,:)=(inWag3(2:end-1,2:end-1,1,:)==1)&((inWag3(1:end-2,2:end-1,1,:)==0)|(inWag3(3:end,2:end-1,1,:)==0)...
%     |(inWag3(2:end-1,1:end-2,1,:)==0)|(inWag3(2:end-1,3:end,1,:)==0)|(inWag3(2:end-1,2:end-1,2,:)==0));
% 
% load('distancesAreas.mat','RAC','hFacC')
% load('geometrySpinupSteady.mat','dInterface')
% clear dZ
% dZ(1,1,1:16)=diff(dInterface(1:17));
% cellVol=repmat(RAC,[1 1 16]).*hFacC(:,:,1:16).*repmat(dZ,[700 200 1]);
% 
% for i=1:131
%     volErr(i)=sum(cellVol(edgeCell(:,:,:,i)));
%     volErr2(i)=0.1*sum(cellVol(edgeCell(:,:,:,i)));
% end
% volDifErr(1:130)=sqrt(volErr2(1:130).^2+volErr2(2:131).^2);
% 
% volTot=volumeWagS26+volumeWag26265+volumeWag26527+volumeWag27275;
% 
% figure; plot(15:131,volTot(15:131),'k','linewidth',2)
% hold on; plot(15:131,volTot(15:131)+volErr2(15:131).','b--')
% plot(15:131,volTot(15:131)-volErr2(15:131).','b--')
% 
% figure; plot(15:130,diff(volTot(15:131))./86400,'k','linewidth',2)
% hold on; plot(15:130,diff(volTot(15:131))./86400+volDifErr(15:130).'./86400,'b--')
% hold on; plot(15:130,diff(volTot(15:131))./86400-volDifErr(15:130).'./86400,'b--')
% plot(gateFtot,'r','linewidth',2)
% NB='volErr2 uses 1/10 edge cell volumes';
% save('volumeWagLagrangeErr2.mat','vol*')