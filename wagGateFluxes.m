runthis1=false;
if runthis1
%old! load('wagAreaAndFluxGate2675v2filled.mat','xGate*','yGate*','gateFluxSign*')
load('geometrySpinupSteady','*C*')
xmin=min(min(XC));
ymin=min(min(YC));
xcm=111111*cosd(YC).*(XC-xmin*ones(size(XC)));
ycm=111111*(YC-ymin*ones(size(YC)));
%% recalculating xGate,yGate
%load('manifoldsIso263Int10NF2hourlyGate')
fillneeded2=[];
for i=120:168%hourly, daily was 14:146
    igate=find(lonWagClosed(i,:)==-4.4);
    if length(igate)==2
        lonGate(i,:)=lonWagClosed(i,igate);
        latGate(i,:)=latWagClosed(i,igate);
    else
        disp(num2str(i))
%        fillneeded2=[fillneeded2 i];
%    elseif ismember(i,[143 145])
%         igate=95:96;%surface
%         lonGate(i,:)=lonWagClosed(i,igate);
%         latGate(i,:)=latWagClosed(i,igate);     
     end
end
%
lonGate(latGate==0)=nan;
latGate(latGate==0)=nan;
yG=111111*(latGate-ymin*ones(size(latGate)));
clear yGate
for i=120:168%1:146
    if yG(i,2)<yG(i,1)
        holdvar=[ yG(i,2):1000:yG(i,1) yG(i,1)];
        gateFluxSign(i)=1;
    else
        holdvar=[yG(i,1):1000:yG(i,2) yG(i,2)];
        gateFluxSign(i)=-1;
    end
    yGate(i,1:length(holdvar))=holdvar; 
    
end
yGate(yGate==0)=nan;
latG=yGate./111111 +ymin*ones(size(yGate));
lonG=-4.4*ones(size(latG)); lonG(isnan(latG))=nan;
xGate=111111*cosd(latG).*(lonG-xmin*ones(size(lonG)));

% for ii=1:length(fillneeded2)%2 for sigma27.5
%     i=fillneeded(ii);
%     [~,n]=size(xGate27);
%     xGate(i,1:n)=xGate27(i,:);
%     yGate(i,1:n)=yGate27(i,:);
%     gateFluxSign(i)=gateFluxSign27(i);
% end

%xGate263=xGate; yGate263=yGate; gateFluxSign263=gateFluxSign;
save('surfaceGateHourly.mat','xGate','yGate','lonGate','latGate','*G','gateFluxSign')%sigma263

end %if runthis1
%% 
%% integrate gate fluxes better using triscatteredinterp
%S263
disp('start S')
%xGate=xGateS; yGate=yGateS; gateFluxSign=gateFluxSignS;
load('isoDepthsNF')
load('manifoldsSurfaceInt14NFanalyzed')
load('surfaceGate1')
%load('isoDepthsNF2.mat');
%load('manifoldsSurfInt9hourlyGate.mat');
%load('surfaceGateHourly')

%load('advTeastdailyNF.mat')
%load vorta from vorticity.mat
load('ueastdailyNF')
%load('uHourly15to30rotatedNF2.mat')
Urot=Urot(:,:,1:30,1:146);
%Urot=Urot.*vorta; clear vorta

load('distancesAreas','dyg','hFacW')
DYG=reshape(dyg,[700 200]);
load('geometrySpinupSteady','dInterface');
dZ(1,1,1:30)=diff(dInterface(1:31));
wface=hFacW(:,:,1:30).*repmat(dZ,[700 200 1]).*repmat(DYG,[1 1 30]);
%
%AdvEt=AdvEt./repmat(wface,[1 1 1 146]);
%
hours=(15*24):24*30;%for velocity
days=hours/24;
%%
%indicesKeep=121:169;%days 20 and 21, last one is day 22
%Urot2=cat(3,Urot(:,:,1,indicesKeep),Urot(:,:,:,indicesKeep));
Urot2=cat(3,Urot(:,:,1,:),Urot);
clear Urot AdvEs AdvEt
Urot2=Urot2(250:450,:,:,:);
 load('uvwSSHDailyDepth1rotated148F.mat','SSHa')
 ssh148=SSHa; clear SSHa
%load('sshThourly.mat', 'SSH')
%SSH=SSH(:,:,360:720);
%ssh148=SSH(:,:,indicesKeep); clear SSH
zu(1)=0.3;
zu(2:31)=-0.5*(dInterface(1:30)+dInterface(2:31));

xu3=repmat(xcm(250:450,:),[1 1 31]);
yu3=repmat(ycm(250:450,:),[1 1 31]);
holdvar(1,1,1:31)=zu;
zu3=repmat(holdvar,[201 200 1]); clear holdvar
%%
for i=14:146 %1:48
    i
    i2=i;%119+i;
    if yGate(i2,1)>0
        np=find(yGate(i2,:)>0,1,'last')
    end
    if np>1
        %
        %uGate=griddata(xcm,ycm,Urot(:,:,i+1),xGate(i,1:np),yGate(i,1:np));
        
        %surface to sigma=26.3
        %dGate=griddata(xcm,ycm,isoDepth263(:,:,i2+1),xGate(i2,1:np),yGate(i2,1:np));
        %dGate=max(dGate./2,5*ones(size(dGate)));
        
        %surface to sigma=27.5
        dGate=griddata(xcm,ycm,isoDepth275(:,:,i2+1),xGate(i2,1:np),yGate(i2,1:np));
        dGate=max(dGate./2,5*ones(size(dGate)));
        
        %all
        hGate=griddata(xcm,ycm,ssh148(:,:,i+1),xGate(i2,1:np),yGate(i2,1:np));
        
        zv=[max(hGate) 0:-1:min(-dGate) min(-dGate)];
        polyz=[-dGate hGate(end:-1:1) -dGate(1)];
        polyy=[yGate(i2,1:np) yGate(i2,np:-1:1) yGate(i2,1)];
        [xv,zv2]=meshgrid(xGate(i2,1:np),zv);
        [yv,~]=meshgrid(yGate(i2,1:np),zv);
        interpolateU=TriScatteredInterp(xu3(:),yu3(:),zu3(:),reshape(Urot2(:,:,:,i+1),[],1));
        uv=interpolateU(xv,yv,zv2);
        ingate=inpolygon(yv,zv2,polyy,polyz);
        uv(~ingate)=0;
        gateFluxS(i+1)=gateFluxSign(i2)*abs(trapz(zv,trapz(yGate(i2,np:-1:1),uv,2)));

    end
end
mean(gateFluxS(15:146))
mean(abs(gateFluxS(15:146)))
Note='gate volume flux using surface gate and ssh+sigma=27.5 height';
disp('saving')
save('gateAdvectionBarotropic1.mat','*lux*','xGate','yGate','gateFluxSign')
%save('gateAdvectionHourlySurface.mat','*lux*','xGate','yGate','gateFluxSign')

%%
%
disp('start 263')
load('sigma263GateHourly.mat')
% load('isoDepthsNF.mat');
% load('ueastdailyNF.mat','Urot')
% Urot2=cat(3,Urot(:,:,1,:),Urot);
% clear Urot
% Urot2=Urot2(250:450,:,:,:);
for i=1:48
    i
    i2=119+i;
    if yGate(i2,1)>0
        np=find(yGate(i2,:)>0,1,'last')
    end
    if np>1
        %uGate=griddata(xcm,ycm,Urot(:,:,i+1),xGate(i,1:np),yGate(i,1:np));
        dGate1=griddata(xcm,ycm,isoDepth265(:,:,i2+1),xGate(i2,1:np),yGate(i2,1:np));
        hGate1=griddata(xcm,ycm,isoDepth263(:,:,i2+1),xGate(i2,1:np),yGate(i2,1:np));
        hGate1(isnan(hGate1))=5;
        dGate1(isnan(dGate1))=5;
        hGate=max(hGate1./2,5*ones(size(hGate1)));
        dGate=max(0.5*(dGate1+hGate1),5*ones(size(hGate)));
        
        zv=[max(-hGate):-1:min(-dGate) min(-dGate)];
        polyz=[-dGate -hGate(end:-1:1) -dGate(1)];
        polyy=[yGate(i2,1:np) yGate(i2,np:-1:1) yGate(i2,1)];
        [xv,zv2]=meshgrid(xGate(i2,1:np),zv);
        [yv,~]=meshgrid(yGate(i2,1:np),zv);
        interpolateU=TriScatteredInterp(xu3(:),yu3(:),zu3(:),reshape(Urot2(:,:,:,i+1),[],1));
        uv=interpolateU(xv,yv,zv2);
        ingate=inpolygon(yv,zv2,polyy,polyz);
        uv(~ingate)=0;
        gateFlux263(i+1)=gateFluxSign(i2)*abs(trapz(zv,trapz(yGate(i2,np:-1:1),uv,2)));

    end
end
save('gateAdvectionSigma263hourly2.mat','*lux*','xGate','yGate','gateFluxSign')
%%
disp('start 265')
load('sigma265Gate1.mat')
%load('sigma265GateHourly.mat')
% load('isoDepthsNF.mat');
% load('ueastdailyNF.mat','Urot')
% Urot2=cat(3,Urot(:,:,1,:),Urot);
% clear Urot
% Urot2=Urot2(250:450,:,:,:);
for i=1:146%48
    i
    i2=i;
    np=0;
    %i2=i+119;
    if yGate(i2,1)>0
        np=find(yGate(i2,:)>0,1,'last')
    end
    if np>1
        %uGate=griddata(xcm,ycm,Urot(:,:,i+1),xGate(i,1:np),yGate(i,1:np));
        dGate1=griddata(xcm,ycm,isoDepth2675(:,:,i2+1),xGate(i2,1:np),yGate(i2,1:np));
        hGate1=griddata(xcm,ycm,isoDepth265(:,:,i2+1),xGate(i2,1:np),yGate(i2,1:np));
        hGate2=griddata(xcm,ycm,isoDepth263(:,:,i2+1),xGate(i2,1:np),yGate(i2,1:np));
        hGate1(isnan(hGate1))=5;
        dGate1(isnan(dGate1))=5;
        hGate2(isnan(hGate2))=5;
        hGate=max(0.5*(hGate1+hGate2),5*ones(size(hGate1)));
        dGate=max(0.5*(hGate1+dGate1),hGate);
        
        zv=[max(-hGate):-1:min(-dGate) min(-dGate)];
        polyz=[-dGate -hGate(end:-1:1) -dGate(1)];
        polyy=[yGate(i2,1:np) yGate(i2,np:-1:1) yGate(i2,1)];
        [xv,zv2]=meshgrid(xGate(i2,1:np),zv);
        [yv,~]=meshgrid(yGate(i2,1:np),zv);
        interpolateU=TriScatteredInterp(xu3(:),yu3(:),zu3(:),reshape(Urot2(:,:,:,i+1),[],1));
        uv=interpolateU(xv,yv,zv2);
        ingate=inpolygon(yv,zv2,polyy,polyz);
        uv(~ingate)=0;
        gateFlux265(i+1)=gateFluxSign(i2)*abs(trapz(zv,trapz(yGate(i2,np:-1:1),uv,2)));

    end
end
%%
figure; plot(gateFlux265);
save('gateAdvection265april.mat','*lux*','xGate','yGate','gateFluxSign')
%save('gateAdvectionSigma265hourly2.mat','*lux*','xGate','yGate','gateFluxSign')
%
disp('start 2675')
load('sigma2675Gate1.mat')
%load('sigma2675GateHourly.mat')
% load('isoDepthsNF.mat');
% load('ueastdailyNF.mat','Urot')
% Urot2=cat(3,Urot(:,:,1,:),Urot);
% clear Urot
% Urot2=Urot2(250:450,:,:,:);
for i=1:145%48
    i
    i2=i;
    %i2=119+i;
    np=0;
    if yGate(i2,1)>0
        np=find(yGate(i2,:)>0,1,'last')
    end
    if np>1
        %uGate=griddata(xcm,ycm,Urot(:,:,i+1),xGate(i,1:np),yGate(i,1:np));
        dGate1=griddata(xcm,ycm,isoDepth27(:,:,i2+1),xGate(i2,1:np),yGate(i2,1:np));
        hGate1=griddata(xcm,ycm,isoDepth2675(:,:,i2+1),xGate(i2,1:np),yGate(i2,1:np));
        hGate2=griddata(xcm,ycm,isoDepth265(:,:,i2+1),xGate(i2,1:np),yGate(i2,1:np));
        hGate=max(0.5*(hGate1+hGate2),5*ones(size(hGate1)));
        dGate=max(0.5*(hGate1+dGate1),hGate);
        
        zv=[max(-hGate):-1:min(-dGate) min(-dGate)];
        polyz=[-dGate -hGate(end:-1:1) -dGate(1)];
        polyy=[yGate(i2,1:np) yGate(i2,np:-1:1) yGate(i2,1)];
        [xv,zv2]=meshgrid(xGate(i2,1:np),zv);
        [yv,~]=meshgrid(yGate(i2,1:np),zv);
        interpolateU=TriScatteredInterp(xu3(:),yu3(:),zu3(:),reshape(Urot2(:,:,:,i+1),[],1));
        uv=interpolateU(xv,yv,zv2);
        ingate=inpolygon(yv,zv2,polyy,polyz);
        uv(~ingate)=0;
        gateFlux2675(i+1)=gateFluxSign(i2)*abs(trapz(zv,trapz(yGate(i2,np:-1:1),uv,2)));

    end
end
save('gateAdvectionSigma2675april.mat','*lux*','xGate','yGate','gateFluxSign')
%save('gateAdvectionSigma2675hourly2.mat','*lux*','xGate','yGate','gateFluxSign')
%%
disp('start 27')
load('sigma27Gate1.mat')
%load('sigma27GateHourly.mat')

% load('isoDepthsNF.mat');
% load('ueastdailyNF.mat','Urot')
% Urot2=cat(3,Urot(:,:,1,:),Urot);
% clear Urot
% Urot2=Urot2(250:450,:,:,:);
for i=1:145%48
    i
    i2=i;
    %i2=119+i;
    np=0;
    if yGate(i2,1)>0
        np=find(yGate(i2,:)>0,1,'last')
    end
    if np>1
        %uGate=griddata(xcm,ycm,Urot(:,:,i+1),xGate(i,1:np),yGate(i,1:np));
        dGate1=griddata(xcm,ycm,isoDepth275(:,:,i2+1),xGate(i2,1:np),yGate(i2,1:np));
        hGate1=griddata(xcm,ycm,isoDepth27(:,:,i2+1),xGate(i2,1:np),yGate(i2,1:np));
        hGate2=griddata(xcm,ycm,isoDepth2675(:,:,i2+1),xGate(i2,1:np),yGate(i2,1:np));
        hGate=max(0.5*(hGate1+hGate2),5*ones(size(hGate1)));
        dGate=max(0.5*(hGate1+dGate1),hGate);
        
        zv=[max(-hGate):-1:min(-dGate) min(-dGate)];
        polyz=[-dGate -hGate(end:-1:1) -dGate(1)];
        polyy=[yGate(i2,1:np) yGate(i2,np:-1:1) yGate(i2,1)];
        [xv,zv2]=meshgrid(xGate(i2,1:np),zv);
        [yv,~]=meshgrid(yGate(i2,1:np),zv);
        interpolateU=TriScatteredInterp(xu3(:),yu3(:),zu3(:),reshape(Urot2(:,:,:,i+1),[],1));
        uv=interpolateU(xv,yv,zv2);
        ingate=inpolygon(yv,zv2,polyy,polyz);
        uv(~ingate)=0;
        gateFlux27(i+1)=gateFluxSign(i2)*abs(trapz(zv,trapz(yGate(i2,np:-1:1),uv,2)));

    end
end
save('gateAdvectionSigma27april.mat','*lux*','xGate','yGate','gateFluxSign')
%save('gateAdvectionSigma27hourly2.mat','*lux*','xGate','yGate','gateFluxSign')

%%
disp('start 275')
load('sigma275GateHourly.mat')
% load('isoDepthsNF.mat');
% load('ueastdailyNF.mat','Urot')
% Urot2=cat(3,Urot(:,:,1,:),Urot);
% clear Urot
% Urot2=Urot2(250:450,:,:,:);
for i=1:48
    i
    i2=119+i;
    if yGate(i2,1)>0
        np=find(yGate(i2,:)>0,1,'last')
    end
    if np>1
        %uGate=griddata(xcm,ycm,Urot(:,:,i+1),xGate(i,1:np),yGate(i,1:np));
        dGate1=griddata(xcm,ycm,isoDepth275(:,:,i2+1),xGate(i2,1:np),yGate(i2,1:np));
        hGate1=griddata(xcm,ycm,isoDepth27(:,:,i2+1),xGate(i2,1:np),yGate(i2,1:np));
        %hGate2=griddata(xcm,ycm,isoDepth2675(:,:,i+1),xGate(i,1:np),yGate(i,1:np));
        hGate=max(0.5*(hGate1+dGate1),5*ones(size(hGate1)));
        dGate=max(dGate1,hGate);
        
        zv=[max(-hGate):-1:min(-dGate) min(-dGate)];
        polyz=[-dGate -hGate(end:-1:1) -dGate(1)];
        polyy=[yGate(i2,1:np) yGate(i2,np:-1:1) yGate(i2,1)];
        [xv,zv2]=meshgrid(xGate(i2,1:np),zv);
        [yv,~]=meshgrid(yGate(i2,1:np),zv);
        interpolateU=TriScatteredInterp(xu3(:),yu3(:),zu3(:),reshape(Urot2(:,:,:,i+1),[],1));
        uv=interpolateU(xv,yv,zv2);
        ingate=inpolygon(yv,zv2,polyy,polyz);
        uv(~ingate)=0;
        gateFlux275(i+1)=gateFluxSign(i2)*abs(trapz(zv,trapz(yGate(i2,np:-1:1),uv,2)));

    end
end
save('gateAdvectionSigma275hourly.mat','*lux*','xGate','yGate','gateFluxSign')
%%
load('gateAdvectionSigma27indexedMar.mat','gateFlux263265')
gateFlux27275=gateFlux263265;
load('gateAdvectionSigma2675indexedMar.mat','gateFlux263265')
gateFlux267527=gateFlux263265;
load('gateAdvectionSigma265indexedMar.mat','gateFlux263265')
gateFlux2652675=gateFlux263265;
load('gateAdvectionSigma263indexedMar.mat','gateFlux263265')
load('gateAdvectionS263indexedMar.mat','gateFluxS263')
gateFluxTotMar=gateFluxS263+gateFlux263265+gateFlux2652675+gateFlux267527+gateFlux27275;
figure; plot(gateFluxTotFeb); hold all; plot(gateFluxTotFeb2); plot(gateFluxTotMar,'--','linewidth',2)
save('gateAdvectionIndexedMar.mat','gateFlux*')
%%
runthis=false;
if runthis
%% Feb set below
load('uvwSSHDailyDepth1rotated148F','Urot')
load('isoDepthsNF.mat','isoDepth263')
load('uvwSSHDailyDepth1rotated148F.mat','SSHa')
ssh148=SSHa; clear SSHa
load('tsDepth1native.mat')
%% 
%update 12/20- indexing for velocity fields here and in creating manifolds
%does not match, manifolds indexed from 0 while all saved fields index from
%1, so I need to re-index manifolds to get correct fluxes
%also will have to update contents
%load('wagAreaAndFluxGate2675v2filled.mat','gateFluxSign*')
%% surface to \sigma=26.3 gate fluxes
disp('start S')
xGate=xGateS; yGate=yGateS; gateFluxSign=gateFluxSignS;

for i=14:146
    if yGateS(i,1)>0
        np=find(yGateS(i,:)>0,1,'last');
    end
%     if ismember(i,[38:42 51 52 80 81 88 109 110 145])
%         np=find(yGate263(i,:)>0,1,'last');
%         xGate(i,1:np)=xGate263(i,1:np);
%         yGate(i,1:np)=yGate263(i,1:np);
%         gateFluxSign(i)=gateFluxSign263(i);
%     elseif ismember(i,[82:85 87 106:108])
%         np=find(yGate265(i,:)>0,1,'last');
%         xGate(i,1:np)=xGate265(i,1:np);
%         yGate(i,1:np)=yGate265(i,1:np);
%         gateFluxSign(i)=gateFluxSign265(i);
%     end
    if np>1
        uGate=griddata(xcm,ycm,Urot(:,:,i+1),xGate(i,1:np),yGate(i,1:np));
        dGate=griddata(xcm,ycm,isoDepth263(:,:,i+1),xGate(i,1:np),yGate(i,1:np));
        hGate=griddata(xcm,ycm,ssh148(:,:,i+1),xGate(i,1:np),yGate(i,1:np));
        uGate2=0.5.*uGate(2:end)+0.5.*uGate(1:end-1);
        meanUS(i+1)=mean(uGate)*gateFluxSign(i);
        dGate2=0.5.*dGate(2:end)+0.5.*dGate(1:end-1);
        hGate2=0.5.*hGate(2:end)+0.5.*hGate(1:end-1);
        meanDS(i+1)=mean(dGate+hGate);
        gateFluxS263(i+1)=gateFluxSign(i)*sum(abs(diff(yGate(i,1:np))).*uGate2.*(hGate2+dGate2));
        
        gateFluxBin1(i+1)=gateFluxSign(i)*sum(abs(diff(yGate(i,1:np))).*uGate2.*5);
        
        sGateS=griddata(xcm,ycm,S(:,:,i+1),xGate(i,1:np),yGate(i,1:np));
        tGateS=griddata(xcm,ycm,T(:,:,i+1),xGate(i,1:np),yGate(i,1:np));
        sGate2=0.5.*sGateS(2:end)+0.5.*sGateS(1:end-1);
        tGate2=0.5.*tGateS(2:end)+0.5.*tGateS(1:end-1);
        
        %salinity
        gateSfluxS263(i+1)=gateFluxSign(i)*sum(abs(diff(yGate(i,1:np))).*uGate2.*sGate2.*(hGate2+dGate2));
        %temp
        gateTfluxS263(i+1)=gateFluxSign(i)*sum(abs(diff(yGate(i,1:np))).*uGate2.*tGate2.*(hGate2+dGate2));
    end
end


save('gateAdvectionS263indexedFeb.mat','*lux*','xGate','yGate','gateFluxSign','meanUS','meanDS')

%%
load('geometrySpinupSteady','Angle*')
load('uvwNativeGridIsoDepth263.mat', 'uIso','vIso')
nt=162;
Urot=uIso.*repmat(AngleCS,[1 1 nt]) - vIso.*repmat(AngleSN,[1 1 nt]); clear vIso uIso
load('isoDepthsNF.mat')
load('tsNativeGridIsoDepth263.mat', 'sIso','tIso')
%%
xGate=xGate263; yGate=yGate263; gateFluxSign=gateFluxSign263;
disp('start 263')
for i=50:146
    if yGate(i,1)>0
        np=find(yGate(i,:)>0,1,'last');
    end
%     if ismember(i,[60:63 71])
%         np=find(yGateS(i,:)>0,1,'last');
%         xGate(i,1:np)=xGateS(i,1:np);
%         yGate(i,1:np)=yGateS(i,1:np);
%         gateFluxSign(i)=gateFluxSignS(i);
%     elseif ismember(i,[82:87 106:108])
%         np=find(yGate265(i,:)>0,1,'last');
%         xGate(i,1:np)=xGate265(i,1:np);
%         yGate(i,1:np)=yGate265(i,1:np);
%         gateFluxSign(i)=gateFluxSign265(i);
%     end
    if np>1
        uGate=griddata(xcm,ycm,Urot(:,:,i+1),xGate(i,1:np),yGate(i,1:np));
        uGate2=0.5.*uGate(2:end)+0.5.*uGate(1:end-1);
         meanUg263(i)=gateFluxSign(i)*mean(uGate);
        
        dGate26=griddata(xcm,ycm,isoDepth263(:,:,i+1),xGate(i,1:np),yGate(i,1:np));
        dGate262=0.5.*dGate26(2:end)+0.5.*dGate26(1:end-1);
        dGate265=griddata(xcm,ycm,isoDepth265(:,:,i+1),xGate(i,1:np),yGate(i,1:np));
        dGate2652=0.5.*dGate265(2:end)+0.5.*dGate265(1:end-1);
        meanD263265(i+1)=mean(dGate265-dGate26);
        
        gateFlux263265(i+1)=gateFluxSign(i)*sum(abs(diff(yGate(i,1:np))).*uGate2.*(dGate2652-dGate262));
        
        hGate=griddata(xcm,ycm,ssh148(:,:,i+1),xGate(i,1:np),yGate(i,1:np));
        hGate2=0.5.*hGate(2:end)+0.5.*hGate(1:end-1);
        gateFlux263S(i+1)=gateFluxSign(i)*sum(abs(diff(yGate(i,1:np))).*uGate2.*(hGate2+dGate262));
        
        sGateS=griddata(xcm,ycm,sIso(:,:,i+1),xGate(i,1:np),yGate(i,1:np));
        tGateS=griddata(xcm,ycm,tIso(:,:,i+1),xGate(i,1:np),yGate(i,1:np));
        sGate2=0.5.*sGateS(2:end)+0.5.*sGateS(1:end-1);
        tGate2=0.5.*tGateS(2:end)+0.5.*tGateS(1:end-1);
        
        %salinity
        gateSflux263265(i+1)=gateFluxSign(i)*sum(abs(diff(yGate(i,1:np))).*uGate2.*sGate2.*(dGate2652-dGate262));
        gateSflux263S(i+1)=gateFluxSign(i)*sum(abs(diff(yGate(i,1:np))).*uGate2.*sGate2.*(hGate2+dGate262));
        %temp
        gateTflux263265(i+1)=gateFluxSign(i)*sum(abs(diff(yGate(i,1:np))).*uGate2.*tGate2.*(dGate2652-dGate262));
        gateTflux263S(i+1)=gateFluxSign(i)*sum(abs(diff(yGate(i,1:np))).*uGate2.*sGate2.*(hGate2+dGate262));
    end
end

save('gateAdvection263265indexedFeb.mat','*lux*','xGate','yGate','gateFluxSign','meanU*','meanD*')
%%
load('uvwNativeGridIsoDepth265.mat', 'uIso','vIso')
nt=162;
Urot=uIso.*repmat(AngleCS,[1 1 nt]) - vIso.*repmat(AngleSN,[1 1 nt]); clear vIso uIso
load('isoDepthsNF.mat','isoDepth2675')
load('tsNativeGridIsoDepth265.mat', 'sIso','tIso')
%%
xGate=xGate265; yGate=yGate265; gateFluxSign=gateFluxSign265;
disp('start 265')
for i=14:146
    if yGate(i,1)>0
        np=find(yGate(i,:)>0,1,'last');
    end
%     if ismember(i,[60:63 71])
%         np=find(yGateS(i,:)>0,1,'last');
%         xGate(i,1:np)=xGateS(i,1:np);
%         yGate(i,1:np)=yGateS(i,1:np);
%         gateFluxSign(i)=gateFluxSignS(i);
%     elseif ismember(i,[64 66:70 103])
%         np=find(yGate263(i,:)>0,1,'last');
%         xGate(i,1:np)=xGate263(i,1:np);
%         yGate(i,1:np)=yGate263(i,1:np);
%         gateFluxSign(i)=gateFluxSign263(i);
%     end 
    if np>1
        uGate=griddata(xcm,ycm,Urot(:,:,i+1),xGate(i,1:np),yGate(i,1:np));
        uGate2=0.5.*uGate(2:end)+0.5.*uGate(1:end-1);
         meanUg263(i)=gateFluxSign(i)*mean(uGate);
        
        dGate26=griddata(xcm,ycm,isoDepth265(:,:,i+1),xGate(i,1:np),yGate(i,1:np));
        dGate262=0.5.*dGate26(2:end)+0.5.*dGate26(1:end-1);
        dGate265=griddata(xcm,ycm,isoDepth2675(:,:,i+1),xGate(i,1:np),yGate(i,1:np));
        dGate2652=0.5.*dGate265(2:end)+0.5.*dGate265(1:end-1);
        meanD2652675(i+1)=mean(dGate265-dGate26);
        
        gateFlux2652675(i+1)=gateFluxSign(i)*sum(abs(diff(yGate(i,1:np))).*uGate2.*(dGate2652-dGate262));
        
        %code here!
        hGate=griddata(xcm,ycm,isoDepth263(:,:,i+1),xGate(i,1:np),yGate(i,1:np));
        hGate2=0.5.*hGate(2:end)+0.5.*hGate(1:end-1);
        gateFlux265263(i+1)=gateFluxSign(i)*sum(abs(diff(yGate(i,1:np))).*uGate2.*(dGate262-hGate2));
 
        sGateS=griddata(xcm,ycm,sIso(:,:,i+1),xGate(i,1:np),yGate(i,1:np));
        tGateS=griddata(xcm,ycm,tIso(:,:,i+1),xGate(i,1:np),yGate(i,1:np));
        sGate2=0.5.*sGateS(2:end)+0.5.*sGateS(1:end-1);
        tGate2=0.5.*tGateS(2:end)+0.5.*tGateS(1:end-1);
        
        %salinity
        gateSflux2652675(i+1)=gateFluxSign(i)*sum(abs(diff(yGate(i,1:np))).*uGate2.*sGate2.*(dGate2652-dGate262));
        gateSflux265263(i+1)=gateFluxSign(i)*sum(abs(diff(yGate(i,1:np))).*uGate2.*sGate2.*(dGate262-hGate2));
        %temp
        gateTflux2652675(i+1)=gateFluxSign(i)*sum(abs(diff(yGate(i,1:np))).*uGate2.*tGate2.*(dGate2652-dGate262));
        gateTflux265263(i+1)=gateFluxSign(i)*sum(abs(diff(yGate(i,1:np))).*uGate2.*tGate2.*(dGate262-hGate2));
    end
end

save('gateAdvection2652675indexedFeb.mat','*lux*','xGate','yGate','gateFluxSign','meanU*','meanD*')

%%
load('uvwNativeGridIsoDepth2675.mat', 'uIso','vIso')
nt=162;
Urot=uIso.*repmat(AngleCS,[1 1 nt]) - vIso.*repmat(AngleSN,[1 1 nt]); clear vIso uIso
load('isoDepthsNF.mat','isoDepth27')
load('tsNativeGridIsoDepth2675.mat', 'sIso','tIso')
%%
xGate=xGate2675; yGate=yGate2675; gateFluxSign=gateFluxSign2675;
disp('start 2675')
for i=14:146
    if yGate(i,1)>0
        np=find(yGate(i,:)>0,1,'last');
    end
%     if ismember(i,[62 63])
%         np=find(yGateS(i,:)>0,1,'last');
%         xGate(i,1:np)=xGateS(i,1:np);
%         yGate(i,1:np)=yGateS(i,1:np);
%         gateFluxSign(i)=gateFluxSignS(i);
%     elseif ismember(i,67:70)
%         np=find(yGate263(i,:)>0,1,'last');
%         xGate(i,1:np)=xGate263(i,1:np);
%         yGate(i,1:np)=yGate263(i,1:np);
%         gateFluxSign(i)=gateFluxSign263(i);
%     elseif ismember(i,[15 43 48 51:56 72 112])
%         np=find(yGate265(i,:)>0,1,'last');
%         xGate(i,1:np)=xGate265(i,1:np);
%         yGate(i,1:np)=yGate265(i,1:np);
%         gateFluxSign(i)=gateFluxSign265(i);
%     elseif ismember(i,[61 64])
%         np=find(yGate27(i,:)>0,1,'last');
%         xGate(i,1:np)=xGate27(i,1:np);
%         yGate(i,1:np)=yGate27(i,1:np);
%         gateFluxSign(i)=gateFluxSign27(i);
%     end 
    if np>1
        uGate=griddata(xcm,ycm,Urot(:,:,i),xGate(i,1:np),yGate(i,1:np));
        uGate2=0.5.*uGate(2:end)+0.5.*uGate(1:end-1);
         meanUg2675(i+1)=gateFluxSign(i)*mean(uGate);
        
        dGate26=griddata(xcm,ycm,isoDepth2675(:,:,i+1),xGate(i,1:np),yGate(i,1:np));
        dGate262=0.5.*dGate26(2:end)+0.5.*dGate26(1:end-1);
        dGate265=griddata(xcm,ycm,isoDepth27(:,:,i+1),xGate(i,1:np),yGate(i,1:np));
        dGate2652=0.5.*dGate265(2:end)+0.5.*dGate265(1:end-1);
        meanD267527(i+1)=mean(dGate265-dGate26);
        
        gateFlux267527(i+1)=gateFluxSign(i)*sum(abs(diff(yGate(i,1:np))).*uGate2.*(dGate2652-dGate262));
        
        hGate=griddata(xcm,ycm,isoDepth265(:,:,i+1),xGate(i,1:np),yGate(i,1:np));
        hGate2=0.5.*hGate(2:end)+0.5.*hGate(1:end-1);
        gateFlux2675265(i+1)=gateFluxSign(i)*sum(abs(diff(yGate(i,1:np))).*uGate2.*(dGate262-hGate2));
        
        sGateS=griddata(xcm,ycm,sIso(:,:,i+1),xGate(i,1:np),yGate(i,1:np));
        tGateS=griddata(xcm,ycm,tIso(:,:,i+1),xGate(i,1:np),yGate(i,1:np));
        sGate2=0.5.*sGateS(2:end)+0.5.*sGateS(1:end-1);
        tGate2=0.5.*tGateS(2:end)+0.5.*tGateS(1:end-1);
        
        %salinity
        gateSflux267527(i+1)=gateFluxSign(i)*sum(abs(diff(yGate(i,1:np))).*uGate2.*sGate2.*(dGate2652-dGate262));
        gateSflux2675265(i+1)=gateFluxSign(i)*sum(abs(diff(yGate(i,1:np))).*uGate2.*sGate2.*(dGate262-hGate2));
        %temp
        gateTflux267527(i+1)=gateFluxSign(i)*sum(abs(diff(yGate(i,1:np))).*uGate2.*tGate2.*(dGate2652-dGate262));
        gateTflux2675265(i+1)=gateFluxSign(i)*sum(abs(diff(yGate(i,1:np))).*uGate2.*tGate2.*(dGate262-hGate2));
    end
end

save('gateAdvection267527indexedFeb.mat','*lux*','xGate','yGate','gateFluxSign','meanU*','meanD*')
%%
load('uvwNativeGridIsoDepth27.mat', 'uIso','vIso')
nt=162;
Urot=uIso.*repmat(AngleCS,[1 1 nt]) - vIso.*repmat(AngleSN,[1 1 nt]); clear vIso uIso
load('isoDepthsNF.mat','isoDepth275')
load('tsNativeGridIsoDepth27.mat', 'sIso','tIso')
%%
xGate=xGate27; yGate=yGate27; gateFluxSign=gateFluxSign27;
disp('start 27')
for i=14:146
    if yGate(i,1)>0
        np=find(yGate(i,:)>0,1,'last');
    end
%     if ismember(i,[62 63])
%         np=find(yGateS(i,:)>0,1,'last');
%         xGate(i,1:np)=xGateS(i,1:np);
%         yGate(i,1:np)=yGateS(i,1:np);
%         gateFluxSign(i)=gateFluxSignS(i);
%     elseif ismember(i,67:70)
%         np=find(yGate263(i,:)>0,1,'last');
%         xGate(i,1:np)=xGate263(i,1:np);
%         yGate(i,1:np)=yGate263(i,1:np);
%         gateFluxSign(i)=gateFluxSign263(i);
%     elseif ismember(i,[71 117])
%         np=find(yGate2675(i,:)>0,1,'last');
%         xGate(i,1:np)=xGate2675(i,1:np);
%         yGate(i,1:np)=yGate2675(i,1:np);
%         gateFluxSign(i)=gateFluxSign2675(i);
%     end 
    if np>1
        uGate=griddata(xcm,ycm,Urot(:,:,i),xGate(i,1:np),yGate(i,1:np));
        uGate2=0.5.*uGate(2:end)+0.5.*uGate(1:end-1);
         meanUg27275(i+1)=gateFluxSign(i)*mean(uGate);
        
        dGate26=griddata(xcm,ycm,isoDepth27(:,:,i+1),xGate(i,1:np),yGate(i,1:np));
        dGate262=0.5.*dGate26(2:end)+0.5.*dGate26(1:end-1);
        dGate265=griddata(xcm,ycm,isoDepth275(:,:,i+1),xGate(i,1:np),yGate(i,1:np));
        dGate2652=0.5.*dGate265(2:end)+0.5.*dGate265(1:end-1);
        meanD27275(i+1)=mean(dGate265-dGate26);
        
        gateFlux27275(i+1)=gateFluxSign(i)*sum(abs(diff(yGate(i,1:np))).*uGate2.*(dGate2652-dGate262));
        
        hGate=griddata(xcm,ycm,isoDepth2675(:,:,i+1),xGate(i,1:np),yGate(i,1:np));
        hGate2=0.5.*hGate(2:end)+0.5.*hGate(1:end-1);
        gateFlux272675(i+1)=gateFluxSign(i)*sum(abs(diff(yGate(i,1:np))).*uGate2.*(dGate262-hGate2));
        
        sGateS=griddata(xcm,ycm,sIso(:,:,i+1),xGate(i,1:np),yGate(i,1:np));
        tGateS=griddata(xcm,ycm,tIso(:,:,i+1),xGate(i,1:np),yGate(i,1:np));
        sGate2=0.5.*sGateS(2:end)+0.5.*sGateS(1:end-1);
        tGate2=0.5.*tGateS(2:end)+0.5.*tGateS(1:end-1);
        
        %salinity
        gateSflux27275(i+1)=gateFluxSign(i)*sum(abs(diff(yGate(i,1:np))).*uGate2.*sGate2.*(dGate2652-dGate262));
        gateSflux272675(i+1)=gateFluxSign(i)*sum(abs(diff(yGate(i,1:np))).*uGate2.*sGate2.*(dGate262-hGate2));
        %temp
        gateTflux27275(i+1)=gateFluxSign(i)*sum(abs(diff(yGate(i,1:np))).*uGate2.*tGate2.*(dGate2652-dGate262));
        gateTflux272675(i+1)=gateFluxSign(i)*sum(abs(diff(yGate(i,1:np))).*uGate2.*tGate2.*(dGate262-hGate2));
    end
end

save('gateAdvection27275indexedFeb.mat','*lux*','xGate','yGate','gateFluxSign','meanU*','meanD*')

%%
load('uvwNativeGridIsoDepth275.mat', 'uIso','vIso')
nt=148;
Urot=uIso.*repmat(AngleCS,[1 1 nt]) - vIso.*repmat(AngleSN,[1 1 nt]); clear vIso uIso
load('isoDepthsNF.mat','isoDepth275')
load('tsNativeGridIsoDepth275.mat', 'sIso','tIso')
%%
xGate=xGate275; yGate=yGate275; gateFluxSign=gateFluxSign275;
disp('start 275')
for i=14:146
    if yGate(i,1)>0
        np=find(yGate(i,:)>0,1,'last');
    end

%     if ismember(i,67:70)
%         np=find(yGate263(i,:)>0,1,'last');
%         xGate(i,1:np)=xGate263(i,1:np);
%         yGate(i,1:np)=yGate263(i,1:np);
%         gateFluxSign(i)=gateFluxSign263(i);
%     elseif ismember(i,[71 117])
%         np=find(yGate2675(i,:)>0,1,'last');
%         xGate(i,1:np)=xGate2675(i,1:np);
%         yGate(i,1:np)=yGate2675(i,1:np);
%         gateFluxSign(i)=gateFluxSign2675(i);
%     elseif ismember(i,[39 56 66 81 118 121 146])
%         np=find(yGate27(i,:)>0,1,'last');
%         xGate(i,1:np)=xGate27(i,1:np);
%         yGate(i,1:np)=yGate27(i,1:np);
%         gateFluxSign(i)=gateFluxSign27(i);
%     end 
    if np>1
        uGate=griddata(xcm,ycm,Urot(:,:,i+1),xGate(i,1:np),yGate(i,1:np));
        uGate2=0.5.*uGate(2:end)+0.5.*uGate(1:end-1);
         meanUg275(i+1)=gateFluxSign(i)*mean(uGate);
        
        dGate26=griddata(xcm,ycm,isoDepth27(:,:,i+1),xGate(i,1:np),yGate(i,1:np));
        dGate262=0.5.*dGate26(2:end)+0.5.*dGate26(1:end-1);
        dGate265=griddata(xcm,ycm,isoDepth275(:,:,i+1),xGate(i,1:np),yGate(i,1:np));
        dGate2652=0.5.*dGate265(2:end)+0.5.*dGate265(1:end-1);
        meanD27527(i+1)=mean(dGate265-dGate26);
        
        gateFlux27527(i+1)=gateFluxSign(i)*sum(abs(diff(yGate(i,1:np))).*uGate2.*(dGate2652-dGate262));
        
                
        sGateS=griddata(xcm,ycm,sIso(:,:,i+1),xGate(i,1:np),yGate(i,1:np));
        tGateS=griddata(xcm,ycm,tIso(:,:,i+1),xGate(i,1:np),yGate(i,1:np));
        sGate2=0.5.*sGateS(2:end)+0.5.*sGateS(1:end-1);
        tGate2=0.5.*tGateS(2:end)+0.5.*tGateS(1:end-1);
        
        %salinity
        gateSflux27527(i+1)=gateFluxSign(i)*sum(abs(diff(yGate(i,1:np))).*uGate2.*sGate2.*(dGate2652-dGate262));
        %temp
        gateTflux27527(i+1)=gateFluxSign(i)*sum(abs(diff(yGate(i,1:np))).*uGate2.*tGate2.*(dGate2652-dGate262));
    end
        
end

save('gateAdvection275indexedFeb.mat','*lux*','xGate','yGate','gateFluxSign','meanU*','meanD*')

disp('done gate fluxes')
end
%% cross-isopycnal t,s
runthis=false;
if runthis
load('vertVolIso275filled.mat', 'vertVol')
load('tsNativeGridIsoDepth275.mat', 'sIso')
load('tsNativeGridIsoDepth275.mat', 'tIso')
load('lagrangeWAGboundaryPortionFilled.mat', 'inWAG275f', 'inWAG27f')
%%
crossIsoTb=squeeze(nansum(nansum(vertVol.*tIso(:,:,1:146).*inWAG275f./5)));
crossIsoSb=squeeze(nansum(nansum(vertVol.*sIso(:,:,1:146).*inWAG275f./5)));

crossIsoTt=squeeze(nansum(nansum(vertVol.*tIso(:,:,1:146).*inWAG27f./5)));
crossIsoSt=squeeze(nansum(nansum(vertVol.*sIso(:,:,1:146).*inWAG27f./5)));

crossIsoErrSTB=abs(crossIsoSb-crossIsoSt)/2;
crossIsoS=0.5*(crossIsoSb+crossIsoSt);
crossIsoErrTTB=abs(crossIsoTb-crossIsoTt)/2;
crossIsoT=0.5*(crossIsoTb+crossIsoTt);
%%
load('gateAdvection275.mat','gate*flux*')
load('gateAdvection27275.mat','gate*flux*')
load('gateAdvection267527.mat','gate*flux*')
load('gateAdvection2652675.mat','gate*flux*')
load('gateAdvection263265.mat','gate*flux*')
load('gateAdvectionS263.mat','gate*flux*')
%%
gateFluxS1=0.5*gateSfluxS263+0.5*gateSflux263S;
gateErrS1=abs(gateSfluxS263-gateSflux263S)/2;
gateFluxS2=0.5*gateSflux263265+0.5*gateSflux265263;
gateErrS2=abs(gateSflux263265-gateSflux265263)/2;
gateFluxS3=0.5*gateSflux2652675+0.5*gateSflux2675265;
gateErrS3=abs(gateSflux2652675-gateSflux2675265)/2;
gateFluxS4=0.5*gateSflux267527+0.5*gateSflux272675;
gateErrS4=abs(gateSflux267527-gateSflux272675)/2;
gateFluxS5=0.5*gateSflux27275+0.5*gateSflux27527;
gateErrS5=abs(gateSflux27275-gateSflux27527)/2;
gateFluxST=gateSfluxS263+gateSflux263265+gateSflux2652675+gateSflux267527+gateSflux27275;
gateFluxSB=gateSflux263S+gateSflux265263+gateSflux2675265+gateSflux272675+gateSflux27527;
gateErrSTB=abs(gateFluxST-gateFluxSB)/2;
gateErrSSum=gateErrS1+gateErrS2+gateErrS3+gateErrS4+gateErrS5;
%%
gateFluxT1=0.5*gateTfluxS263+0.5*gateTflux263S;
gateErrT1=abs(gateTfluxS263-gateTflux263S)/2;
gateFluxT2=0.5*gateTflux263265+0.5*gateTflux265263;
gateErrT2=abs(gateTflux263265-gateTflux265263)/2;
gateFluxT3=0.5*gateTflux2652675+0.5*gateTflux2675265;
gateErrT3=abs(gateTflux2652675-gateTflux2675265)/2;
gateFluxT4=0.5*gateTflux267527+0.5*gateTflux272675;
gateErrT4=abs(gateTflux267527-gateTflux272675)/2;
gateFluxT5=0.5*gateTflux27275+0.5*gateTflux27527;
gateErrT5=abs(gateTflux27275-gateTflux27527)/2;
gateFluxTT=gateTfluxS263+gateTflux263265+gateTflux2652675+gateTflux267527+gateTflux27275;
gateFluxTB=gateTflux263S+gateTflux265263+gateTflux2675265+gateTflux272675+gateTflux27527;
gateErrTTB=abs(gateFluxTT-gateFluxTB)/2;
gateErrTSum=gateErrT1+gateErrT2+gateErrT3+gateErrT4+gateErrT5;
%%
load('geometrySpinupSteady','dInterface')
load('distancesAreas','RAC','hFacC')
dZ(1,1,1:46)=diff(dInterface);
cellVol=repmat(RAC,[1 1 46]).*hFacC.*repmat(dZ,[700 200 1]);
nt=146;
%cellVol4=repmat(cellVol,[1 1 1 nt]);
%cellVol4=cellVol4(:,:,1:20,:);

load('lagrangeWAGboundaryPortionFilledTB.mat', 'inWag3portionF*')
load('uvwSSHDailyDepth1rotated148F.mat', 'SSHa')
load('tsPress162NF.mat','PotTave')
PotTave=PotTave(:,:,1:20,1:146);
%potTcontentTave=squeeze(sum(sum(sum(PotTave.*cellVol4.*inWag3portionFT./5))));
%potTcontentBave=squeeze(sum(sum(sum(PotTave.*cellVol4.*inWag3portionFB./5)))); 
PotTave=squeeze(PotTave(:,:,1,:));
potTsshTave=squeeze(sum(sum(PotTave.*repmat(RAC,[1 1 nt]).*squeeze(inWag3portionFT(:,:,1,:))./5)));
potTsshBave=squeeze(sum(sum(PotTave.*repmat(RAC,[1 1 nt]).*squeeze(inWag3portionFB(:,:,1,:))./5)));
clear PotTave
load('tsSigmaSnapshots162NF.mat', 'PotT')
PotT=PotT(:,:,1:20,1:146);
%potTcontentB=squeeze(sum(sum(sum(PotT.*cellVol4.*inWag3portionFB./5))));
%potTcontentT=squeeze(sum(sum(sum(PotT.*cellVol4.*inWag3portionFT./5)))); 
PotT=squeeze(PotT(:,:,1,:));
potTsshT=squeeze(sum(sum(PotT.*repmat(RAC,[1 1 nt]).*squeeze(inWag3portionFT(:,:,1,:))./5)));
potTsshB=squeeze(sum(sum(PotT.*repmat(RAC,[1 1 nt]).*squeeze(inWag3portionFB(:,:,1,:))./5)));
clear PotT 
contentT1=[potTcontentB potTcontentBave potTcontentT potTcontentTave];
contentTssh=[potTsshB potTsshBave potTsshT potTsshTave];
dTdt1=diff(contentT1,1,1)./86400;
dTdt1(13,:)=0;
dTdt=diff(contentT1+contentTssh,1,1)./86400;
dTdt(13,:)=0;
figure; plot(1:145,dTdt1); hold all; plot(1:145,dTdt,'--')
%%
load('tsPress162NF.mat','PractSave')
PractSave=PractSave(:,:,1:20,1:146);
%ScontentTave=squeeze(sum(sum(sum(PractSave.*cellVol4.*inWag3portionFT./5))));
%ScontentBave=squeeze(sum(sum(sum(PractSave.*cellVol4.*inWag3portionFB./5)))); 
PractSave=squeeze(PractSave(:,:,1,:));
SsshTave=squeeze(sum(sum((PractSave.*repmat(RAC,[1 1 nt]).*squeeze(inWag3portionFT(:,:,1,:))./5))));
SsshBave=squeeze(sum(sum((PractSave.*repmat(RAC,[1 1 nt]).*squeeze(inWag3portionFB(:,:,1,:))./5)))); 
clear PractSave
load('tsSigmaSnapshots162NF.mat', 'PractS')
PractS=PractS(:,:,1:20,1:146);
%ScontentB=squeeze(sum(sum(sum(PractS.*cellVol4.*inWag3portionFB./5))));
%ScontentT=squeeze(sum(sum(sum(PractS.*cellVol4.*inWag3portionFT./5)))); 
PractS=squeeze(PractS(:,:,1,:));
SsshT=squeeze(sum(sum((PractS.*repmat(RAC,[1 1 nt]).*squeeze(inWag3portionFT(:,:,1,:))./5))));
SsshB=squeeze(sum(sum((PractS.*repmat(RAC,[1 1 nt]).*squeeze(inWag3portionFB(:,:,1,:))./5)))); 
clear PractS
contentS1=[ScontentB ScontentBave ScontentT ScontentTave];
contentSssh=[SsshB SsshBave SsshT SsshTave];
dSdt1=diff(contentS1,1,1)./86400;
dSdt=diff(contentS1+contentSssh,1,1)./86400;
dSdt(13,:)=0;
dSdt1(13,:)=0;
clear inWag3portionF*
%%
load('lagrangeWAGadvectionTS.mat')
load('temperatureLagrangeWAGdif.mat', 'difTZopen','difTHopenS')
load('tsLagrangeWAGcontent.mat')
%%
figure; plot(gateSfluxS263+gateSflux263265+gateSflux2652675+gateSflux267527+gateSflux27275); hold all; plot(crossIsoS)
plot(gateSfluxS263+gateSflux263265+gateSflux2652675+gateSflux267527+gateSflux27275+crossIsoS.')
legend('gate','crossIso','sum')
xlabel('days'); ylabel('psu m^3/s')
title('Salinity advective fluxes')
%%
dTdtM=mean(dTdt,2);
errdTdtU=max(dTdt,[],2)-dTdtM;
errdTdtL=dTdtM-min(dTdt,[],2);
errTotTU=errdTdtL+gateErrTSum(1:145).'+crossIsoErrTTB(1:145);
errTotTL=errdTdtU+gateErrTSum(1:145).'+crossIsoErrTTB(1:145);
%%
figure; plot(gateFluxTB); hold all; plot(crossIsoT)
plot(difTZopen); plot(difTHopenS)
plot(-dTdtM);
plot(gateFluxTB(1:145).'+crossIsoT(1:145)+difTZopen(1:145)+difTHopenS(1:145)-dTdtM,'k--','linewidth',2)
errorbar(1:145,gateFluxTB(1:145).'+crossIsoT(1:145)+difTZopen(1:145)+difTHopenS(1:145)-dTdtM,errTotTL,errTotTU,'k.')
legend('gate','crossIso','diffusion through bottom','diffusion through side','-dT/dt','sum')
xlabel('days'); ylabel('\circ C m^3/s')
%%
figure; plot(gateFluxTB(1:145).'+crossIsoT(1:145)-dTdtM)
hold all;
plot(difTZopen); plot(difTHopenS)
plot(gateFluxTB(1:145).'+crossIsoT(1:145)+difTZopen(1:145)+difTHopenS(1:145)-dTdtM,'k--','linewidth',2)
errorbar(1:145,gateFluxTB(1:145).'+crossIsoT(1:145)+difTZopen(1:145)+difTHopenS(1:145)-dTdtM,errTotTL,errTotTU,'k.')
legend('gate+crossIso-dT/dt','diffusion through bottom','diffusion through side','sum')
xlabel('days'); ylabel('\circ C m^3/s')
title('Temperature fluxes')
%%
figure; 
plot(gateFluxTB(1:145).'+crossIsoT(1:145)+difTZopen(1:145)+difTHopenS(1:145),'linewidth',2)
hold all;
errorbar(1:145,gateFluxTB(1:145).'+crossIsoT(1:145)+difTZopen(1:145)+difTHopenS(1:145),gateErrTSum(1:145).'+crossIsoErrTTB(1:145),'b.')
plot(1:145,dTdt)

end