%% load for circulation
load('edgesWAGeuler2017NF.mat', 'inWag','XC','YC')
load('uvwSSHDailyDepth1rotated148F.mat','Urot','Vrot','SSHa')

load('distancesAreas')
load('geometrySpinupSteady')
xmin=min(XC(:)); ymin=min(YC(:));
xm=111111*cosd(YC).*(XC-xmin*ones(size(XC)));
ym=111111*(YC-ymin*ones(size(YC)));
xcoast=111111*cosd(latCoast).*(lonCoast-xmin*ones(size(lonCoast)));
ycoast=111111*(latCoast-ymin*ones(size(latCoast)));
xum=111111*cosd(YU).*(XU-xmin*ones(size(XU)));
yum=111111*(YU-ymin*ones(size(YU)));
xvm=111111*cosd(YV).*(XV-xmin*ones(size(XV)));
yvm=111111*(YV-ymin*ones(size(YV)));
%% square edge
xi1=375;%300
xi2=425;%350
yi1=100;%60
yi2=150;%100
x1=xvm(xi1,yi1:yi2);%left bottom to top
y1=yvm(xi1,yi1:yi2);
x2=xum(xi1:xi2,yi2);%top left to right
y2=yum(xi1:xi2,yi2);
x3=xvm(xi2,yi1:yi2);%right bottom to top
y3=yvm(xi2,yi1:yi2);
x4=xum(xi1:xi2,yi1);%bottom left to right
y4=yum(xi1:xi2,yi1);
x2=x2(:); y2=y2(:); x1=x1(:); y1=y1(:); 
x4=x4(:); y4=y4(:); x3=x3(:); y3=y3(:); 
xSquare=[x4(:); x3(:);x2(end:-1:1);x1(end:-1:1)];
ySquare=[y4(:); y3(:);y2(end:-1:1);y1(end:-1:1)];
ulogic2=false([700 200]); ulogic4=ulogic2; ulogic4(xi1:xi2,yi1)=true; ulogic2(xi1:xi2,yi2)=true;
vlogic1=false([700 200]); vlogic3=vlogic1; vlogic1(xi1,yi1:yi2)=true; vlogic3(xi2,yi1:yi2)=true;
DYU=reshape(dyu,[700 200]);
DXV=reshape(dxv,[700 200]);
DYC=reshape(dyc,[700 200]);
DXC=reshape(dxc,[700 200]);
%earthR=6371.0*1e3;
ds1=DYC(xi1,yi1:yi2);%sqrt(diff([x4(1); x1(:)]).^2 +diff([y4(1); y1(:)]).^2);%
ds2=DXC(xi1:xi2,yi2);%sqrt(diff([x1(end); x2(:)]).^2 +diff([y1(end);y2(:)]).^2);%
ds3=DYC(xi2,yi1:yi2);%sqrt(diff([x3(:); x2(end)]).^2 +diff([y3(:); y2(end)]).^2);%
ds4=DXC(xi1:xi2,yi1);%sqrt(diff([x4(:); x3(1)]).^2 +diff([y4(:); y3(1)]).^2);%
ds1=ds1(:);ds2=ds2(:);ds3=ds3(:);ds4=ds4(:);
%[~,~,nt]=size(Uext);
dx=diff([xSquare(:); xSquare(1)]); dy=diff([ySquare(:); ySquare(1)]);
vortLogic=false([700 200]); vortLogic(xi1:xi2-1,yi1:yi2-1)=true;
%% vorticity budget
xSquare=xSquare.';
ySquare=ySquare.';
[ integral1, sign1 ] = circulationInt( xum,yum,Urot,xvm,yvm,Vrot,xSquare,ySquare );
%%
load('momentumWindForcing.mat','*ext')
%rotate to east-north
load('geometrySpinupSteady','Angle*')
[~,~,nt]=size(Uext);
windE=Uext.*repmat(AngleCS,[1 1 nt]) - Vext.*repmat(AngleSN,[1 1 nt]);  
windN=Uext.*repmat(AngleSN,[1 1 nt]) + Vext.*repmat(AngleCS,[1 1 nt]); 
[ wind, sign2 ] = circulationInt( xum,yum,windE,xvm,yvm,windN,xSquare,ySquare );
[ windc, ~ ] = circulationInt3( xum,yum,windE,xvm,yvm,windN,xSquare,ySquare );

load('momentumDiffusionSurface.mat')
%remove nans (where there is land)
UDif2a(isnan(UDif2a))=0;
VDif2a(isnan(VDif2a))=0;
UDif2b(isnan(UDif2b))=0;
VDif2b(isnan(VDif2b))=0;
%rotate to east-north
diffE=UDif2a.*repmat(AngleCS,[1 1 nt]) - VDif2a.*repmat(AngleSN,[1 1 nt]);  
diffN=UDif2a.*repmat(AngleSN,[1 1 nt]) + VDif2a.*repmat(AngleCS,[1 1 nt]); 
[ diffusion, sign3 ] = circulationInt( xum,yum,diffE,xvm,yvm,diffN,xSquare,ySquare );
[ diffusionc, sign3 ] = circulationInt3( xum,yum,diffE,xvm,yvm,diffN,xSquare,ySquare );

diffE2=UDif2b.*repmat(AngleCS,[1 1 nt]) - VDif2b.*repmat(AngleSN,[1 1 nt]);  
diffN2=UDif2b.*repmat(AngleSN,[1 1 nt]) + VDif2b.*repmat(AngleCS,[1 1 nt]); 
[ diffusion2, sign30 ] = circulationInt( xum,yum,diffE2,xvm,yvm,diffN2,xSquare,ySquare );
[ diffusion2c, sign30 ] = circulationInt3( xum,yum,diffE2,xvm,yvm,diffN2,xSquare,ySquare );

load('momentumAdvDissSurface.mat','*Diss')
%rotate to east-north
dissE=UDiss.*repmat(AngleCS,[1 1 nt]) - VDiss.*repmat(AngleSN,[1 1 nt]);  
dissN=UDiss.*repmat(AngleSN,[1 1 nt]) + VDiss.*repmat(AngleCS,[1 1 nt]); 
[ dissipation, sign4 ] = circulationInt( xum,yum,dissE,xvm,yvm,dissN,xSquare,ySquare );
[ dissipationc, ~ ] = circulationInt3( xum,yum,dissE,xvm,yvm,dissN,xSquare,ySquare );

load('momentumTendAbSurface.mat','Ab*')
%rotate to east-north
abE=AbU.*repmat(AngleCS,[1 1 nt]) - AbV.*repmat(AngleSN,[1 1 nt]);  
abN=AbU.*repmat(AngleSN,[1 1 nt]) + AbV.*repmat(AngleCS,[1 1 nt]); 
[ timestep, sign5 ] = circulationInt( xum,yum,abE,xvm,yvm,abN,xSquare,ySquare );
[ timestepc, ~ ] = circulationInt3( xum,yum,abE,xvm,yvm,abN,xSquare,ySquare );

load('momentumAdvDissSurface.mat','Adv*')
%rotate to east-north
advE=AdvU.*repmat(AngleCS,[1 1 nt]) - AdvV.*repmat(AngleSN,[1 1 nt]);  
advN=AdvU.*repmat(AngleSN,[1 1 nt]) + AdvV.*repmat(AngleCS,[1 1 nt]); 
[ advection, sign6 ] = circulationInt( xum,yum,advE,xvm,yvm,advN,xSquare,ySquare );

[ advectionc, ~ ] = circulationInt3( xum,yum,advE,xvm,yvm,advN,xSquare,ySquare );

load('momentumTendAbSurface.mat','*tend')
nt=148;
%rotate to east-north
tendE=Utend.*repmat(AngleCS,[1 1 nt])./86400 - Vtend.*repmat(AngleSN,[1 1 nt])./86400;  
tendN=Utend.*repmat(AngleSN,[1 1 nt])./86400 + Vtend.*repmat(AngleCS,[1 1 nt])./86400; 
[ tendency, sign8 ] = circulationInt( xum,yum,tendE,xvm,yvm,tendN,xSquare,ySquare );
disp('1')
[ tendencyC, sign8 ] = circulationInt3( xum,yum,tendE,xvm,yvm,tendN,xSquare,ySquare );

load('momentumPressureSurface.mat')
pressUssh=Utend/86400-AdvU-AbU-UDiss-UDif2a-Uext-UdPdx;
pressVssh=Vtend/86400-AdvV-AbV-VDiss-VDif2a-Vext-VdPdy;
pressUssh2=Utend/86400-AdvU-AbU-UDiss-UDif2b-Uext-UdPdx;
pressVssh2=Vtend/86400-AdvV-AbV-VDiss-VDif2b-Vext-VdPdy;
clear Adv* Ab* *Diss *Dif2* *ext 
%rotate to east-north
press1E=UdPdx.*repmat(AngleCS,[1 1 nt]) - VdPdy.*repmat(AngleSN,[1 1 nt]);  
press1N=UdPdx.*repmat(AngleSN,[1 1 nt]) + VdPdy.*repmat(AngleCS,[1 1 nt]); 
press2E=pressUssh.*repmat(AngleCS,[1 1 nt]) - pressVssh.*repmat(AngleSN,[1 1 nt]);  
press2N=pressUssh.*repmat(AngleSN,[1 1 nt]) + pressVssh.*repmat(AngleCS,[1 1 nt]); 
pressE=(pressUssh+UdPdx).*repmat(AngleCS,[1 1 nt]) - (pressVssh+VdPdy).*repmat(AngleSN,[1 1 nt]);  
pressN=(pressUssh+UdPdx).*repmat(AngleSN,[1 1 nt]) + (pressVssh+VdPdy).*repmat(AngleCS,[1 1 nt]); 
pressE2=(pressUssh2+UdPdx).*repmat(AngleCS,[1 1 nt]) - (pressVssh2+VdPdy).*repmat(AngleSN,[1 1 nt]);  
pressN2=(pressUssh2+UdPdx).*repmat(AngleSN,[1 1 nt]) + (pressVssh2+VdPdy).*repmat(AngleCS,[1 1 nt]); 
% pressure terms
[ pressureTrop, sign7 ] = circulationInt( xum,yum,press2E,xvm,yvm,press2N,xSquare,ySquare );
[ pressureClin, sign9 ] = circulationInt( xum,yum,press1E,xvm,yvm,press1N,xSquare,ySquare );
[ pressureClinc, ~ ] = circulationInt3( xum,yum,press1E,xvm,yvm,press1N,xSquare,ySquare );

[ pressure, sign10 ] = circulationInt( xum,yum,pressE,xvm,yvm,pressN,xSquare,ySquare );
[ pressure2, sign20 ] = circulationInt( xum,yum,pressE2,xvm,yvm,pressN2,xSquare,ySquare );%negligible differences!
%%
figure; plot(-tendency,'k','linewidth',2); hold all; plot(wind,'linewidth',2);plot(diffusion,'linewidth',2);plot(dissipation,'linewidth',2);
plot(advection,'linewidth',2);plot(pressureClin,'linewidth',2);plot(timestep,'linewidth',2);
plot(-tendency+wind+dissipation+advection+timestep+pressureClin+diffusion,'g--','linewidth',2)
legend('-d/dt','wind','vertical diffusion','dissipation','advection','baroclinic pressure','timestep','total')

figure; plot(-tendency,'k','linewidth',2); hold all; plot(wind,'linewidth',2);plot(diffusion,'linewidth',2);plot(dissipation,'linewidth',2);
plot(advection,'linewidth',2);plot(pressure,'linewidth',2);plot(timestep,'linewidth',2);
plot(-tendency+wind+dissipation+advection+timestep+pressure+diffusion,'g--','linewidth',2)
legend('-d/dt','wind','vertical diffusion','dissipation','advection','full pressure','timestep','total')
%%
load('vorticitySurface.mat', 'vorta','vort')
load('distancesAreas')
rAz=reshape(raz,[700 200]);
vortInt=squeeze(nansum(nansum(vort(:,:,1:149).*repmat(vortLogic.*rAz.*hFacC(:,:,1),[1 1 149]))));
vortaInt=squeeze(nansum(nansum(vorta(:,:,1:149).*repmat(vortLogic.*rAz.*hFacC(:,:,1),[1 1 149]))));

dcircdt=tendency; dcircdt(2:end-1)=(integral1(3:end)-integral1(1:end-2))./(2*86400);
dvortdt=tendency; dvortdt(2:end)=(vortInt(3:end)-vortInt(1:end-2))./(2*86400);
dvortadt=tendency; dvortadt(2:end)=(vortaInt(3:end)-vortaInt(1:end-2))./(2*86400);

figure; plot(-tendency+wind+dissipation+advection+timestep+pressureClin+diffusion,'--','linewidth',2)
hold all; plot(-dcircdt+wind+dissipation+advection+timestep+pressureClin+diffusion,'--','linewidth',2)
hold all; plot(-dvortadt+wind+dissipation+advection+timestep+pressureClin+diffusion,'--','linewidth',2)
hold all; plot(-dvortdt+wind+dissipation+advection+timestep+pressureClin+diffusion,'--','linewidth',2)
legend('tendency','circ','daily vort','vort snap')
title('vorticity budget, new square, remainder')
%% sets for error
windset=[wind windc];%wind100 windc wind100c];
windS=std(windset,0,2);

tendset=[tendency tendencyC dcircdt dvortdt];%tend100 tendencyC tend100C];
tendS=std(tendset,0,2);

dissset=[dissipation dissipationc];
dissS=std(dissset,0,2);

advset=[advection advectionc];
advS=std(advset,0,2);

timeset=[timestep timestepc];
timeS=std(timeset,0,2);

clinset=[pressureClin pressureClinc];
clinS=std(clinset,0,2);

diffset=[diffusion diffusionc diffusion2 diffusion2c];
diffS=std(diffset,0,2);

totS=tendS+windS+dissS+advS+timeS+clinS+diffS;
figure; plot(-tendency+wind+dissipation+advection+timestep+pressureClin+diffusion)
hold all; plot(totS,'k--'); plot(-totS,'k--')
