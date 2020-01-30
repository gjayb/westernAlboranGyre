addpath('../mStuff')

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

load('wagSurfaceClosedWorking.mat','*Closed')
for i=14:146
    np=find(latWagClosed(i,:)>0,1,'last');
    if np>10
        lonWagClosedS(i,1:np)=lonWagClosed(i,1:np);
        latWagClosedS(i,1:np)=latWagClosed(i,1:np);
    end
end
latWagClosed=latWagClosedS; clear latWagClosedS
lonWagClosed=lonWagClosedS; clear lonWagClosedS

latWagClosed(latWagClosed==0)=nan;
lonWagClosed(latWagClosed==0)=nan;
 xm=111000*cosd(latWagClosed).*(lonWagClosed-xmin*ones(size(lonWagClosed)));
 ym=111000*(latWagClosed-ymin*ones(size(latWagClosed)));
 xcm=111000*cosd(YC).*(XC-xmin*ones(size(XC)));
 ycm=111000*(YC-ymin*ones(size(YC)));
 xcoast=111000*cosd(latCoast).*(lonCoast-xmin*ones(size(lonCoast)));
 ycoast=111000*(latCoast-ymin*ones(size(latCoast)));
 %% remove coast points
 coastFlag=ismember(lonWagClosed,lonCoast)&ismember(latWagClosed,latCoast);
 lonWagClosed2=lonWagClosed.*double(~coastFlag);
 latWagClosed2=latWagClosed.*double(~coastFlag);
for j=14:146
    np1=find(latWagClosed2(j,:)>0,1,'first');
    np2=find(latWagClosed2(j,:)>0,1,'last');
   lonWagClosed3(j,1:(np2-np1+1))=lonWagClosed2(j,np1:np2);
   latWagClosed3(j,1:(np2-np1+1))=latWagClosed2(j,np1:np2);
end
xm=111000*cosd(latWagClosed3).*(lonWagClosed3-xmin*ones(size(lonWagClosed3)));
 ym=111000*(latWagClosed3-ymin*ones(size(latWagClosed3)));
 xm(xm<min(xcm(:)))=NaN;
 xm(xm>max(xcm(:)))=NaN;
 ym(ym<min(ycm(:)))=NaN;
 ym(ym>max(ycm(:)))=NaN;
%%
load('uvwDailyDepth1interpolatedNF.mat', 'XC','YC','u','v','xvel','yvel','tvel','*Coast')
tvel

sizeU=size(u)
%%
dlmin=2000; %RR/3; %10^(-3);
for j=14:145
    j
    np=find(ym(j,:)>0,1,'last');
    if np>1
        xunst=xm(j,1:np);%[xm(j,1:np) xm(j,1)]; %end of previous integration, make a loop
        yunst=ym(j,1:np);%[ym(j,1:np) ym(j,1)]; %no loop when no coast!

        tunst=(1:1:length(xunst)); t1unst=(1:0.005:length(xunst));
        x1unst=interp1(tunst,xunst,t1unst); y1unst=interp1(tunst,yunst,t1unst); %interpolate to 200x as many points

        dxunst=diff(x1unst); dyunst=diff(y1unst); dlunst=sqrt(dxunst.^2+dyunst.^2); %distances between points

        nunst=[];
        ii=1; 
        while ii<(length(dlunst)-1) %while index ii is less than the number of points
            if dlunst(ii)>dlmin
                nunst=[nunst ii]; ii=ii+1; %if this distance is bigger than the minimum, add ii to list nunst
            else 
                jj=1;
                while ((dlunst(ii)+dlunst(ii+jj)<dlmin)&&((ii+jj)<length(dlunst)))
                      dlunst(ii)=dlunst(ii)+dlunst(ii+jj); jj=jj+1;  %otherwise, find how many points forward you can move until it gets too far
                 end
                nunst=[nunst ii+jj]; ii=ii+jj; %add the next index, where the distance just goes over dlmin, to nunst
            end
        %     if mod(ii,10000)==0
        %         disp(num2str(ii))
        %     end
        end
        x0b=x1unst(nunst); y0b=y1unst(nunst); 
        %set input to next integration to all the interpolated points that are needed to have separations no larger than dlmin
            z0f=[x0b;y0b];
        disp('interp done')
        t1=86400*(j):86400:(j+1)*86400;
            options=odeset('RelTol',10^(-10),'AbsTol',10^(-13));
            [~,zz]=ode45(@HamEqSolver_BiLin_Irina,t1,z0f(:),options,u,v,xvel,yvel,tvel);%expects 1D xvel yvel tvel where u,v are defined on the meshgrid
        xtr(j+1,1:length(zz(end,1:2:end)))=zz(end,1:2:end-1); %end of previous integration
        ytr(j+1,1:length(zz(end,1:2:end)))=zz(end,2:2:end);
        x0(j+1,1:length(zz(1,1:2:end)))=zz(1,1:2:end-1); %start of previous integration
        y0(j+1,1:length(zz(1,1:2:end)))=zz(1,2:2:end);   
    end
end
lattrF=ones(size(ytr)).*ymin+ytr./111000;
lontrF=ones(size(xtr)).*xmin+xtr./(111000.*cosd(lattrF));
lattr0=ones(size(y0)).*ymin+y0./111000;
lontr0=ones(size(x0)).*xmin+x0./(111000.*cosd(lattr0));
save('surfaceEdgeCheckNocoast.mat','*Closed','*tr*','x*','y*')

%% size of errors
xtr(xtr<min(xcm(:)))=NaN;
 xtr(xtr>max(xcm(:)))=NaN;
 ytr(ytr<min(ycm(:)))=NaN;
 ytr(ytr>max(ycm(:)))=NaN;
 
dx=diff(xtr,1,2);
dy=diff(ytr,1,2);
ds=sqrt(dx.^2+dy.^2);
xc=0.5*xtr(:,1:end-1)+0.5*xtr(:,2:end);
yc=0.5*ytr(:,1:end-1)+0.5*ytr(:,2:end);
lineX1=xc-8640.*dy./ds;
lineX2=xc+8640.*dy./ds;
lineY1=yc+8640.*dx./ds;
lineY2=yc-8640.*dx./ds;
manDist=zeros(size(xc));
multipleManDist=manDist;
manSign=manDist;
%%
for i=15:146
    i
    np=find(xc(i,:)>0,1,'last');
    np2=find(xm(i,:)>0,1,'last');
    for k=1:np
        [xhold,yhold]=polyxpoly([lineX1(i,k) lineX2(i,k)],[lineY1(i,k) lineY2(i,k)],xm(i,1:np2),ym(i,1:np2));
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
dz=5;
vol=nansum(dz.*manDist.*manSign.*ds./86400,2); %dz!!!!
%save('surfaceEdgeCheck2nocoast.mat','vol','manDist','manSign','multipleManDist','d*','x*','y*')
%% find coast part
for i=15:146
    np=find(x0(i,:)>0,1,'last');
    if np>1
    for k=1:np
        dxhold=xcoast(500:1300)-x0(i,k);
        dyhold=ycoast(500:1300)-y0(i,k);
        [coastDist(i,k),~]=min(sqrt(dxhold.^2+dyhold.^2));
    end
    end
end
isCoast=(coastDist<2000.1).*(x0>0); isCoast=logical(isCoast);
volNoCoast=nansum(dz.*manDist.*manSign.*ds.*double(~(isCoast(:,1:end-1)&isCoast(:,2:end)))./86400,2);
volCoast=nansum(dz.*manDist.*manSign.*ds.*double(isCoast(:,1:end-1)&isCoast(:,2:end))./86400,2);
%% find gate part
%xm,ym and xc,yc are coast, backwards manifold, forwards manifold
gateAffected=false(size(x0));
gate1=zeros(1,146); gate2=gate1;
x0(x0==0)=nan;
y0(y0==0)=nan;
xtr(xtr==0)=nan;
ytr(ytr==0)=nan;
x0(x0<min(xcm(:)))=NaN;
 x0(x0>max(xcm(:)))=NaN;
 y0(y0<min(ycm(:)))=NaN;
 y0(y0>max(ycm(:)))=NaN;
lattrF=ones(size(ytr)).*ymin+ytr./111000;
lontrF=ones(size(xtr)).*xmin+xtr./(111000.*cosd(lattrF));
lattr0=ones(size(y0)).*ymin+y0./111000;
lontr0=ones(size(x0)).*xmin+x0./(111000.*cosd(lattr0));
for i=15:145
    holdvar=find(lontr0(i+1,:)>-4.4,1,'last');
    holdvar2=find(lontrF(i,:)>-4.4,1,'last')+1;
    if ~isempty(holdvar)
        gate1(i)=holdvar;
        if ~isempty(holdvar2)
            gate2(i)=holdvar2;
            if gate2(i)>gate1(i)
                i
                gateAffected(i,gate1(i):gate2(i))=1;
            else
                gateAffected(i,gate2(i):gate1(i))=1;
            end
        end
    elseif ~isempty(holdvar2)
        gate2(i)=holdvar2;
    end
end
volNoGate=nansum(dz.*manDist.*manSign.*ds.*double(~gateAffected(:,2:end))./86400,2);
volGate=nansum(dz.*manDist.*manSign.*ds.*double(gateAffected(:,2:end))./86400,2);
%%
volErr=nansum(dz.*manDist.*manSign.*ds.*double(~(gateAffected(:,2:end)|isCoast(:,1:end-1)|isCoast(:,2:end)))./86400,2);
%% 0.01m/s error equivalent
volErr1cms=nansum(dz.*864.*manSign.*ds.*double(~(gateAffected(:,2:end)))./86400,2);
volErrn1cms=nansum(dz.*864.*-manSign.*ds.*double(~(gateAffected(:,2:end)))./86400,2);
vol1cmsSignedP=max(volErr1cms,volErrn1cms);
vol1cmsSignedN=min(volErr1cms,volErrn1cms);
volErrP1cms=nansum(dz.*864.*1.*ds.*double(~(gateAffected(:,2:end)))./86400,2);
volErrM1cms=nansum(dz.*864.*-1.*ds.*double(~(gateAffected(:,2:end)))./86400,2);
%
%% actual depths
load('isoDepthsNF.mat', 'isoDepth263')
load('uvwSSHDailyDepth1rotated148F.mat', 'SSHa')
for i=15:146
    h=griddata(xcm,ycm,SSHa(:,:,i),xc(i,:),yc(i,:));
    d=griddata(xcm,ycm,isoDepth263(:,:,i),xc(i,:),yc(i,:));
    dzAll(i,:)=h+d;
end
%%
vol2=nansum(dzAll.*manDist.*manSign.*ds./86400,2); 
volNoGate2=nansum(dzAll.*manDist.*manSign.*ds.*double(~gateAffected(:,2:end))./86400,2);
volGate2=nansum(dzAll.*manDist.*manSign.*ds.*double(gateAffected(:,2:end))./86400,2);
volNoCoast2=nansum(dzAll.*manDist.*manSign.*ds.*double(~(isCoast(:,1:end-1)&isCoast(:,2:end)))./86400,2);
volCoast2=nansum(dzAll.*manDist.*manSign.*ds.*double(isCoast(:,1:end-1)&isCoast(:,2:end))./86400,2);
volErr2=nansum(dzAll.*manDist.*manSign.*ds.*double(~(gateAffected(:,2:end)|isCoast(:,1:end-1)|isCoast(:,2:end)))./86400,2);
%%
save('surfaceEdgeCheck3nocoastJan18.mat','vol*','dzAll','gateAffected','isCoast')