%allows remaking of lonWagClosed,latWagClosed by hand
%% closing surface curves- get curves
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
load('wagSurfaceClosedWorking.mat')
for i=14:146
   np=find(latWagClosed(i,:)>0,1,'last');
   if np>10
       disp(num2str(i))
       latWagClosedS(i,1:np)=latWagClosed(i,1:np);
       lonWagClosedS(i,1:np)=lonWagClosed(i,1:np);
   end
end
latWagClosed=latWagClosedS; 
lonWagClosed=lonWagClosedS;
latWagClosed(latWagClosed==0)=nan;
lonWagClosed(latWagClosed==0)=nan;
%%
iscoast=ismember(lonWagClosed,lonCoast)&ismember(latWagClosed,latCoast);
figure; pcolor(double(iscoast)); shading 'flat'
%%
a1=190;
a2=390;
lonHold=lonWagClosed(i,a1:a2);
latHold=latWagClosed(i,a1:a2);
np=length(a1:a2);
lonWagClosed(i,1:np)=lonHold;
latWagClosed(i,1:np)=latHold;
lonWagClosed(i,np+1:end)=nan;
latWagClosed(i,1+np:end)=nan;
%%
lonWagClosed(i,np+1)=lonWagClosed(i,1);
latWagClosed(i,np+1)=latWagClosed(i,1);
%%

lonWagClosed(i,1:length(coasti:coastj))=lonCoast(coasti:coastj);
latWagClosed(i,1:length(coasti:coastj))=latCoast(coasti:coastj);
lonWagClosed(i,length(coasti:coastj)+1:length(coasti:coastj)+length(lonHold))=lonHold;
latWagClosed(i,length(coasti:coastj)+1:length(coasti:coastj)+length(lonHold))=latHold;
lonWagClosed(i,210)=lonWagClosed(i,1);
latWagClosed(i,210)=latWagClosed(i,1);
lonWagClosed(i,211:end)=nan;
latWagClosed(i,211:end)=nan;
%% closing surface curves- plot, find nearest coast points
i=41;%
   np=find(latWagClosed(i,:)>0,1,'last')
figure; plot(lonCoast,latCoast,'k'); hold all; plot(lonWagClosed(i,:),latWagClosed(i,:))
scatter(lonWagClosed(i,:),latWagClosed(i,:),25,1:length(lonWagClosed),'filled'); colorbar
[~,coasti]=min(sqrt((lonCoast(1:3000)-lonWagClosed(i,1)).^2+(latCoast(1:3000)-latWagClosed(i,1)).^2))
[~,coastj]=min(sqrt((lonCoast(1:3000)-lonWagClosed(i,np)).^2+(latCoast(1:3000)-latWagClosed(i,np)).^2))
%%
lonWagClosed(i,np+1:length(coastj:coasti)+np)=lonCoast(coastj:coasti);
latWagClosed(i,np+1:length(coastj:coasti)+np)=latCoast(coastj:coasti);
np2=length(coastj:coasti)+np;
lonWagClosed(i,np2+1:end)=nan;
latWagClosed(i,np2+1:end)=nan;
%%
empties=[39];
%%
close all
save('sigma275ClosedNF.mat','lonWagClosed','latWagClosed','empties')
%%
%did surface days 14-68 and 144-146 on Jan 23
%%

load('manifoldsSurfaceInt14NFanalyzed.mat')
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
load('geometrySpinupSteady')
[nt,~]=size(lontrFall)

lontrFall(lontrFall==0)=nan;
lontrBall(lontrBall==0)=nan;
lattrFall(lattrFall==0)=nan;
lattrBall(lattrBall==0)=nan;

xmin=min(min(XC)); ymin=min(min(YC));
lonWagClosed(latWagClosed==0)=nan;
latWagClosed(latWagClosed==0)=nan;
%%
daysToCheck=unique([37 39:48 51 52 63 79:81 84 85 87 88 106 107 109 110 117:121 144:146 29 38 39 52:54 63 70 71 79 81 82 84:89 106:111 145 146 24 29 38 63 71:75 80 81 85 89 97 102 110 117 133 144 146]);

%%
load('wagSurfaceClosedWorking.mat')
%%
i=146
figure; plot(lonCoast,latCoast); hold all
plot(lontrFall(i,:),lattrFall(i,:))
plot(lontrBall(i,:),lattrBall(i,:))
plot(lonWagClosed(i,:),latWagClosed(i,:),'--')
title(num2str(i))
[xif,yif,iif]=polyxpoly(lonCoast(575:700),latCoast(575:700),lontrFall(i,5000:6000),lattrFall(i,5000:6000));
[xib,yib,iib]=polyxpoly(lonCoast(700:1000),latCoast(700:1000),lontrBall(i,5000:6000),lattrBall(i,5000:6000));
[~,~,iif2]=polyxpoly([-4.4 -4.4],[36.6 35.8],lontrFall(i,5000:6000),lattrFall(i,5000:6000));
[~,~,iib2]=polyxpoly([-4.4 -4.4],[37 35.8],lontrBall(i,5000:6000),lattrBall(i,5000:6000));
iifu=unique(iif(:,2))
iif2u=unique(iif2(:,2))
%%
for j=1:length(iifu)
    for k=1:length(iif2u)
        figure; plot(lonCoast(500:1000),latCoast(500:1000),'k')
        hold all; plot(lonWagClosed(i,:),latWagClosed(i,:))
        if iif2u(k)>iifu(j)
            plot(lontrFall(i,5000+iifu(j):5000+iif2u(k)),lattrFall(i,5000+iifu(j):5000+iif2u(k)))
        else
            plot(lontrFall(i,5000+iif2u(k):5000+iifu(j)),lattrFall(i,5000+iif2u(k):5000+iifu(j)))
        end
        tstr=strcat('forward iif ',num2str(j)',', iif2 ',num2str(k));
        title(tstr)
    end
end
%%
iibu=unique(iib(:,2))
iib2u=unique(iib2(:,2))
%%
for j=1:length(iibu)
    for k=1:length(iib2u)
        figure; plot(lonCoast(500:1000),latCoast(500:1000),'k')
        hold all; plot(lonWagClosed(i,:),latWagClosed(i,:))
        if iib2u(k)>iibu(j)
            plot(lontrBall(i,5000+iibu(j):5000+iib2u(k)),lattrBall(i,5000+iibu(j):5000+iib2u(k)))
        else
            plot(lontrBall(i,5000+iib2u(k):5000+iibu(j)),lattrBall(i,5000+iib2u(k):5000+iibu(j)))
        end
        tstr=strcat('backward iib ',num2str(j)',', iib2 ',num2str(k));
        title(tstr)
    end
end

%%
fis=5000+iif2u(2):1:5000+iifu(1);
bis=5000+iibu(1):1:5000+iib2u(2);
lon1=[lontrBall(i,bis) lontrFall(i,fis)];
lat1=[lattrBall(i,bis) lattrFall(i,fis)];
figure; plot(lon1,lat1)
%%
lonWagClosed(i,1:length(lon1))=lon1;
latWagClosed(i,1:length(lat1))=lat1;
lonWagClosed(i,1+length(lon1):end)=nan;
latWagClosed(i,1+length(lat1):end)=nan;
%%
