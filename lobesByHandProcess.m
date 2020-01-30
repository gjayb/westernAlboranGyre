% process for finding and tracking each lobe. based on lobes found using
% lobes2.m, but often these repeat. I want only lobes that exist for at
% least 2 days and that do not self-intersect. This can be run for each
% simulation day.

%% load 
%load('lobes8day148DailyVSurfaceV2.mat');
%load('lobesByHand.mat')

load('analyzedIso27Man14dayFew2.mat','*all')
lontrFall27=lontrFall;
lattrFall27=lattrFall;
lontrBall27=lontrBall;
lattrBall27=lattrBall;

load('analyzedIso26Man8dayV2.mat', '*all')
lontrFall26=lontrFall;
lattrFall26=lattrFall;
lontrBall26=lontrBall;
lattrBall26=lattrBall;

load('lobesByHandSigma275.mat')
% lobeNames27=lobeNames;
% xlobeH27=xlobeH;
% ylobeH27=ylobeH;
% lobeAhand27=lobeAhand;
% clear xlobeHm ylobeHm *LobeH

load('analyzedIso275Man14dayFew2.mat')
%load('lobes14dayIso265days1to50small.mat')
%load('lobes14dayIso265.mat')
clear Sigma CT Angle* density* open* wagR np* notaBene nt2 nt3 perimeters areas
clear delZ dist* clusterorder ans *U *V *G C DRC ii j jc jcN jj kk k
% load('analyzedIso26Man8dayV2.mat')
% clear Sigma CT Angle* density* open* wagR
%% create variables- only needed 1st time
% clear lobeNames lobeAhand xlobeH ylobeH insideLobeH
% lobeNames(1:141,1)=cellstr('not set'); %to get back string, char(lobeNames(1,1))
% lobeAhand(1:141,1)=nan;
% xlobeH(1:141,1)=lobeAhand;
% ylobeH=xlobeH;
% insideLobeH=zeros([700 200 141]);

%% setup for the given day
i=20; %simulation day is i+8? doesn't seem like it for analyzedIso_Man8day.mat
        %but it IS for V2!! NO it isnt?!
for ii=i-1:i+1
    figure; plot(lontrFall(ii,:),lattrFall(ii,:),'r')
    hold on; plot(lontrBall(ii,:),lattrBall(ii,:),'b')
    plot(lonCoast,latCoast,'k')
    axis([-6.5 -2 35 37])
    title(num2str(ii))
end

% xlobeHS(xlobeHS==0)=nan;
% ylobeHS(ylobeHS==0)=nan;
% figure; plot(squeeze(xlobeHS(i-8,:,:)).',squeeze(ylobeHS(i-8,:,:)).')
% title('surface lobes')

xlobeH26(xlobeH26==0)=nan;
ylobeH26(ylobeH26==0)=nan;
figure; plot(squeeze(xlobeH26(i,:,:)).',squeeze(ylobeH26(i,:,:)).')
title('\sigma=26')

% xlobeH265(xlobeH265==0)=nan;
% ylobeH265(ylobeH265==0)=nan;
% figure; plot(squeeze(xlobeH265(i,:,:)).',squeeze(ylobeH265(i,:,:)).')
% title('\sigma=26.5')

xlobeH265(xlobeH27==0)=nan;
ylobeH265(ylobeH27==0)=nan;
figure; plot(squeeze(xlobeH27(i,:,:)).',squeeze(ylobeH27(i,:,:)).')
title('\sigma=27')

figure; plot(lontrFall27(i,:),lattrFall27(i,:),'r')
    hold on; plot(lontrBall27(i,:),lattrBall27(i,:),'b')
    plot(lonCoast,latCoast,'k')
    axis([-6.5 -2 35 37])
    title('\sigma=27 manifolds')
    
    figure; plot(lontrFall26(i,:),lattrFall26(i,:),'r')
    hold on; plot(lontrBall26(i,:),lattrBall26(i,:),'b')
    plot(lonCoast,latCoast,'k')
    axis([-6.5 -2 35 37])
    title('\sigma=26 manifolds')

nji=length(find(lobAfew(i,:)>0));

clear xlobes ylobes xlobems ylobems

if nji>0
    for j=1:nji
        if fiintFew(i,j,1)>0 & compactness(i,j)<0.1
            if fiintFew(i,j,2)>fiintFew(i,j,1)
                    if biintFew(i,j,2)>biintFew(i,j,1)
                        xlobe=[lontrFall(i,fiintFew(i,j,1):fiintFew(i,j,2)) lontrBall(i,biintFew(i,j,2):-1:biintFew(i,j,1))];
                        ylobe=[lattrFall(i,fiintFew(i,j,1):fiintFew(i,j,2)) lattrBall(i,biintFew(i,j,2):-1:biintFew(i,j,1))];
                    else
                        xlobe=[lontrFall(i,fiintFew(i,j,1):fiintFew(i,j,2)) lontrBall(i,biintFew(i,j,2):biintFew(i,j,1))];
                        ylobe=[lattrFall(i,fiintFew(i,j,1):fiintFew(i,j,2)) lattrBall(i,biintFew(i,j,2):biintFew(i,j,1))];
                    end
            elseif biintFew(i,j,2)>biintFew(i,j,1)
                    xlobe=[lontrFall(i,fiintFew(i,j):-1:fiintFew(i,j,2)) lontrBall(i,biintFew(i,j,2):-1:biintFew(i,j,1))];
                    ylobe=[lattrFall(i,fiintFew(i,j):-1:fiintFew(i,j,2)) lattrBall(i,biintFew(i,j,2):-1:biintFew(i,j,1))];
            else
                    xlobe=[lontrFall(i,fiintFew(i,j):-1:fiintFew(i,j,2)) lontrBall(i,biintFew(i,j,2):biintFew(i,j,1))];
                    ylobe=[lattrFall(i,fiintFew(i,j):-1:fiintFew(i,j,2)) lattrBall(i,biintFew(i,j,2):biintFew(i,j,1))];
            end
            %inside(:,:,i)=inside(:,:,i)+inpolygon(XC,YC,xlobe,ylobe);

        end
 
       figure; plot(xlobe,ylobe); title(num2str(j))
       np1=length(xlobe);
       xlobes(j,1:np1)=xlobe;
       ylobes(j,1:np1)=ylobe;
    end
    xlobems=111000*cosd(ylobes).*(xlobes-xmin*ones(size(xlobes)));
    ylobems=111000*(ylobes-ymin*ones(size(ylobes)));
end

%% alternate
% clear xlobes ylobes
np1=length([2825:2939 191:248]);
xlobes(9,1:np1)=[lontrBall(i,2825:2939) lontrFall(i,191:248)];%lontrBall(i,781:856)];
ylobes(9,1:np1)=[lattrBall(i,2825:2939) lattrFall(i,191:248)];%lattrBall(i,781:856)];

% np1=length(473:538);
% xlobes(9,1:np1)=lontrBall(i,473:538);
% ylobes(9,1:np1)=lattrBall(i,473:538);
% 
xlobems=111000*cosd(ylobes).*(xlobes-xmin*ones(size(xlobes)));
    ylobems=111000*(ylobes-ymin*ones(size(ylobes)));
%% working to set results

j=2;
j2=8;
lobeNames(i,j)=cellstr('L0'); %to get back string, char(lobeNames(1,1))

% indices1=[77:87];
% np=length(indices1); %specified
% xlobeHm(i,j,:)=0; ylobeHm(i,j,:)=0;
% xlobeHm=xlobems(j2,indices1); ylobeHm=ylobems(j2,indices1);
% xlobeH(i,j,:)=0; ylobeH(i,j,:)=0;
% xlobeH(i,j,1:np)=xlobes(j2,indices1);%1:np);
% ylobeH(i,j,1:np)=ylobes(j2,indices1);

% j3=13;
% indices2=[80:147];
% np2=length(indices2); %specified
% xlobeHm(i,j,np+1:np+np2)=xlobems(j3,indices2); 
% ylobeHm(i,j,np+1:np+np2)=ylobems(j3,indices2);
% xlobeH(i,j,np+1:np+np2)=xlobes(j3,indices2);%1:np);
% ylobeH(i,j,np+1:np+np2)=ylobes(j3,indices2);

%full
xlobeHm=xlobems(j2,:); ylobeHm=ylobems(j2,:);
np=find(ylobes(j2,:)>0,1,'last'); %full, found by lobes2.m
xlobeH(i,j,:)=0;
ylobeH(i,j,:)=0;
xlobeH(i,j,1:np)=xlobes(j2,1:np);%1:np);
ylobeH(i,j,1:np)=ylobes(j2,1:np);


lobeAhand(i,j)=polyarea(xlobeHm(1:np),ylobeHm(1:np));

%probably need to redo this once done with each curve of *lobeH
insideLobeH(:,:,i)=insideLobeH(:,:,i)+inpolygon(XC,YC,squeeze(xlobeH(i,j,1:np)),squeeze(ylobeH(i,j,1:np)));
%% by hand
clear xlobes ylobes
%i=23;
j=4;
iF=1953:-1:1856;
iB=978:1086;
np1=length([iB iF]);
xlobes(j,1:np1)=[lontrBall(i,iB) lontrFall(i,iF)];
ylobes(j,1:np1)=[lattrBall(i,iB) lattrFall(i,iF)];
xlobeH(i,j,:)=0;
ylobeH(i,j,:)=0;
xlobeH(i,j,1:np1)=xlobes(j,1:np1);
ylobeH(i,j,1:np1)=ylobes(j,1:np1);
xlobems=111000*cosd(ylobes).*(xlobes-xmin*ones(size(xlobes)));
    ylobems=111000*(ylobes-ymin*ones(size(ylobes)));
xlobeHm(i,j,:)=0;
ylobeHm(i,j,:)=0;
xlobeHm(i,j,1:np1)=xlobems(j,1:np1);
ylobeHm(i,j,1:np1)=ylobems(j,1:np1);
lobeAhand(i,j)=polyarea(xlobeHm(1:np1),ylobeHm(1:np1));
insideLobeH(:,:,i)=insideLobeH(:,:,i)+inpolygon(XC,YC,squeeze(xlobeH(i,j,1:np1)),squeeze(ylobeH(i,j,1:np1)));
lobeNames(i,j)=cellstr('b');

%% plot lobe edges
xlobeH(xlobeH==0)=nan;
ylobeH(ylobeH==0)=nan;
figure; plot(squeeze(xlobeH(i,:,:)).',squeeze(ylobeH(i,:,:)).')
%% close
close all
%% saving
lobeAhand(lobeAhand==0)=nan;
insideLobeH2=double(sign(insideLobeH));
figure; pcolor(XC,YC,sum(insideLobeH2,3)); shading 'flat'
nLobeH=nansum(sign(lobeAhand),2);
areaTotLobeH=nansum((lobeAhand),2);
save('lobesByHandSigma275.mat','lobeNames*','lobeAhand*','*obeH*')

%% redo insideLobeH
insideLobeH=zeros([700 200 138]);
for i=1:18
    i
    for j=1:4
        np=find(ylobeH(i,j,:)>0,1,'last');
        if np>1
        insideLobeH(:,:,i)=insideLobeH(:,:,i)+inpolygon(XC,YC,squeeze(xlobeH(i,j,1:np)),squeeze(ylobeH(i,j,1:np)));
        end%if np
    end%for j
end%for i