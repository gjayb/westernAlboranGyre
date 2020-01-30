%lagrange WAG volume
% location
%load('wagAreaAndFluxGate.mat', 'inWAG')%adjust by adding 26, 265, or 27 to filename
load('geometrySpinupSteady.mat','XC','YC','XG','YG','*Coast','d')
%load('wagAreaAndFluxGate275v2.mat', '*Closed')
load('surfaceClosedNF.mat')
xmin=min(XC(:)); ymin=min(YC(:));
xcm=111000*cosd(YC).*(XC-xmin*ones(size(XC)));
 ycm=111000*(YC-ymin*ones(size(YC)));
 xcoast=111000*cosd(latCoast).*(lonCoast-xmin*ones(size(lonCoast)));
 ycoast=111000*(latCoast-ymin*ones(size(latCoast)));
xgm=111000*cosd(YG).*(XG-xmin*ones(size(XG)));
 ygm=111000*(YG-ymin*ones(size(YG)));
 load('distancesAreas','RAC')
gridPolygonsX(:,:,1)=xgm(1:end-1,1:end-1);
gridPolygonsX(:,:,2)=xgm(1:end-1,2:end);
gridPolygonsX(:,:,3)=xgm(2:end,2:end);
gridPolygonsX(:,:,4)=xgm(2:end,1:end-1);
gridPolygonsX(:,:,5)=xgm(1:end-1,1:end-1);
gridPolygonsY(:,:,1)=ygm(1:end-1,1:end-1);
gridPolygonsY(:,:,2)=ygm(1:end-1,2:end);
gridPolygonsY(:,:,3)=ygm(2:end,2:end);
gridPolygonsY(:,:,4)=ygm(2:end,1:end-1);
gridPolygonsY(:,:,5)=ygm(1:end-1,1:end-1);
 %%
 lonWagClosed(latWagClosed==0)=nan;
 latWagClosed(latWagClosed==0)=nan;
 xClosed=111000*cosd(latWagClosed).*(lonWagClosed-xmin*ones(size(lonWagClosed)));
 yClosed=111000*(latWagClosed-ymin*ones(size(latWagClosed)));
 inwagSareas=zeros(700,200,146);
figure;
for i=14:146
    [x2,y2]=poly2cw(xClosed(i,:),yClosed(i,:));
    c=[];
    if latWagClosed(i,1)>0
        np=find(latWagClosed(i,:)>0,1,'last');
        inWAG=inpolygon(XC,YC,lonWagClosed(i,1:np),latWagClosed(i,1:np));
        inWAGll=inpolygon(XG(1:end-1,1:end-1),YG(1:end-1,1:end-1),lonWagClosed(i,1:np),latWagClosed(i,1:np));
        inWAGlu=inpolygon(XG(1:end-1,2:end),YG(1:end-1,2:end),lonWagClosed(i,1:np),latWagClosed(i,1:np));
        inWAGrl=inpolygon(XG(2:end,1:end-1),YG(2:end,1:end-1),lonWagClosed(i,1:np),latWagClosed(i,1:np));
        inWAGru=inpolygon(XG(2:end,2:end),YG(2:end,2:end),lonWagClosed(i,1:np),latWagClosed(i,1:np));
        inWAGScorners(:,:,i)=double(inWAGll)+double(inWAGlu)+double(inWAGrl)+double(inWAGru);
        [c,h]=contour(XC,YC,double(inWAG),[1 1]);
        [c2,h2]=contour(xcm,ycm,double(inWAG),[1 1]);
        
        for j=230:600
        for k=1:199
            if ispolycw(squeeze(gridPolygonsX(j,k,:)),squeeze(gridPolygonsY(j,k,:)))&d(j,k)>0   
                [holdx,holdy]=polybool('intersection',squeeze(gridPolygonsX(j,k,:)),squeeze(gridPolygonsY(j,k,:)),x2,y2);
                if ~isempty(holdx)
                    inwagSareas(j,k,i)=polyarea(holdx,holdy);
                elseif inWAGScorners(j,k,i)==4
                    inwagSareas(j,k,i)=RAC(j,k);
                end
            elseif d(j,k)>0
                [x1,y1]=poly2cw(squeeze(gridPolygonsX(j,k,:)),squeeze(gridPolygonsY(j,k,:)));
                [holdx,holdy]=polybool('intersection',squeeze(gridPolygonsX(j,k,:)),squeeze(gridPolygonsY(j,k,:)),x2,y2);
                if ~isempty(holdx)
                    inwagSareas(j,k,i)=polyarea(holdx,holdy);
                elseif inWAGScorners(j,k,i)==4
                    inwagSareas(j,k,i)=RAC(j,k);
                end
            end
        end
        end
        areaHold=inwagSareas(:,:,i);
        [p1,p2]=find(isnan(inwagSareas(:,:,i))==1);
        for index1=1:length(p1); 
            areaHold(p1(index1),p2(index1))=nanmean([areaHold(p1(index1)-1,p2(index1)) areaHold(p1(index1)+1,p2(index1)) areaHold(p1(index1),p2(index1)+1) areaHold(p1(index1),p2(index1)-1)]); 
        end
        inwagSareas(:,:,i)=areaHold;
    end
%     if ~isempty(c)
%     npointsi=find(c(2,:)==floor(c(2,:)));
%     varhold=max(c(2,npointsi));
%     i1=1+find(c(2,:)==varhold);
%     iend=varhold+i1-1;
%     inWAG275(:,:,i)=inpolygon(XC,YC,c(1,i1:iend),c(2,i1:iend));
%     inWAGll=inpolygon(XG(1:end-1,1:end-1),YG(1:end-1,1:end-1),c(1,i1:iend),c(2,i1:iend));
%     inWAGlu=inpolygon(XG(1:end-1,2:end),YG(1:end-1,2:end),c(1,i1:iend),c(2,i1:iend));
%     inWAGrl=inpolygon(XG(2:end,1:end-1),YG(2:end,1:end-1),c(1,i1:iend),c(2,i1:iend));
%     inWAGru=inpolygon(XG(2:end,2:end),YG(2:end,2:end),c(1,i1:iend),c(2,i1:iend));
%     inWAG275corners2(:,:,i)=double(inWAGll)+double(inWAGlu)+double(inWAGrl)+double(inWAGru);
%     
%     npointsi=find(c2(2,:)==floor(c2(2,:)));
%     varhold=max(c2(2,npointsi));
%     i1=1+find(c2(2,:)==varhold);
%     iend=varhold+i1-1;
%     edgeWagX(i,1:length(i1:iend))=c2(1,i1:iend);
%     edgeWagY(i,1:length(i1:iend))=c2(2,i1:iend);
%     end
end
%inWAG275raw=inWAG;
%figure; pcolor(xcm,ycm,inwagSareas(:,:,16)); shading 'flat'; colorbar
%hold on; plot(xClosed(16,:),yClosed(16,:),'k')

%this works but I need to clean up the lon/latWagClosed curves so they are
%smoother, otherwise too many intersections gives nan area

%clear inWAG inWAGl* inWAGr*

save('surfaceWAGareas.mat','inwagSareas','inWAGScorners','lonWagClosed','latWagClosed')
%% from areas to volume, uses SSH and isoDepth of sigma=26.3 rev method


%%
runthis=false;
if runthis
load('wagAreaAndFluxGate27v2.mat', '*Closed')
figure;
for i=1:146
    c=[];
    if latWagClosed(i,1)>0
        np=find(latWagClosed(i,:)>0,1,'last');
        inWAG=inpolygon(XC,YC,lonWagClosed(i,1:np),latWagClosed(i,1:np));
        [c,h]=contour(XC,YC,double(inWAG),[1 1]);
    end
    if ~isempty(c)
    npointsi=find(c(2,:)==floor(c(2,:)));
    varhold=max(c(2,npointsi));
    i1=1+find(c(2,:)==varhold);
    iend=varhold+i1-1;
    inWAG27(:,:,i)=inpolygon(XC,YC,c(1,i1:iend),c(2,i1:iend));
    inWAGll=inpolygon(XG(1:end-1,1:end-1),YG(1:end-1,1:end-1),c(1,i1:iend),c(2,i1:iend));
    inWAGlu=inpolygon(XG(1:end-1,2:end),YG(1:end-1,2:end),c(1,i1:iend),c(2,i1:iend));
    inWAGrl=inpolygon(XG(2:end,1:end-1),YG(2:end,1:end-1),c(1,i1:iend),c(2,i1:iend));
    inWAGru=inpolygon(XG(2:end,2:end),YG(2:end,2:end),c(1,i1:iend),c(2,i1:iend));
    inWAG27corners(:,:,i)=double(inWAGll)+double(inWAGlu)+double(inWAGrl)+double(inWAGru);
    end
end
clear inWAG inWAGl* inWAGr*
%%
load('wagAreaAndFluxGate2675v2.mat', '*Closed')
figure;
for i=1:146
    c=[];
    if latWagClosed(i,1)>0
        np=find(latWagClosed(i,:)>0,1,'last');
        inWAG=inpolygon(XC,YC,lonWagClosed(i,1:np),latWagClosed(i,1:np));
        [c,h]=contour(XC,YC,double(inWAG),[1 1]);
    end
    if ~isempty(c)
    npointsi=find(c(2,:)==floor(c(2,:)));
    varhold=max(c(2,npointsi));
    i1=1+find(c(2,:)==varhold);
    iend=varhold+i1-1;
    inWAG2675(:,:,i)=inpolygon(XC,YC,c(1,i1:iend),c(2,i1:iend));
    inWAGll=inpolygon(XG(1:end-1,1:end-1),YG(1:end-1,1:end-1),c(1,i1:iend),c(2,i1:iend));
    inWAGlu=inpolygon(XG(1:end-1,2:end),YG(1:end-1,2:end),c(1,i1:iend),c(2,i1:iend));
    inWAGrl=inpolygon(XG(2:end,1:end-1),YG(2:end,1:end-1),c(1,i1:iend),c(2,i1:iend));
    inWAGru=inpolygon(XG(2:end,2:end),YG(2:end,2:end),c(1,i1:iend),c(2,i1:iend));
    inWAG2675corners(:,:,i)=double(inWAGll)+double(inWAGlu)+double(inWAGrl)+double(inWAGru);
    end
end
clear inWAG inWAGl* inWAGr*
%
load('wagAreaAndFluxGate265v2.mat', '*Closed')
figure;
for i=1:146
    c=[];
    if latWagClosed(i,1)>0
        np=find(latWagClosed(i,:)>0,1,'last');
        inWAG=inpolygon(XC,YC,lonWagClosed(i,1:np),latWagClosed(i,1:np));
        [c,h]=contour(XC,YC,double(inWAG),[1 1]);
    end
    if ~isempty(c)
    npointsi=find(c(2,:)==floor(c(2,:)));
    varhold=max(c(2,npointsi));
    i1=1+find(c(2,:)==varhold);
    iend=varhold+i1-1;
    inWAG265(:,:,i)=inpolygon(XC,YC,c(1,i1:iend),c(2,i1:iend));
    inWAGll=inpolygon(XG(1:end-1,1:end-1),YG(1:end-1,1:end-1),c(1,i1:iend),c(2,i1:iend));
    inWAGlu=inpolygon(XG(1:end-1,2:end),YG(1:end-1,2:end),c(1,i1:iend),c(2,i1:iend));
    inWAGrl=inpolygon(XG(2:end,1:end-1),YG(2:end,1:end-1),c(1,i1:iend),c(2,i1:iend));
    inWAGru=inpolygon(XG(2:end,2:end),YG(2:end,2:end),c(1,i1:iend),c(2,i1:iend));
    inWAG265corners(:,:,i)=double(inWAGll)+double(inWAGlu)+double(inWAGrl)+double(inWAGru);
    end
end
clear inWAG inWAGl* inWAGr*
%
load('wagAreaAndFluxGate263v2.mat', '*Closed')
figure;
for i=1:146
    i
    c=[];
    if latWagClosed(i,1)>0
        np=find(latWagClosed(i,:)>0,1,'last');
        inWAG=inpolygon(XC,YC,lonWagClosed(i,1:np),latWagClosed(i,1:np));
        [c,h]=contour(XC,YC,double(inWAG),[1 1]);
    end
    if ~isempty(c)
    npointsi=find(c(2,:)==floor(c(2,:)));
    varhold=max(c(2,npointsi));
    i1=1+find(c(2,:)==varhold);
    iend=varhold+i1-1;
    inWAG263(1:700,1:200,i)=inpolygon(XC,YC,c(1,i1:iend),c(2,i1:iend));
    inWAGll=inpolygon(XG(1:end-1,1:end-1),YG(1:end-1,1:end-1),c(1,i1:iend),c(2,i1:iend));
    inWAGlu=inpolygon(XG(1:end-1,2:end),YG(1:end-1,2:end),c(1,i1:iend),c(2,i1:iend));
    inWAGrl=inpolygon(XG(2:end,1:end-1),YG(2:end,1:end-1),c(1,i1:iend),c(2,i1:iend));
    inWAGru=inpolygon(XG(2:end,2:end),YG(2:end,2:end),c(1,i1:iend),c(2,i1:iend));
    inWAG263corners(:,:,i)=double(inWAGll)+double(inWAGlu)+double(inWAGrl)+double(inWAGru);
    else
        inWAG263(:,:,i)=0;
    end
end
clear inWAG inWAGl* inWAGr*
%
load('wagAreaAndFluxGateS2.mat', '*Closed')
figure;
for i=1:146
    c=[];
    if latWagClosed(i,1)>0
        np=find(latWagClosed(i,:)>0,1,'last');
        inWAG=inpolygon(XC,YC,lonWagClosed(i,1:np),latWagClosed(i,1:np));
        [c,h]=contour(XC,YC,double(inWAG),[1 1]);
    end
    if ~isempty(c)
    npointsi=find(c(2,:)==floor(c(2,:)));
    varhold=max(c(2,npointsi));
    i1=1+find(c(2,:)==varhold);
    iend=varhold+i1-1;
    inWAGS(:,:,i)=inpolygon(XC,YC,c(1,i1:iend),c(2,i1:iend));
    inWAGll=inpolygon(XG(1:end-1,1:end-1),YG(1:end-1,1:end-1),c(1,i1:iend),c(2,i1:iend));
    inWAGlu=inpolygon(XG(1:end-1,2:end),YG(1:end-1,2:end),c(1,i1:iend),c(2,i1:iend));
    inWAGrl=inpolygon(XG(2:end,1:end-1),YG(2:end,1:end-1),c(1,i1:iend),c(2,i1:iend));
    inWAGru=inpolygon(XG(2:end,2:end),YG(2:end,2:end),c(1,i1:iend),c(2,i1:iend));
    inWAGScorners(:,:,i)=double(inWAGll)+double(inWAGlu)+double(inWAGrl)+double(inWAGru);
    end
end
clear inWAG inWAGl* inWAGr*

save('inWagManifoldsCorners.mat','inWAG*')
%
inWAG275corners2(700,200,146)=0;
inWAG275points=double(inWAG275)+inWAG275corners2;
clear inWAG275corners* inWAG275
inWAG27corners(700,200,146)=0;
inWAG27points=double(inWAG27)+inWAG27corners;
clear inWAG27 inWAG27corners
inWAG2675corners(700,200,146)=0;
inWAG2675points=double(inWAG2675)+inWAG2675corners;
clear inWAG2675 inWAG2675corners
inWAG265corners(700,200,146)=0;
inWAG265points=double(inWAG265)+inWAG265corners;
clear inWAG265 inWAG265corners
inWAG263corners(700,200,146)=0;
inWAG263points=double(inWAG263)+inWAG263corners;
clear inWAG263 inWAG263corners
inWAGScorners(700,200,146)=0;
inWAGSpoints=double(inWAGS)+inWAGScorners;
clear inWAGS inWAGScorners
end
%% new volume: areas
inWag3area=zeros([700 200 20 146]);
load('saTCtCpSigma162NF.mat','Sigma')
Sigma=Sigma(:,:,1:20,1:146);
load('surfaceWAGareasOld.mat', 'inwagSareas')
inwagSareas1=inwagSareas;
load('sigma263WAGareasOld.mat', 'inwagSareas')
inwag263areas=inwagSareas;
load('sigma265WAGareasOld.mat', 'inwagSareas')
inwag265areas=inwagSareas;
load('sigma2675WAGareasOld.mat', 'inwagSareas')
inwag2675areas=inwagSareas;
load('sigma27WAGareasOld.mat', 'inwagSareas')
inwag27areas=inwagSareas;
load('sigma275WAGareasOld.mat', 'inwagSareas')
inwag275areas=inwagSareas;
%%
%load('sigma263ClosedNF.mat', 'fillneeded')
%fill263=fillneeded;
%load('sigma265ClosedNF.mat', 'fillneeded')
%fill265=fillneeded;
%load('sigma2675ClosedNF.mat', 'fillneeded')
%fill2675=fillneeded;
%load('sigma27ClosedNF.mat', 'fillneeded')
%fill27=fillneeded;
fill275=find(squeeze(sum(sum(sign(inwag275areas))))==0);
fill275=fill275(14:end);
fill27=find(squeeze(sum(sum(sign(inwag27areas))))==0);
fill27=fill27(14:end);
fill2675=find(squeeze(sum(sum(sign(inwag2675areas))))==0);
fill2675=fill2675(14:end);
fill265=find(squeeze(sum(sum(sign(inwag265areas))))==0);
fill265=fill265(14:end);
fill263=find(squeeze(sum(sum(sign(inwag263areas))))==0);
fill263=fill263(14:end);
fillany=unique([fill263; fill265; fill2675; fill27; fill275]);
%%
for ii=1:length(fill263)
    i=fill263(ii);
    if sum(sum(sign(inwagSareas1(:,:,i))))>0
        inwag263areas(:,:,i)=inwagSareas1(:,:,i);
    elseif sum(sum(sign(inwag265areas(:,:,i))))>0
        inwag263areas(:,:,i)=inwag265areas(:,:,i);
    elseif sum(sum(sign(inwag2675areas(:,:,i))))>0
        inwag263areas(:,:,i)=inwag2675areas(:,:,i);
    elseif sum(sum(sign(inwag27areas(:,:,i))))>0
        inwag263areas(:,:,i)=inwag27areas(:,:,i);
    end
end
for ii=1:length(fill265)
    i=fill265(ii);
    if sum(sum(sign(inwag263areas(:,:,i))))>0
        inwag265areas(:,:,i)=inwag263areas(:,:,i);
    elseif sum(sum(sign(inwag2675areas(:,:,i))))>0
        inwag265areas(:,:,i)=inwag2675areas(:,:,i);
    elseif sum(sum(sign(inwagSareas1(:,:,i))))>0
        inwag265areas(:,:,i)=inwagSareas1(:,:,i);
    elseif sum(sum(sign(inwag27areas(:,:,i))))>0
        inwag265areas(:,:,i)=inwag27areas(:,:,i);
    end
end
for ii=1:length(fill2675)
    i=fill2675(ii);
    if sum(sum(sign(inwag265areas(:,:,i))))>0
        inwag2675areas(:,:,i)=inwag265areas(:,:,i);
    elseif sum(sum(sign(inwag27areas(:,:,i))))>0
        inwag2675areas(:,:,i)=inwag27areas(:,:,i);
    elseif sum(sum(sign(inwag263areas(:,:,i))))>0
        inwag2675areas(:,:,i)=inwag263areas(:,:,i);
    elseif sum(sum(sign(inwagSareas1(:,:,i))))>0
        inwag2675areas(:,:,i)=inwagSareas1(:,:,i);
    elseif sum(sum(sign(inwag275areas(:,:,i))))>0
        inwag2675areas(:,:,i)=inwag275areas(:,:,i);
    end
end
for ii=1:length(fill27)
    i=fill27(ii);
    if sum(sum(sign(inwag2675areas(:,:,i))))>0
        inwag27areas(:,:,i)=inwag2675areas(:,:,i);
    elseif sum(sum(sign(inwag275areas(:,:,i))))>0
        inwag27areas(:,:,i)=inwag275areas(:,:,i);
    elseif sum(sum(sign(inwag265areas(:,:,i))))>0
        inwag27areas(:,:,i)=inwag265areas(:,:,i);
    elseif sum(sum(sign(inwag263areas(:,:,i))))>0
        inwag27areas(:,:,i)=inwag263areas(:,:,i);
    elseif sum(sum(sign(inwagSareas1(:,:,i))))>0
        inwag27areas(:,:,i)=inwagSareas1(:,:,i);
    end
end
for ii=1:length(fill275)
    i=fill275(ii);
    if sum(sum(sign(inwag27areas(:,:,i))))>0
        inwag275areas(:,:,i)=inwag27areas(:,:,i);
        %disp('1')
    elseif sum(sum(sign(inwag2675areas(:,:,i))))>0
        inwag275areas(:,:,i)=inwag2675areas(:,:,i);
        %disp('2')
    elseif sum(sum(sign(inwag265areas(:,:,i))))>0
        inwag275areas(:,:,i)=inwag265areas(:,:,i);
        %disp('3')
    elseif sum(sum(sign(inwag263areas(:,:,i))))>0
        inwag275areas(:,:,i)=inwag263areas(:,:,i);
        %disp('4')
    elseif sum(sum(sign(inwagSareas1(:,:,i))))>0
        inwag275areas(:,:,i)=inwagSareas1(:,:,i);
        %disp('5')
    end
end
disp('fill done')
figure; plot(squeeze(sum(sum(sign(inwag27areas)))))
%%
k=1;
inWag3area(:,:,k,:)=inwagSareas1(:,:,1:146);
for k=2:20%through 16 is enough!
    k
inWag3area(:,:,k,1:146)=double(squeeze(Sigma(:,:,k,:))<26.3).*double(squeeze(Sigma(:,:,k,:))>1).*inwagSareas1(:,:,1:146)...
    +double(squeeze(Sigma(:,:,k,:))<26.5).*double(squeeze(Sigma(:,:,k,:))>=26.3).*inwag263areas(:,:,1:146)...
    +double(squeeze(Sigma(:,:,k,:))<26.75).*double(squeeze(Sigma(:,:,k,:))>=26.5).*inwag265areas(:,:,1:146)...
    +double(squeeze(Sigma(:,:,k,:))<27).*double(squeeze(Sigma(:,:,k,:))>=26.75).*inwag2675areas(:,:,1:146)...
    +double(squeeze(Sigma(:,:,k,:))<27.5).*double(squeeze(Sigma(:,:,k,:))>=27).*inwag27areas(:,:,1:146);
end

inWag3areaB=zeros([700 200 20 146]);
for k=1:20%through 16 is enough!
    k
inWag3areaB(:,:,k,1:146)=double(squeeze(Sigma(:,:,k,:))<26.3).*double(squeeze(Sigma(:,:,k,:))>1).*inwag263areas(:,:,1:146)...
    +double(squeeze(Sigma(:,:,k,:))<26.5).*double(squeeze(Sigma(:,:,k,:))>=26.3).*inwag265areas(:,:,1:146)...
    +double(squeeze(Sigma(:,:,k,:))<26.75).*double(squeeze(Sigma(:,:,k,:))>=26.5).*inwag2675areas(:,:,1:146)...
    +double(squeeze(Sigma(:,:,k,:))<27).*double(squeeze(Sigma(:,:,k,:))>=26.75).*inwag27areas(:,:,1:146)...
    +double(squeeze(Sigma(:,:,k,:))<27.5).*double(squeeze(Sigma(:,:,k,:))>=27).*inwag275areas(:,:,1:146);%...
    %+double(squeeze(Sigma(:,:,k,:))<28).*double(squeeze(Sigma(:,:,k,:))>=27.5).*inwag275areas(:,:,1:131);
end
%save('lagrangeWagAreasFilledOld.mat','inW*','inw*','volW*','dVdt*','-v7.3')
%% volume centered
k=1;
inWag3areaCs(:,:,k,:)=inwagSareas(:,:,1:146);
for k=2:30%through 16 is enough!
    k
inWag3areaCs(:,:,k,1:146)=double(squeeze(Sigma(:,:,k,:))<26.2).*double(squeeze(Sigma(:,:,k,:))>1).*inwagSareas(:,:,1:146)...
    +double(squeeze(Sigma(:,:,k,:))<26.4).*double(squeeze(Sigma(:,:,k,:))>=26.2).*inwag263areas(:,:,1:146)...
    +double(squeeze(Sigma(:,:,k,:))<26.625).*double(squeeze(Sigma(:,:,k,:))>=26.4).*inwag265areas(:,:,1:146)...
    +double(squeeze(Sigma(:,:,k,:))<26.875).*double(squeeze(Sigma(:,:,k,:))>=26.625).*inwag2675areas(:,:,1:146)...
    +double(squeeze(Sigma(:,:,k,:))<27.25).*double(squeeze(Sigma(:,:,k,:))>=26.875).*inwag27areas(:,:,1:146)...
    +double(squeeze(Sigma(:,:,k,:))<27.5).*double(squeeze(Sigma(:,:,k,:))>=27.25).*inwag275areas(:,:,1:146);
end
%% volume centered Hourly
k=1;
inWag3areaCH(:,:,k,:)=inwagSareas(:,:,1:146);
for k=2:30%through 16 is enough!
    k
inWag3areaCH(:,:,k,1:49)=double(squeeze(Sigma(:,:,k,1:49))<26.2).*double(squeeze(Sigma(:,:,k,1:49))>1).*inwagSareas(:,:,1:49)...
    +double(squeeze(Sigma(:,:,k,1:49))<26.4).*double(squeeze(Sigma(:,:,k,1:49))>=26.2).*inwag263areas(:,:,1:49)...
    +double(squeeze(Sigma(:,:,k,1:49))<26.625).*double(squeeze(Sigma(:,:,k,1:49))>=26.4).*inwag265areas(:,:,1:49)...
    +double(squeeze(Sigma(:,:,k,1:49))<26.875).*double(squeeze(Sigma(:,:,k,1:49))>=26.625).*inwag2675areas(:,:,1:49)...
    +double(squeeze(Sigma(:,:,k,1:49))<27.25).*double(squeeze(Sigma(:,:,k,1:49))>=26.875).*inwag27areas(:,:,1:49)...
    +double(squeeze(Sigma(:,:,k,1:49))<27.5).*double(squeeze(Sigma(:,:,k,1:49))>=27.25).*inwag275areas(:,:,1:49);
end
save('lagrangeAreasHourlyDays2021.mat','*area*','ssh148','-v7.3')
%% hourly volume continued
load('geometrySpinupSteady','dInterface')
dZ=diff(dInterface);
k=1;
volLagrangeCH(1,1:49)=squeeze(nansum(nansum((5+ssh148).*squeeze(inWag3areaCH(:,:,1,1:49)))));
for k=2:30
    volLagrangeCH(k,1:49)=squeeze(nansum(nansum(dZ(k).*inWag3areaCH(:,:,k,1:49))));
end
volWagCH=squeeze(nansum(volLagrangeCH));
dVdt=diff(volWagCH)./86400;
figure; plot(dVdt); axis tight

%% new volume: volumes, dV/dt, NB:these are from the above areas using daily-averaged sigma
load('geometrySpinupSteady','dInterface')
dZ=diff(dInterface);
%for k=1:20
%volLagrange(k,1:146)=squeeze(nansum(nansum(dZ(k).*inWag3area(:,:,k,1:146))));
%volLagrangeB(k,1:146)=squeeze(nansum(nansum(dZ(k).*inWag3areaB(:,:,k,1:146))));
%volLagrangeC(k,1:146)=squeeze(nansum(nansum(dZ(k).*inWag3areaCs(:,:,k,1:146))));
%end
for k=1:30
    volLagrangeCs(k,1:146)=squeeze(nansum(nansum(dZ(k).*inWag3areaCs(:,:,k,1:146))));
end
%volWag=squeeze(nansum(volLagrange));
%volWagB=squeeze(nansum(volLagrangeB));
volWagCs=squeeze(nansum(volLagrangeCs));
%dVdt1=diff(volWag)./86400;
%volWagB=squeeze(nansum(volLagrangeB));
%dVdt2=diff(volWagB)./86400;
dVdt3s=diff(volWagCs)./86400;
figure; plot(dVdt1); hold all; plot(dVdt2); plot(dVdt3); plot(dVdt3s); axis tight


%% new volume; plot gate flux and dVdt
%dVdt1(13)=0; dVdt2(13)=0;

%may want 'indexed' versions, not sure
load('gateAdvectionS263indexedFeb.mat','gateFluxS263','yGate','meanDS')%indexed
yGateS=yGate;
load('gateAdvection267527indexedFeb.mat','gateFlux267527','gateFlux2675265','yGate','meanD267527')
yGate2675=yGate;
load('gateAdvection27275indexedFeb.mat','gateFlux27275','gateFlux272675','yGate','meanD27275')
yGate27=yGate;
load('gateAdvection2652675indexedFeb.mat','gateFlux2652675','gateFlux265263','yGate','meanD2652675')
yGate265=yGate;
load('gateAdvection263265indexedFeb.mat','gateFlux263265','gateFlux263S','yGate','meanD263265')
yGate263=yGate;
load('gateAdvection275indexedFeb.mat', 'gateFlux27527')

gateFluxTotFeb=gateFlux263265(2:147)+gateFlux2652675(2:147)+gateFlux267527(2:147)+gateFlux27275(2:147)+gateFluxS263(2:147);
% errS=abs(gateFluxS263).*sqrt((1./meanDS).^2+(1./range(yGateS.')).^2);
% errS(isnan(errS))=0;
% err263=abs(gateFlux263265).*sqrt((1./meanD263265).^2+(1./range(yGate263.')).^2);%...
% err263(isnan(err263))=0;
% err265=abs(gateFlux2652675).*sqrt((1./meanD2652675).^2+(1./range(yGate265.')).^2);
% err265(isnan(err265))=0;
% err2675=abs(gateFlux267527).*sqrt((1./meanD267527).^2+(1./range(yGate2675.')).^2);%...
% err2675(isnan(err2675))=0;
% err27=abs(gateFlux27275).*sqrt((1./meanD27275).^2+(1./range(yGate27.')).^2);
% err27(isnan(err27))=0;
% errGateTot=errS+err263+err265+err2675+err27;
gateFluxTot2Feb=gateFlux265263(1:146)+gateFlux2675265(1:146)+gateFlux272675(1:146)+gateFlux27527(1:146)+gateFlux263S(1:146);
%%
%figure; plot(dVdt1,'linewidth',2); hold all; plot(0:145,gateFluxTot,'linewidth',2); title('WAG top down')
%figure; plot(gateFluxTot(1:145)-dVdt1); hold all; plot(errGateTot,'k--'); plot(-errGateTot,'k--')
%figure; plot(dVdt2,'linewidth',2); hold all; plot(0:145,gateFluxTot2,'linewidth',2); title('WAG bottom up')
figure; plot(dVdt1,'b','linewidth',2); hold all; plot(dVdt1,'c','linewidth',2); plot(gateFluxTot,'r','linewidth',2); plot(0:145,gateFluxToti,'m','linewidth',2)
figure; plot(dVdt2,'b','linewidth',2); hold all; plot(dVdt2,'c','linewidth',2); plot(gateFluxTot2,'r','linewidth',2); plot(0:145,gateFluxTot2i,'m','linewidth',2)
%%
%% diffusive motion of isopycnals at 'steps'
load('diffusiveMovedIsopycnal263.mat')
isoNew=isoDepth;
load('iso263depthNFrev')
dvdt263StepTd=squeeze(nansum(nansum( (inwagSareas1-inwag263areas).*(isoNew(:,:,1:146)-isoDepth(:,:,1:146)) )))./86400;
dvdt263StepBu=squeeze(nansum(nansum( (inwag263areas-inwag265areas).*(isoNew(:,:,1:146)-isoDepth(:,:,1:146)) )))./86400;

load('diffusiveMovedIsopycnal265.mat')
isoNew=isoDepth;
load('iso265depthNFrev')
dvdt265StepTd=squeeze(nansum(nansum( (inwag263areas-inwag265areas).*(isoNew(:,:,1:146)-isoDepth(:,:,1:146)) )))./86400;
dvdt265StepBu=squeeze(nansum(nansum( (inwag265areas-inwag2675areas).*(isoNew(:,:,1:146)-isoDepth(:,:,1:146)) )))./86400;

load('diffusiveMovedIsopycnal2675.mat')
isoNew=isoDepth;
load('iso2675depthNFrev')
dvdt2675StepTd=squeeze(nansum(nansum( (inwag265areas-inwag2675areas).*(isoNew(:,:,1:146)-isoDepth(:,:,1:146)) )))./86400;
dvdt2675StepBu=squeeze(nansum(nansum( (inwag2675areas-inwag27areas).*(isoNew(:,:,1:146)-isoDepth(:,:,1:146)) )))./86400;

load('diffusiveMovedIsopycnal27.mat')
isoNew=isoDepth;
load('iso27depthNFrev')
dvdt27StepTd=squeeze(nansum(nansum( (inwag2675areas-inwag27areas).*(isoNew(:,:,1:146)-isoDepth(:,:,1:146)) )))./86400;
dvdt27StepBu=squeeze(nansum(nansum( (inwag27areas-inwag275areas).*(isoNew(:,:,1:146)-isoDepth(:,:,1:146)) )))./86400;
%% SSH component of volume, dV/dt
% load('uvwSSHDailyDepth1rotated148F.mat', 'SSHa')
% load('sshSnap.mat', 'sshSnap')
volssh1=squeeze(nansum(nansum(SSHa(:,:,1:146).*inwagSareas1)));
volssh1s=squeeze(nansum(nansum(sshSnap(:,:,1:146).*inwagSareas1)));
volssh2=squeeze(nansum(nansum(SSHa(:,:,1:146).*inwag263areas)));
volssh2s=squeeze(nansum(nansum(sshSnap(:,:,1:146).*inwag263areas)));
%% 3 day running mean gate flux and dVdt (incl ssh)
dVdt11=dVdt1+diff(volssh1).'/86400; dVdt11(13)=0;
dVdt3day1=[0 dVdt11(1:end-2)+dVdt11(2:end-1)+dVdt11(3:end) 0]./3;
dVdt1s1=dVdt1+diff(volssh1s).'/86400; dVdt1s1(13)=0;
dVdt3day1s=[0 dVdt1s1(1:end-2)+dVdt1s1(2:end-1)+dVdt1s1(3:end) 0]./3;

dVdt12=dVdt2+diff(volssh2).'/86400; dVdt12(13)=0;
dVdt3day2=[0 dVdt12(1:end-2)+dVdt12(2:end-1)+dVdt12(3:end) 0]./3;
dVdt2s1=dVdt2+diff(volssh2s).'/86400; dVdt2s1(13)=0;
dVdt3day2s=[0 dVdt2s1(1:end-2)+dVdt2s1(2:end-1)+dVdt2s1(3:end) 0]./3;

gateFluxTot3day=[0 gateFluxTot(1:end-2)+gateFluxTot(2:end-1)+gateFluxTot(3:end) 0]./3;
gateFluxTot3dayi=[0 gateFluxToti(1:end-2)+gateFluxToti(2:end-1)+gateFluxToti(3:end) 0]./3;
gateFluxTot3day2=[0 gateFluxTot2(1:end-2)+gateFluxTot2(2:end-1)+gateFluxTot2(3:end) 0]./3;
gateFluxTot3day2i=[0 gateFluxTot2i(1:end-2)+gateFluxTot2i(2:end-1)+gateFluxTot2i(3:end) 0]./3;
gateFluxTot3dayF=[0 gateFluxTotFeb(1:end-2)+gateFluxTotFeb(2:end-1)+gateFluxTotFeb(3:end) 0]./3;
%% plot gate-dVdt, small fluxes
load('precip.mat')
for i=1:146
   vfluxsS(i)=nansum(nansum(fwflux(:,:,i).*inwagSareas1(:,:,i)./1000)); %PmE kg/m^2/s, so for volume flux we divide by 1000kg/m^3 and multiply by RAC (m^2)
   vfluxs263(i)=nansum(nansum(fwflux(:,:,i).*inwag263areas(:,:,i)./1000));
end
clear fwflux
%%
totexplained=errGateTot(1:145)+abs(dVdt11.*2e-3)+abs(volCrossMan(1:145)).'+abs(dvdtStepTd(1:145)).'+abs(vfluxsS(1:145))+abs(dvdtDiff275in27(1:145)).';
gateFluxTdSet=[gateFluxTot(1:145);gateFluxToti(2:146);gateFluxTotFeb(1:145)];
dVdtSetTd=[dVdt11;dVdt1s1];
totexplained1vi=range(gateFluxTdSet)+(max(abs(dVdtSetTd)).*2e-3)+abs(volCrossMan(1:145)).'+abs(dvdtStepTd(1:145)).'+abs(vfluxsS(1:145))+abs(dvdtDiff275in27(1:145)).';
errgdv=max(range(gateFluxTdSet),errGateTot(1:145))+max(range(dVdtSetTd),(max(abs(dVdtSetTd)).*2e-3));
figure; plot(totexplained); hold all; plot(totexplained1vi); plot(errgdv)
errgate=max(range(gateFluxTdSet),errGateTot(1:145));
errvol=max(range(dVdtSetTd),(max(abs(dVdtSetTd)).*2e-3));
errsmalls=abs(0.02*volCrossMan(1:145)).'+abs(0.02*dvdtStepTd(1:145)).'+abs(0.02*vfluxsS(1:145))+abs(0.02*dvdtDiff275in27(1:145)).';
err10=abs(0.1*dVdt1s1)+abs(0.1*gateFluxTotFeb(1:145))+abs(0.1*volCrossMan(1:145)).'+abs(0.1*dvdtStepTd(1:145)).'+abs(0.1*vfluxsS(1:145))+abs(0.1*dvdtDiff275in27(1:145)).';

figure; errorbar(gateFluxTotFeb(1:145),errgate); hold all; errorbar(dVdt1s1,errvol); errorbar(smalltermssum,errsmalls)
figure; plot(gateFluxTotFeb(1:145)-dVdt1s1+smalltermssum,errgate+errvol+errsmalls)
%%
figure; plot(gateFluxTot(1:145)-dVdt11,'linewidth',2); hold all; plot(volCrossMan); plot(dvdtDiff275in27); plot(dvdtStepTd); plot(vfluxsS)
plot(errGateTot(1:145)+abs(dVdt11.*2e-3),'k--'); plot(-errGateTot(1:145)-abs(dVdt11.*2e-3),'k--')
legend('gate-dV/dt','cross manifold','bottom diffusion','step diffusion','precip-evap','uncertainty in gate and dVdt')
title('WAG Volume Budget Top-down')
figure; plot(gateFluxTot3day(1:145)-dVdt3day1,'linewidth',2); hold all; plot(volCrossMan); plot(dvdtDiff275in27); plot(dvdtStepTd); plot(vfluxsS)
plot(errGateTot(1:145)+abs(dVdt11.*2e-3),'k--'); plot(-errGateTot(1:145)-abs(dVdt11.*2e-3),'k--')
legend('3day gate-dV/dt','cross manifold','bottom diffusion','step diffusion','precip-evap','uncertainty in gate and dVdt')
title('WAG Volume Budget Top-down')
figure; plot(gateFluxTot3day(1:145)-dVdt3day1,'linewidth',2); hold all; plot(volCrossMan); plot(dvdtDiff275in27); plot(dvdtStepTd); plot(vfluxsS)
%totexplained=errGateTot(1:145)+abs(dVdt11.*2e-3)+abs(volCrossMan(1:145)).'+abs(dvdtStepTd(1:145)).'+abs(vfluxsS(1:145))+abs(dvdtDiff275in27(1:145)).';
plot(totexplained,'k--'); plot(-totexplained,'k--')
legend('3day gate-dV/dt','cross manifold','bottom diffusion','step diffusion','precip-evap','uncertainty in gate and dVdt plus magnitude of other terms')
title('WAG Volume Budget Top-down')
%%
figure; plot(gateFluxTot(1:145)-dVdt11,'linewidth',2); hold all;
plot(gateFluxTot3day(1:145)-dVdt3day1,'--','linewidth',2);
plot(gateFluxToti(2:146)-dVdt1s1,'linewidth',2);
plot(gateFluxTot3dayi(2:146)-dVdt3day1s,'--','linewidth',2);
plot(gateFluxTotFeb(1:145)-dVdt1s1,'linewidth',2);
plot(gateFluxTot3dayF(1:145)-dVdt3day1s,'--','linewidth',2);
%plot(gateFluxTot2(1:145)-dVdt12,'linewidth',2);
%plot(gateFluxTot3day2(1:145)-dVdt3day2,'--','linewidth',2);
%plot(gateFluxTot2i(2:146)-dVdt2s1,'linewidth',2);
%plot(gateFluxTot3day2i(2:146)-dVdt3day2s,'--','linewidth',2);
plot(totexplained1vi,'k:'); plot(-totexplained1vi,'k:')
title('Top-down gate-dV/dt, solid raw, dashed 3-day running mean, dotted is explainable')
%%
smalltermssum=volCrossMan(1:145).'+dvdtStepTd(1:145).'+vfluxsS(1:145)+dvdtDiff275in27(1:145).';
possibleSums=[gateFluxTotFeb(1:145)-dVdt1s1+smalltermssum;gateFluxTotFeb(1:145)-dVdt11+smalltermssum;gateFluxTot(1:145)-dVdt1s1+smalltermssum;gateFluxTot(1:145)-dVdt11+smalltermssum;...
            gateFluxToti(2:146)-dVdt1s1+smalltermssum;gateFluxToti(2:146)-dVdt11+smalltermssum];
%%
figure; plot(-dVdt1s1,'linewidth',2); hold all; plot(gateFluxTotFeb,'linewidth',2)
plot(volCrossMan,'linewidth',2); plot(dvdtStepTd,'linewidth',2); plot(dvdtDiff275in27,'linewidth',2); plot(vfluxsS,'linewidth',2)
%plot(1:145,possibleSums(1,:),'--','linewidth',2)
plot(notclosed1,'--','linewidth',2)
plot(abs(dVdt1s1)*0.2+0.1*abs(gateFluxTotFeb(1:145)),'k:')
plot(-abs(dVdt1s1)*0.2-0.1*abs(gateFluxTotFeb(1:145)),'k:')
legend('-dV/dt','gate','cross-manifold','diff at steps','diff at bottom','precip-evap','total','error')
%% 
figure; plot(possibleSums(1,:)); hold on; plot(possibleSums(1,:)+errgate+errvol+50*errsmalls,'k--'); plot(possibleSums(1,:)-errgate-errvol-50*errsmalls,'k--')
%%
figure; plot(dVdt1s1); hold all; plot(gateFluxTotFeb(1:145)+smalltermssum);
fluxTotF=gateFluxTotFeb(1:145)+smalltermssum;
fluxTot3day=[0 fluxTotF(1:end-2)+fluxTotF(2:end-1)+fluxTotF(3:end) 0]./3;
%%
%logibad=(abs(gateFluxTot(1:145)-dVdt11)>1.8e6)&(abs(gateFluxTotFeb(1:145)-dVdt1s1)>1.8e6);
%logigood=(abs(gateFluxTot(1:145)-dVdt11)<0.9*totexplained)&(abs(gateFluxTotFeb(1:145)-dVdt1s1)<0.9*totexplained);
%load('wagAreaAndFluxGateS2.mat','lat*','lon*')
load('surfaceClosedNF')
lonS=lonWagClosed;
latS=latWagClosed;
lonS(latS==0)=nan;
latS(latS==0)=nan;
%load('wagAreaAndFluxGate263v2.mat','*Closed')
load('sigma263ClosedNF')
lon263=lonWagClosed;
lat263=latWagClosed;
lon263(lat263==0)=nan;
%load('wagAreaAndFluxGate265v2.mat','*Closed')
load('sigma265ClosedNF')
lon265=lonWagClosed;
lat265=latWagClosed;
lon265(lat265==0)=nan;
%load('wagAreaAndFluxGate2675v2.mat','*Closed')
load('sigma2675ClosedNF')
lon2675=lonWagClosed;
lat2675=latWagClosed;
lon2675(lat2675==0)=nan;
%load('wagAreaAndFluxGate27v2.mat','*Closed')
load('sigma27ClosedNF')
lon27=lonWagClosed;
lat27=latWagClosed;
lon27(lat27==0)=nan;
%load('wagAreaAndFluxGate275v2.mat','*Closed')
load('sigma27ClosedNF')
lon275=lonWagClosed;
lat275=latWagClosed;
lon275(lat275==0)=nan;

for i=121:128
   %if logibad(i)
       figure; plot(lonCoast,latCoast,'k'); hold all
       plot(lon275(i,:),lat275(i,:))
       plot(lon27(i,:),lat27(i,:))
       plot(lon2675(i,:),lat2675(i,:))
       plot(lon265(i,:),lat265(i,:))
       plot(lon263(i,:),lat263(i,:))
       plot(lonS(i,:),latS(i,:))
       t1=num2str(i);
   %    t1=strcat('bad, i=',num2str(i));
       title(t1)
       axis([-5.5 -2 35 37])
       legend('coast','27.5','27','26.75','26.5','26.3','S')
%    elseif logigood(i)
%        figure; plot(lonCoast,latCoast,'k'); hold all
%        plot(lon275(i,:),lat275(i,:))
%        plot(lon27(i,:),lat27(i,:))
%        plot(lon2675(i,:),lat2675(i,:))
%        plot(lon265(i,:),lat265(i,:))
%        plot(lon263(i,:),lat263(i,:))
%        plot(lonS(i,:),latS(i,:))
%        t1=strcat('good, i=',num2str(i));
%        title(t1)
%        axis([-5.5 -2 35 37])
%        legend('coast','27.5','27','26.75','26.5','26.3','S')
%    end 
end
%%
inWag3=zeros([700 200 20 146]);
inWag3portion=inWag3;
load('saTCtCpSigma162NF.mat','Sigma')

Sigma=Sigma(:,:,1:20,1:146);
k=1;
inWag3(:,:,k,:)=sign(inWAGSpoints(:,:,1:146));
inWag3portion(:,:,k,:)=inWAGSpoints(:,:,1:146)./5;
for k=2:20
inWag3(:,:,k,1:146)=double(squeeze(Sigma(:,:,k,:))<26.3).*double(squeeze(Sigma(:,:,k,:))>1).*double(inWAGS(:,:,1:146))...
    +double(squeeze(Sigma(:,:,k,:))<26.5).*double(squeeze(Sigma(:,:,k,:))>=26.3).*double(inWAG263(:,:,1:146))...
    +double(squeeze(Sigma(:,:,k,:))<26.75).*double(squeeze(Sigma(:,:,k,:))>=26.5).*double(inWAG265(:,:,1:146))...
    +double(squeeze(Sigma(:,:,k,:))<27).*double(squeeze(Sigma(:,:,k,:))>=26.75).*double(inWAG2675(:,:,1:146))...
    +double(squeeze(Sigma(:,:,k,:))<27.5).*double(squeeze(Sigma(:,:,k,:))>=27).*double(inWAG27(:,:,1:146));
inWag3portion(:,:,k,1:146)=double(squeeze(Sigma(:,:,k,:))<26.3).*double(squeeze(Sigma(:,:,k,:))>1).*(inWAGSpoints(:,:,1:146)./5)...
    +double(squeeze(Sigma(:,:,k,:))<26.5).*double(squeeze(Sigma(:,:,k,:))>=26.3).*(inWAG263points(:,:,1:146)./5)...
    +double(squeeze(Sigma(:,:,k,:))<26.75).*double(squeeze(Sigma(:,:,k,:))>=26.5).*(inWAG265points(:,:,1:146)./5)...
    +double(squeeze(Sigma(:,:,k,:))<27).*double(squeeze(Sigma(:,:,k,:))>=26.75).*(inWAG2675points(:,:,1:146)./5)...
    +double(squeeze(Sigma(:,:,k,:))<27.5).*double(squeeze(Sigma(:,:,k,:))>=27).*(inWAG27points(:,:,1:146)./5);
end
save('lagrangeWAGboundaryPoints.mat','inW*','-v7.3')
%
load('lagrangeWAGboundaryNF2.mat')
inWag3filled=zeros([700 200 20 146]);
inWAGSf=inWAGSpoints;
inWAGSf(:,:,[38:42 51 52 80 81 88 109 110 145])=inWAG263points(:,:,[38:42 51 52 80 81 88 109 110 145]);
inWAGSf(:,:,[82:85 87 106:108])=inWAG265points(:,:,[82:85 87 106:108]);
inWAG263f=inWAG263points;
inWAG263f(:,:,[60:63 71])=inWAGSpoints(:,:,[60:63 71]);
inWAG263f(:,:,[82:87 106:108])=inWAG265points(:,:,[82:87 106:108]);
inWAG265f=inWAG265points;
inWAG265f(:,:,[60:63 71])=inWAGSpoints(:,:,[60:63 71]);
inWAG265f(:,:,[64 66:70 103])=inWAG263points(:,:,[64 66:70 103]);
inWAG2675f=inWAG2675points;
inWAG2675f(:,:,62:63)=inWAGSpoints(:,:,62:63);
inWAG2675f(:,:,67:70)=inWAG263points(:,:,67:70);
inWAG2675f(:,:,[15 43 48 51:56 72 112])=inWAG265points(:,:,[15 43 48 51:56 72 112]);
inWAG2675f(:,:,[61 64])=inWAG27points(:,:,[61 64]);
inWAG27f=inWAG27points;
inWAG27f(:,:,62:63)=inWAGSpoints(:,:,62:63);
inWAG27f(:,:,67:70)=inWAG263points(:,:,67:70);
inWAG27f(:,:,[71 117])=inWAG2675points(:,:,[71 117]);
inWAG275f=inWAG275points;
inWAG275f(:,:,67:70)=inWAG263points(:,:,67:70);
inWAG275f(:,:,[71 117])=inWAG2675points(:,:,[71 117]);
inWAG275f(:,:,[39 56 66 81 118 121 146])=inWAG27points(:,:,[39 56 66 81 118 121 146]);
%
load('saTCtCpSigma162NF.mat','Sigma')
load('lagrangeWAGboundaryNF2filled.mat')
load('tsSigmaSnapshots162NF.mat','Sigma')
Sigma=Sigma(:,:,1:20,1:146);
%
k=1;
inWag3portionFT=zeros([700 200 20 146]);
inWag3portionFT(:,:,k,1:146)=inWAGSf;
inWag3portionFB=zeros([700 200 20 146]);
inWag3portionFB(:,:,k,1:146)=inWAG263f;
for k=2:20%through 16 is enough!
inWag3snap2(:,:,k,1:145)=double(squeeze(Sigma(:,:,k,:))<26.3).*double(squeeze(Sigma(:,:,k,:))>1).*inWAGSf(:,:,2:146)...
    +double(squeeze(Sigma(:,:,k,:))<26.5).*double(squeeze(Sigma(:,:,k,:))>=26.3).*inWAG263f(:,:,2:146)...
    +double(squeeze(Sigma(:,:,k,:))<26.75).*double(squeeze(Sigma(:,:,k,:))>=26.5).*inWAG265f(:,:,2:146)...
    +double(squeeze(Sigma(:,:,k,:))<27).*double(squeeze(Sigma(:,:,k,:))>=26.75).*inWAG2675f(:,:,2:146)...
    +double(squeeze(Sigma(:,:,k,:))<27.5).*double(squeeze(Sigma(:,:,k,:))>=27).*inWAG27f(:,:,2:146);
inWag3portionFT(:,:,k,1:146)=double(squeeze(Sigma(:,:,k,:))<26.3).*double(squeeze(Sigma(:,:,k,:))>1).*(inWAGSf(:,:,1:146)./5)...
    +double(squeeze(Sigma(:,:,k,:))<26.5).*double(squeeze(Sigma(:,:,k,:))>=26.3).*(inWAG263f(:,:,1:146)./5)...
    +double(squeeze(Sigma(:,:,k,:))<26.75).*double(squeeze(Sigma(:,:,k,:))>=26.5).*(inWAG265f(:,:,1:146)./5)...
    +double(squeeze(Sigma(:,:,k,:))<27).*double(squeeze(Sigma(:,:,k,:))>=26.75).*(inWAG2675f(:,:,1:146)./5)...
    +double(squeeze(Sigma(:,:,k,:))<27.5).*double(squeeze(Sigma(:,:,k,:))>=27).*(inWAG27f(:,:,1:146)./5);
inWag3portionFB(:,:,k,1:146)=double(squeeze(Sigma(:,:,k,:))<26.3).*double(squeeze(Sigma(:,:,k,:))>1).*(inWAG263f(:,:,1:146)./5)...
    +double(squeeze(Sigma(:,:,k,:))<26.5).*double(squeeze(Sigma(:,:,k,:))>=26.3).*(inWAG265f(:,:,1:146)./5)...
    +double(squeeze(Sigma(:,:,k,:))<26.75).*double(squeeze(Sigma(:,:,k,:))>=26.5).*(inWAG2675f(:,:,1:146)./5)...
    +double(squeeze(Sigma(:,:,k,:))<27).*double(squeeze(Sigma(:,:,k,:))>=26.75).*(inWAG27f(:,:,1:146)./5)...
    +double(squeeze(Sigma(:,:,k,:))<27.5).*double(squeeze(Sigma(:,:,k,:))>=27).*(inWAG275f(:,:,1:146)./5);
end
save('lagrangeWAGboundaryPortionFilledTB.mat','inWag3p*','-v7.3')
clear Sigma
%
volPortionT=squeeze(nansum(nansum(squeeze(inWag3portionFT(:,:,1,:)).*SSHa(:,:,1:146).*repmat(RAC,[1 1 146]))));
volPortionT=volPortionT+squeeze(nansum(nansum(nansum(repmat(cellVol(:,:,1:20),[1 1 1 146]).*inWag3portionFT))));
volPortionB=squeeze(nansum(nansum(squeeze(inWag3portionFB(:,:,1,:)).*SSHa(:,:,1:146).*repmat(RAC,[1 1 146]))));
volPortionB=volPortionB+squeeze(nansum(nansum(nansum(repmat(cellVol(:,:,1:20),[1 1 1 146]).*inWag3portionFB))));
errVol(14:146)=abs(volPortionB(14:146)-volPortionT(14:146));
dVdt=diff(0.5*volPortionT+0.5*volPortionB)./86400;
dVdt(13)=0;
%
for k=1:20
inWag3S(:,:,k,1:146)=double(squeeze(Sigma(:,:,k,:))<27.5).*double(squeeze(Sigma(:,:,k,:))>1).*(inWAGSf(:,:,1:146)./5);
inWag3263(:,:,k,1:146)=double(squeeze(Sigma(:,:,k,:))<27.5).*double(squeeze(Sigma(:,:,k,:))>1).*(inWAG263f(:,:,1:146)./5);
inWag3265(:,:,k,1:146)=double(squeeze(Sigma(:,:,k,:))<27.5).*double(squeeze(Sigma(:,:,k,:))>1).*(inWAG265f(:,:,1:146)./5);
inWag32675(:,:,k,1:146)=double(squeeze(Sigma(:,:,k,:))<27.5).*double(squeeze(Sigma(:,:,k,:))>1).*(inWAG2675f(:,:,1:146)./5);
inWag327(:,:,k,1:146)=double(squeeze(Sigma(:,:,k,:))<27.5).*double(squeeze(Sigma(:,:,k,:))>1).*(inWAG27f(:,:,1:146)./5);
inWag3275(:,:,k,1:146)=double(squeeze(Sigma(:,:,k,:))<27.5).*double(squeeze(Sigma(:,:,k,:))>1).*(inWAG275f(:,:,1:146)./5);
end
save('lagrangeWAGboundariesSingleLayer.mat','inWag3*','-v7.3')
('lagrangeWAGvolumesLayers.mat')
load('uvwSSHDailyDepth1rotated148F.mat', 'SSHa')
volS=squeeze(nansum(nansum(squeeze(inWag3S(:,:,1,:)).*SSHa(:,:,1:146).*repmat(RAC,[1 1 146]))))+squeeze(nansum(nansum(nansum(repmat(cellVol,[1 1 1 146]).*inWag3S))));
vol263=squeeze(nansum(nansum(squeeze(inWag3263(:,:,1,:)).*SSHa(:,:,1:146).*repmat(RAC,[1 1 146]))))+squeeze(nansum(nansum(nansum(repmat(cellVol,[1 1 1 146]).*inWag3263))));
vol265=squeeze(nansum(nansum(squeeze(inWag3265(:,:,1,:)).*SSHa(:,:,1:146).*repmat(RAC,[1 1 146]))))+squeeze(nansum(nansum(nansum(repmat(cellVol,[1 1 1 146]).*inWag3265))));
vol2675=squeeze(nansum(nansum(squeeze(inWag32675(:,:,1,:)).*SSHa(:,:,1:146).*repmat(RAC,[1 1 146]))))+squeeze(nansum(nansum(nansum(repmat(cellVol,[1 1 1 146]).*inWag32675))));
vol27=squeeze(nansum(nansum(squeeze(inWag327(:,:,1,:)).*SSHa(:,:,1:146).*repmat(RAC,[1 1 146]))))+squeeze(nansum(nansum(nansum(repmat(cellVol,[1 1 1 146]).*inWag327))));
vol275=squeeze(nansum(nansum(squeeze(inWag3275(:,:,1,:)).*SSHa(:,:,1:146).*repmat(RAC,[1 1 146]))))+squeeze(nansum(nansum(nansum(repmat(cellVol,[1 1 1 146]).*inWag3275))));

% better lagrangeWAGboundary
inWag3b=zeros([700 200 20 146]);
for i=14:146
   if squeeze(sum(sum(inWAGS(:,:,i))))==0
       if squeeze(sum(sum(inWAG263(:,:,i))))>0
            inWAGS(:,:,i)=inWAG263(:,:,i);
       elseif squeeze(sum(sum(inWAG265(:,:,i))))>0
            inWAGS(:,:,i)=inWAG265(:,:,i);
            inWAG263(:,:,i)=inWAG265(:,:,i);
       end
   end
   
   if squeeze(sum(sum(inWAG263(:,:,i))))==0
       if squeeze(sum(sum(inWAGS(:,:,i))))>0
            inWAG263(:,:,i)=inWAGS(:,:,i);
       elseif squeeze(sum(sum(inWAG265(:,:,i))))>0
            inWAGS(:,:,i)=inWAG265(:,:,i);
            inWAG263(:,:,i)=inWAG265(:,:,i);
       end
   end
   
   if squeeze(sum(sum(inWAG265(:,:,i))))==0
       if squeeze(sum(sum(inWAG263(:,:,i))))>0
            inWAG265(:,:,i)=inWAG263(:,:,i);
       elseif squeeze(sum(sum(inWAG2675(:,:,i))))>0
            inWAG265(:,:,i)=inWAG2675(:,:,i);
            inWAG263(:,:,i)=inWAG2675(:,:,i);
       end
   end
   
   if squeeze(sum(sum(inWAG2675(:,:,i))))==0
       if squeeze(sum(sum(inWAG265(:,:,i))))>0
            inWAG2675(:,:,i)=inWAG265(:,:,i);
       elseif squeeze(sum(sum(inWAG27(:,:,i))))>0
            inWAG265(:,:,i)=inWAG27(:,:,i);
            inWAG2675(:,:,i)=inWAG27(:,:,i);
       end
   end
   
   if squeeze(sum(sum(inWAG27(:,:,i))))==0
       if squeeze(sum(sum(inWAG2675(:,:,i))))>0
            inWAG27(:,:,i)=inWAG2675(:,:,i);
       elseif squeeze(sum(sum(inWAG275(:,:,i))))>0
            inWAG2675(:,:,i)=inWAG275(:,:,i);
            inWAG27(:,:,i)=inWAG275(:,:,i);
       end
   end
   
   if squeeze(sum(sum(inWAG275(:,:,i))))==0
       if squeeze(sum(sum(inWAG27(:,:,i))))>0
            inWAG275(:,:,i)=inWAG27(:,:,i);
       end
   end

end
for k=1:20%through 16 is enough!
inWag3(:,:,k,1:146)=double(squeeze(Sigma(:,:,k,:))<26.3).*double(squeeze(Sigma(:,:,k,:))>1).*inWAGS(:,:,1:146)...
    +double(squeeze(Sigma(:,:,k,:))<26.5).*double(squeeze(Sigma(:,:,k,:))>=26.3).*inWAG263(:,:,1:146)...
    +double(squeeze(Sigma(:,:,k,:))<26.75).*double(squeeze(Sigma(:,:,k,:))>=26.5).*inWAG265(:,:,1:146)...
    +double(squeeze(Sigma(:,:,k,:))<27).*double(squeeze(Sigma(:,:,k,:))>=26.75).*inWAG2675(:,:,1:146)...
    +double(squeeze(Sigma(:,:,k,:))<27.5).*double(squeeze(Sigma(:,:,k,:))>=27).*inWAG27(:,:,1:146);%...
    %+double(squeeze(Sigma(:,:,k,:))<28).*double(squeeze(Sigma(:,:,k,:))>=27.5).*inWAG275(:,:,1:131);
end
save('lagrangeWAGboundaryNF3.mat','inW*','-v7.3')
