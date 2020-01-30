%% explanation
%% load for nice sections plots
section1=17:23; section2=30:36; section3=40:50; section4=88:98; section5=134:142; section3issues=[43 48];
load('gateAdvectionIndexedMar.mat', 'gateFluxTotMar')
load('lagrangeVolBudgetWorkingMar1.mat', 'dvdtDiff275in27', 'vfluxsS', 'volCrossMan', 'dVdt1s1','dvdtStepTd')
%%
gateFlux=gateFluxMarC(2:end);%gateFluxTotMar(2:end);%-4e5*ones(size(gateFluxTotMar(2:end)));
dVdt=dVdt3;%dVdt1s1;%+3e5*ones(size(dVdt1s1));
sumtermserrs=-dVdt(1:145)+gateFlux(1:145)+dvdtDiff275in27(1:145).'+vfluxsS(1:145)+volCrossMan(1:145).'+dvdtStepTd(1:145).';
sumterms=-dVdt(1:145)+gateFlux(1:145)+dvdtDiff275in27(1:145).'+vfluxsS(1:145);
sumtermsold=-dVdt1s1(1:145)+gateFluxTotMar(2:146)+dvdtDiff275in27(1:145).'+vfluxsS(1:145);
figure; plot(section1,-dVdt(section1),'b','linewidth',2);hold all; plot(section1,gateFlux(section1),'r','linewidth',2);
plot(section1,dvdtDiff275in27(section1),'m');
plot(section1,vfluxsS(section1),'c');
plot(section1,sumterms(section1),'k','linewidth',2);
plot(section1,volCrossMan(section1),'Color',[0 0.5 0]);
plot(section1,dvdtStepTd(section1),'Color',[0.75 0 0.75]);
legend('-dV/dt','gate','diffusion through bottom','surface P-E','sum','error cross-manifold','error diffusion steps')
plot(section2,-dVdt(section2),'b','linewidth',2);plot(section3,-dVdt(section3),'b','linewidth',2);plot(section4,-dVdt(section4),'b','linewidth',2);plot(section5,-dVdt(section5),'b','linewidth',2); 
plot(section2,gateFlux(section2),'r','linewidth',2);plot(section3,gateFlux(section3),'r','linewidth',2);plot(section4,gateFlux(section4),'r','linewidth',2);plot(section5,gateFlux(section5),'r','linewidth',2);
plot(section2,dvdtDiff275in27(section2),'m');plot(section3,dvdtDiff275in27(section3),'m');plot(section4,dvdtDiff275in27(section4),'m');plot(section5,dvdtDiff275in27(section5),'m');
plot(section2,vfluxsS(section2),'c');plot(section3,vfluxsS(section3),'c');plot(section4,vfluxsS(section4),'c');plot(section5,vfluxsS(section5),'c');
plot(section2,volCrossMan(section2),'Color',[0 0.5 0]);plot(section3,volCrossMan(section3),'Color',[0 0.5 0]);plot(section4,volCrossMan(section4),'Color',[0 0.5 0]);plot(section5,volCrossMan(section5),'Color',[0 0.5 0]);
plot(section2,dvdtStepTd(section2),'Color',[0.75 0 0.75]);plot(section3,dvdtStepTd(section3),'Color',[0.75 0 0.75]);plot(section4,dvdtStepTd(section4),'Color',[0.75 0 0.75]);plot(section5,dvdtStepTd(section5),'Color',[0.75 0 0.75]);
plot(section2,sumterms(section2),'k','linewidth',2);plot(section3,sumterms(section3),'k','linewidth',2);plot(section4,sumterms(section4),'k','linewidth',2);plot(section5,sumterms(section5),'k','linewidth',2);
%% location
%load('wagAreaAndFluxGate.mat', 'inWAG')%adjust by adding 26, 265, or 27 to filename
load('geometrySpinupSteady.mat','XC','YC','XG','YG')
load('wagAreaAndFluxGate275v2.mat', '*Closed')
figure;
for i=1:146
    c=[];
    if latWagClosed(i,1)>0
        np=find(latWagClosed(i,:)>0,1,'last');
        inWAG=inpolygon(XC,YC,lonWagClosed(i,1:np),latWagClosed(i,1:np));
        %inWAGll=inpolygon(XG(1:end-1,1:end-1),YG(1:end-1,1:end-1),lonWagClosed(i,1:np),latWagClosed(i,1:np));
        %inWAGlu=inpolygon(XG(1:end-1,2:end),YG(1:end-1,2:end),lonWagClosed(i,1:np),latWagClosed(i,1:np));
        %inWAGrl=inpolygon(XG(2:end,1:end-1),YG(2:end,1:end-1),lonWagClosed(i,1:np),latWagClosed(i,1:np));
        %inWAGru=inpolygon(XG(2:end,2:end),YG(2:end,2:end),lonWagClosed(i,1:np),latWagClosed(i,1:np));
        [c,h]=contour(XC,YC,double(inWAG),[1 1]);
    end
    if ~isempty(c)
    npointsi=find(c(2,:)==floor(c(2,:)));
    varhold=max(c(2,npointsi));
    i1=1+find(c(2,:)==varhold);
    iend=varhold+i1-1;
    inWAG275(:,:,i)=inpolygon(XC,YC,c(1,i1:iend),c(2,i1:iend));
    inWAGll=inpolygon(XG(1:end-1,1:end-1),YG(1:end-1,1:end-1),c(1,i1:iend),c(2,i1:iend));
    inWAGlu=inpolygon(XG(1:end-1,2:end),YG(1:end-1,2:end),c(1,i1:iend),c(2,i1:iend));
    inWAGrl=inpolygon(XG(2:end,1:end-1),YG(2:end,1:end-1),c(1,i1:iend),c(2,i1:iend));
    inWAGru=inpolygon(XG(2:end,2:end),YG(2:end,2:end),c(1,i1:iend),c(2,i1:iend));
    inWAG275corners2(:,:,i)=double(inWAGll)+double(inWAGlu)+double(inWAGrl)+double(inWAGru);
    end
end
clear inWAG inWAGl* inWAGr*
%%
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
%%
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
%%
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
%%
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
%%
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
%%
inWag3=zeros([700 200 20 146]);
load('saTCtCpSigma162NF.mat','Sigma')

Sigma=Sigma(:,:,1:20,1:146);
k=1;
inWag3(:,:,k,:)=inWAGS(:,:,1:146);
for k=2:20%through 16 is enough!
inWag3(:,:,k,1:146)=double(squeeze(Sigma(:,:,k,:))<26.3).*double(squeeze(Sigma(:,:,k,:))>1).*inWAGS(:,:,1:146)...
    +double(squeeze(Sigma(:,:,k,:))<26.5).*double(squeeze(Sigma(:,:,k,:))>=26.3).*inWAG263(:,:,1:146)...
    +double(squeeze(Sigma(:,:,k,:))<26.75).*double(squeeze(Sigma(:,:,k,:))>=26.5).*inWAG265(:,:,1:146)...
    +double(squeeze(Sigma(:,:,k,:))<27).*double(squeeze(Sigma(:,:,k,:))>=26.75).*inWAG2675(:,:,1:146)...
    +double(squeeze(Sigma(:,:,k,:))<27.5).*double(squeeze(Sigma(:,:,k,:))>=27).*inWAG27(:,:,1:146);%...
    %+double(squeeze(Sigma(:,:,k,:))<28).*double(squeeze(Sigma(:,:,k,:))>=27.5).*inWAG275(:,:,1:131);
end
%save('lagrangeWAGboundaryNF2.mat','inW*','-v7.3')
%%
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
%%
%load('lagrangeWAGboundaryNF2.mat')
%inWag3filled=zeros([700 200 20 146]);
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
%%
load('saTCtCpSigma162NF.mat','Sigma')
%load('lagrangeWAGboundaryNF2filled.mat')
%load('tsSigmaSnapshots162NF.mat','Sigma')
Sigma=Sigma(:,:,1:20,1:146);
%%
k=1;
inWag3portionFT=zeros([700 200 20 146]);
inWag3portionFT(:,:,k,1:146)=inWAGSf;
inWag3portionFB=zeros([700 200 20 146]);
inWag3portionFB(:,:,k,1:146)=inWAG263f;
for k=2:20%through 16 is enough!
% inWag3snap2(:,:,k,1:145)=double(squeeze(Sigma(:,:,k,:))<26.3).*double(squeeze(Sigma(:,:,k,:))>1).*inWAGSf(:,:,2:146)...
%     +double(squeeze(Sigma(:,:,k,:))<26.5).*double(squeeze(Sigma(:,:,k,:))>=26.3).*inWAG263f(:,:,2:146)...
%     +double(squeeze(Sigma(:,:,k,:))<26.75).*double(squeeze(Sigma(:,:,k,:))>=26.5).*inWAG265f(:,:,2:146)...
%     +double(squeeze(Sigma(:,:,k,:))<27).*double(squeeze(Sigma(:,:,k,:))>=26.75).*inWAG2675f(:,:,2:146)...
%     +double(squeeze(Sigma(:,:,k,:))<27.5).*double(squeeze(Sigma(:,:,k,:))>=27).*inWAG27f(:,:,2:146);
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
%%
volPortionT=squeeze(nansum(nansum(squeeze(inWag3portionFT(:,:,1,:)).*SSHa(:,:,1:146).*repmat(RAC,[1 1 146]))));
volPortionT=volPortionT+squeeze(nansum(nansum(nansum(repmat(cellVol(:,:,1:20),[1 1 1 146]).*inWag3portionFT))));
volPortionB=squeeze(nansum(nansum(squeeze(inWag3portionFB(:,:,1,:)).*SSHa(:,:,1:146).*repmat(RAC,[1 1 146]))));
volPortionB=volPortionB+squeeze(nansum(nansum(nansum(repmat(cellVol(:,:,1:20),[1 1 1 146]).*inWag3portionFB))));
errVol(14:146)=abs(volPortionB(14:146)-volPortionT(14:146));
dVdt=diff(0.5*volPortionT+0.5*volPortionB)./86400;
dVdt(13)=0;
%%
for k=1:20
inWag3S(:,:,k,1:146)=double(squeeze(Sigma(:,:,k,:))<27.5).*double(squeeze(Sigma(:,:,k,:))>1).*(inWAGSf(:,:,1:146)./5);
inWag3263(:,:,k,1:146)=double(squeeze(Sigma(:,:,k,:))<27.5).*double(squeeze(Sigma(:,:,k,:))>1).*(inWAG263f(:,:,1:146)./5);
inWag3265(:,:,k,1:146)=double(squeeze(Sigma(:,:,k,:))<27.5).*double(squeeze(Sigma(:,:,k,:))>1).*(inWAG265f(:,:,1:146)./5);
inWag32675(:,:,k,1:146)=double(squeeze(Sigma(:,:,k,:))<27.5).*double(squeeze(Sigma(:,:,k,:))>1).*(inWAG2675f(:,:,1:146)./5);
inWag327(:,:,k,1:146)=double(squeeze(Sigma(:,:,k,:))<27.5).*double(squeeze(Sigma(:,:,k,:))>1).*(inWAG27f(:,:,1:146)./5);
inWag3275(:,:,k,1:146)=double(squeeze(Sigma(:,:,k,:))<27.5).*double(squeeze(Sigma(:,:,k,:))>1).*(inWAG275f(:,:,1:146)./5);
end
%save('lagrangeWAGboundariesSingleLayer.mat','inWag3*','-v7.3')
%('lagrangeWAGvolumesLayers.mat')
%load('uvwSSHDailyDepth1rotated148F.mat', 'SSHa')
volS=squeeze(nansum(nansum(squeeze(inWag3S(:,:,1,:)).*SSHa(:,:,1:146).*repmat(RAC,[1 1 146]))))+squeeze(nansum(nansum(nansum(repmat(cellVol,[1 1 1 146]).*inWag3S))));
vol263=squeeze(nansum(nansum(squeeze(inWag3263(:,:,1,:)).*SSHa(:,:,1:146).*repmat(RAC,[1 1 146]))))+squeeze(nansum(nansum(nansum(repmat(cellVol,[1 1 1 146]).*inWag3263))));
vol265=squeeze(nansum(nansum(squeeze(inWag3265(:,:,1,:)).*SSHa(:,:,1:146).*repmat(RAC,[1 1 146]))))+squeeze(nansum(nansum(nansum(repmat(cellVol,[1 1 1 146]).*inWag3265))));
vol2675=squeeze(nansum(nansum(squeeze(inWag32675(:,:,1,:)).*SSHa(:,:,1:146).*repmat(RAC,[1 1 146]))))+squeeze(nansum(nansum(nansum(repmat(cellVol,[1 1 1 146]).*inWag32675))));
vol27=squeeze(nansum(nansum(squeeze(inWag327(:,:,1,:)).*SSHa(:,:,1:146).*repmat(RAC,[1 1 146]))))+squeeze(nansum(nansum(nansum(repmat(cellVol,[1 1 1 146]).*inWag327))));
vol275=squeeze(nansum(nansum(squeeze(inWag3275(:,:,1,:)).*SSHa(:,:,1:146).*repmat(RAC,[1 1 146]))))+squeeze(nansum(nansum(nansum(repmat(cellVol,[1 1 1 146]).*inWag3275))));

%% better lagrangeWAGboundary
% inWag3b=zeros([700 200 20 146]);
% for i=14:146
%    if squeeze(sum(sum(inWAGS(:,:,i))))==0
%        if squeeze(sum(sum(inWAG263(:,:,i))))>0
%             inWAGS(:,:,i)=inWAG263(:,:,i);
%        elseif squeeze(sum(sum(inWAG265(:,:,i))))>0
%             inWAGS(:,:,i)=inWAG265(:,:,i);
%             inWAG263(:,:,i)=inWAG265(:,:,i);
%        end
%    end
%    
%    if squeeze(sum(sum(inWAG263(:,:,i))))==0
%        if squeeze(sum(sum(inWAGS(:,:,i))))>0
%             inWAG263(:,:,i)=inWAGS(:,:,i);
%        elseif squeeze(sum(sum(inWAG265(:,:,i))))>0
%             inWAGS(:,:,i)=inWAG265(:,:,i);
%             inWAG263(:,:,i)=inWAG265(:,:,i);
%        end
%    end
%    
%    if squeeze(sum(sum(inWAG265(:,:,i))))==0
%        if squeeze(sum(sum(inWAG263(:,:,i))))>0
%             inWAG265(:,:,i)=inWAG263(:,:,i);
%        elseif squeeze(sum(sum(inWAG2675(:,:,i))))>0
%             inWAG265(:,:,i)=inWAG2675(:,:,i);
%             inWAG263(:,:,i)=inWAG2675(:,:,i);
%        end
%    end
%    
%    if squeeze(sum(sum(inWAG2675(:,:,i))))==0
%        if squeeze(sum(sum(inWAG265(:,:,i))))>0
%             inWAG2675(:,:,i)=inWAG265(:,:,i);
%        elseif squeeze(sum(sum(inWAG27(:,:,i))))>0
%             inWAG265(:,:,i)=inWAG27(:,:,i);
%             inWAG2675(:,:,i)=inWAG27(:,:,i);
%        end
%    end
%    
%    if squeeze(sum(sum(inWAG27(:,:,i))))==0
%        if squeeze(sum(sum(inWAG2675(:,:,i))))>0
%             inWAG27(:,:,i)=inWAG2675(:,:,i);
%        elseif squeeze(sum(sum(inWAG275(:,:,i))))>0
%             inWAG2675(:,:,i)=inWAG275(:,:,i);
%             inWAG27(:,:,i)=inWAG275(:,:,i);
%        end
%    end
%    
%    if squeeze(sum(sum(inWAG275(:,:,i))))==0
%        if squeeze(sum(sum(inWAG27(:,:,i))))>0
%             inWAG275(:,:,i)=inWAG27(:,:,i);
%        end
%    end
% 
% end
% for k=1:20%through 16 is enough!
% inWag3(:,:,k,1:146)=double(squeeze(Sigma(:,:,k,:))<26.3).*double(squeeze(Sigma(:,:,k,:))>1).*inWAGS(:,:,1:146)...
%     +double(squeeze(Sigma(:,:,k,:))<26.5).*double(squeeze(Sigma(:,:,k,:))>=26.3).*inWAG263(:,:,1:146)...
%     +double(squeeze(Sigma(:,:,k,:))<26.75).*double(squeeze(Sigma(:,:,k,:))>=26.5).*inWAG265(:,:,1:146)...
%     +double(squeeze(Sigma(:,:,k,:))<27).*double(squeeze(Sigma(:,:,k,:))>=26.75).*inWAG2675(:,:,1:146)...
%     +double(squeeze(Sigma(:,:,k,:))<27.5).*double(squeeze(Sigma(:,:,k,:))>=27).*inWAG27(:,:,1:146);%...
%     %+double(squeeze(Sigma(:,:,k,:))<28).*double(squeeze(Sigma(:,:,k,:))>=27.5).*inWAG275(:,:,1:131);
% end
% save('lagrangeWAGboundaryNF3.mat','inW*','-v7.3')
%% volume budget precip
%load('lagrangeWAGboundaryNF3.mat')
%load('volumeWagLagrange.mat')
%load('volumeWagLagrangeErr2.mat')
%load('oceFWflx.mat')
load('precip.mat')
%load('gateAdvection26527.mat')
%load('gateAdvection26265.mat')
%load('gateAdvection27275.mat')
%load('gateAdvectionSurface.mat')
load('distancesAreas.mat')
%load('lagrangeWAGboundary.mat','inW*')
for i=1:146
   vfluxsPF(i)=nansum(nansum(fwflux(:,:,i).*RAC.*inWag3portionF(:,:,1,i)./5000)); %PmE kg/m^2/s, so for volume flux we divide by 1000kg/m^3 and multiply by RAC (m^2)
end
clear fwflux
%% volume budget SSH storage
load('uvwSSHDailyDepth1rotated148F.mat', 'SSHa')
for i=1:146
   volSSHpf(:,:,i)=SSHa(:,:,i).*RAC.*double(inWag3portionF(:,:,1,i))./5; %m^3
end
sshFluxPF=diff(squeeze(nansum(nansum(volSSHpf))))/86400;
%% layered budgets

%surface
load('vertVolIso263.mat', 'vertVol263inWagS')
load('wagAreaAndFluxGateS2.mat','gateFluxS26')
load('geometrySpinupSteady.mat','dInterface')
dZ(1,1,1:46)=diff(dInterface);
cellVol=repmat(RAC,[1 1 46]).*hFacC.*repmat(dZ,[700 200 1]);
inWagS3(:,:,1,1:146)=inWAGS(:,:,1:146);
for k=1:20
inWagS3(:,:,k,1:146)=double(squeeze(Sigma(:,:,k,:))<26.3).*double(squeeze(Sigma(:,:,k,:))>1).*inWAGS(:,:,1:146);
end
volWagS=squeeze(sum(sum(volSSH+squeeze(sum(inWagS3.*repmat(cellVol(:,:,1:20),[1 1 1 146]),3)))));
volWagS(volWagS==0)=nan;
fluxVolWagS=diff(volWagS)/86400;
%% surface layer volume budget cont
fluxVolWagS2=fluxVolWagS;
fluxVolWagS2(fluxVolWagS2>2e6)=0;%fluxVolWagS2(fluxVolWagS2>1e6)/10;
fluxVolWagS2(fluxVolWagS2<-2e6)=0;%fluxVolWagS2(fluxVolWagS2<-1e6)/10;
figure; plot(gateFluxS26); hold all; %plot(vertVol263inWagS); 
plot(vfluxs); 
plot(-fluxVolWagS2);
plot(gateFluxS26(1:145)+vfluxs(1:145)+fluxVolWagS2.','k--','linewidth',2)
%plot(gateFluxS26(1:145)+vertVol263inWagS(1:145).'+vfluxs(1:145)-fluxVolWagS(1:145).')
title('Surface to \sigma=26.3 volume fluxes')
legend('gate','cross-\sigma=26.3','P-E','-dV/dt','total')

%save('surfaceVolLagrange.mat','volSSH','vfluxs','volWagS','gateFluxS26','inWagS3','fluxVolWagS','sshFlux')
%% overall volume budget attempt
load('wagAreaAndFluxGateS2filled.mat','gateFluxS26','yGateS','meanDS')
%load('wagAreaAndFluxGate275v2filled.mat','gateFlux27528')
load('wagAreaAndFluxGate2675v2filled.mat','gateFlux267527','yGate2675','meanD267527')
load('wagAreaAndFluxGate27v2filled.mat','gateFlux27275','yGate27','meanD27275')
load('wagAreaAndFluxGate265v2filled.mat','gateFlux2652675','yGate265','meanD2652675')
load('wagAreaAndFluxGate263v2filled.mat','gateFlux26265','yGate263','meanD26265')
gateFluxTot1=gateFlux26265(1:146)+gateFlux2652675+gateFlux267527+gateFlux27275+gateFluxS26(1:146);
errS=abs(gateFluxS26).*sqrt((1./meanDS).^2+(1./range(yGateS.')).^2);
errS(isnan(errS))=0;
err263=abs(gateFlux26265).*sqrt((1./meanD26265).^2+(1./range(yGate263.')).^2);%...
err263(isnan(err263))=0;
err265=abs(gateFlux2652675).*sqrt((1./meanD2652675).^2+(1./range(yGate265.')).^2);
err265(isnan(err265))=0;
err2675=abs(gateFlux267527).*sqrt((1./meanD267527).^2+(1./range(yGate2675.')).^2);%...
err2675(isnan(err2675))=0;
err27=abs(gateFlux27275).*sqrt((1./meanD27275).^2+(1./range(yGate27.')).^2);
err27(isnan(err27))=0;
errGateTot1=errS+err263+err265+err2675+err27;
%% updated gate volume fluxes
load('gateAdvectionS263.mat','gateFluxS263','yGate','meanDS')
yGateS=yGate;
load('gateAdvection267527.mat','gateFlux267527','yGate','meanD267527')
yGate2675=yGate;
load('gateAdvection27275.mat','gateFlux27275','yGate','meanD27275')
yGate27=yGate;
load('gateAdvection2652675.mat','gateFlux2652675','yGate','meanD2652675')
yGate265=yGate;
load('gateAdvection263265.mat','gateFlux263265','yGate','meanD263265')
yGate263=yGate;

gateFluxTot=gateFlux263265(1:146)+gateFlux2652675+gateFlux267527+gateFlux27275+gateFluxS26(1:146);
errS=abs(gateFluxS26).*sqrt((1./meanDS).^2+(1./range(yGateS.')).^2);
errS(isnan(errS))=0;
err263=abs(gateFlux263265).*sqrt((1./meanD26265).^2+(1./range(yGate263.')).^2);%...
err263(isnan(err263))=0;
err265=abs(gateFlux2652675).*sqrt((1./meanD2652675).^2+(1./range(yGate265.')).^2);
err265(isnan(err265))=0;
err2675=abs(gateFlux267527).*sqrt((1./meanD267527).^2+(1./range(yGate2675.')).^2);%...
err2675(isnan(err2675))=0;
err27=abs(gateFlux27275).*sqrt((1./meanD27275).^2+(1./range(yGate27.')).^2);
err27(isnan(err27))=0;
errGateTot=errS+err263+err265+err2675+err27;

%% overall volume budget attempt
volWagPortion=squeeze(sum(sum(volSSHpf(:,:,1:146)+squeeze(sum(inWag3portionF.*repmat(cellVol(:,:,1:20),[1 1 1 146]),3)))));
volWagPortion0=squeeze(sum(sum(volSSHp(:,:,1:146)+squeeze(sum(inWag3portion.*repmat(cellVol(:,:,1:20),[1 1 1 146]),3)))));
volWag=squeeze(sum(sum(volSSH(:,:,1:146)+squeeze(sum(inWag3.*repmat(cellVol(:,:,1:20),[1 1 1 146]),3)))));
volWagF=squeeze(sum(sum(volSSHf(:,:,1:146)+squeeze(sum(inWag3filled.*repmat(cellVol(:,:,1:20),[1 1 1 146]),3)))));

volWagSnap2=squeeze(sum(sum(volSSH(:,:,1:145)+squeeze(sum(inWag3snap2.*repmat(cellVol(:,:,1:20),[1 1 1 145]),3)))));
volWagSnap2(volWagSnap2==0)=nan;
%volWag2=0.5*volWag(1:end-1)+0.5*volWag(2:end);
fluxVolWag22=diff(volWagSnap2)/86400;
%fluxVolWag2=0.5*fluxVolWag(2:end)+0.5*fluxVolWag(1:end-1);

%% volume using each layer horizontal boundaries
% load('lagrangeWAGboundaryPortionFilled.mat', 'inWAG263f')
% load('lagrangeWAGboundaryPortionFilled.mat', 'inWAG265f')
% load('lagrangeWAGboundaryPortionFilled.mat', 'inWAG2675f')
% load('lagrangeWAGboundaryPortionFilled.mat', 'inWAG275f')
% load('lagrangeWAGboundaryPortionFilled.mat', 'inWAG27f')
% load('lagrangeWAGboundaryPortionFilled.mat', 'inWAGSf')
% load('distancesAreas.mat','RAC')
% load('distancesAreas.mat','hFacC')
% load('geometrySpinupSteady.mat','dInterface')
% dZ(1,1,1:46)=diff(dInterface);
% cellVol=repmat(RAC,[1 1 46]).*hFacC.*repmat(dZ,[700 200 1]);
% load('saTCtCpSigma162NF.mat', 'Sigma')
% Sigma=Sigma(:,:,1:20,1:146);
%SigmaRule=double(Sigma<27.5).*double(Sigma>1); clear Sigma
% load('uvwSSHDailyDepth1rotated148F.mat', 'SSHa')
volS1=squeeze(nansum(nansum(inWAGSf.*SSHa(:,:,1:146).*repmat(RAC,[1 1 146])./5)));
vol2631=squeeze(nansum(nansum(inWAG263f.*SSHa(:,:,1:146).*repmat(RAC,[1 1 146])./5)));
vol2651=squeeze(nansum(nansum(inWAG265f.*SSHa(:,:,1:146).*repmat(RAC,[1 1 146])./5)));
vol26751=squeeze(nansum(nansum(inWAG2675f.*SSHa(:,:,1:146).*repmat(RAC,[1 1 146])./5)));
vol271=squeeze(nansum(nansum(inWAG27f.*SSHa(:,:,1:146).*repmat(RAC,[1 1 146])./5)));
vol2751=squeeze(nansum(nansum(inWAG275f.*SSHa(:,:,1:146).*repmat(RAC,[1 1 146])./5)));
for i=1:20
    i
    volS2(i,1:146)=squeeze(nansum(nansum(repmat(cellVol(:,:,i),[1 1 146]).*inWAGSf.*squeeze(SigmaRule(:,:,i,:))./5)));
    vol2632(i,1:146)=squeeze(nansum(nansum(repmat(cellVol(:,:,i),[1 1 146]).*inWAG263f.*squeeze(SigmaRule(:,:,i,:))./5)));
    vol2652(i,1:146)=squeeze(nansum(nansum(repmat(cellVol(:,:,i),[1 1 146]).*inWAG265f.*squeeze(SigmaRule(:,:,i,:))./5)));
    vol26752(i,1:146)=squeeze(nansum(nansum(repmat(cellVol(:,:,i),[1 1 146]).*inWAG2675f.*squeeze(SigmaRule(:,:,i,:))./5)));
    vol272(i,1:146)=squeeze(nansum(nansum(repmat(cellVol(:,:,i),[1 1 146]).*inWAG27f.*squeeze(SigmaRule(:,:,i,:))./5)));
    vol2752(i,1:146)=squeeze(nansum(nansum(repmat(cellVol(:,:,i),[1 1 146]).*inWAG275f.*squeeze(SigmaRule(:,:,i,:))./5)));
end
volS=volS1+sum(volS2,1).';
vol263=vol2631+sum(vol2632,1).';
vol265=vol2651+sum(vol2652,1).';
vol2675=vol26751+sum(vol26752,1).';
vol27=vol271+sum(vol272,1).';
vol275=vol2751+sum(vol2752,1).';
figure; plot(volS); hold all; plot(vol263); plot(vol265); plot(vol2675); plot(vol27); plot(vol275)
load('lagrangeWAGvolumes.mat', 'volWagF')
plot(volWagF,'m--','linewidth',2)
legend('surface','\sigma=26.3','\sigma=26.5','\sigma=26.75','\sigma=27','\sigma=27.5','3D')
volSetL=[volS vol263 vol265 vol2675 vol27 vol275];
dVdtL=diff(volSetL,1,1)./86400;

figure; errorbar(14:146,volWagPortion(14:146).',range(volSetL(14:146,:).')./2)
figure; errorbar(14:145,diff(volWagPortion(14:146))./86400,range(dVdtL(14:145,:).')./2)
errdVdt=range(dVdtL.')./2;
figure; errorbar(14:146,volWagPortion(14:146).',std(volSetL(14:146,:).'))
figure; errorbar(14:145,diff(volWagPortion(14:146))./86400,std(dVdtL(14:145,:).'))
errdVdt2=std(dVdtL.');
%% overall volume budget with errorbars figure
load('lagrangeWAGvolumeBudget.mat')
dVdtSet2=diff(volWagSet2,1,1)/86400;
dVdtMean2=mean(dVdtSet2,2);
dVdtErrb1=max(dVdtSet2,[],2)-dVdtMean;
dVdtErrb2=-min(dVdtSet2,[],2)+dVdtMean;

errTot1=dVdtErrb2+vertVolErr1(1:145)+errGateTot(1:145).';
errTot2=dVdtErrb1+vertVolErr2(1:145)+errGateTot(1:145).';

load('lagrangeWAGvolumesLayers.mat','errdVdt');
errTot1b=errdVdt2.'+abs(vertVolErr1(1:145))+abs(errGateTot(1:145)).';
errTot2b=errdVdt2.'+abs(vertVolErr2(1:145))+abs(errGateTot(1:145)).';

figure;
plot(14:146,gateFluxTot(14:146),'linewidth',2)
hold all
plot(vfluxsPF,'linewidth',2)
plot(14:146,vertVolMean(14:146),'linewidth',2)
plot(14:145,-dVdtMean2(14:145),'c','linewidth',2)
plot(14:145,gateFluxTot(14:145).'-dVdtMean2(14:145)+vfluxsPF(14:145).'+vertVolMean(14:145),'k--','linewidth',2)
legend('fill days','closed days','gate','precip','cross-iso','-dVdt','total')
%errorbar(14:146,gateFluxTot(14:146),errGateTot(14:146),'b.')
%errorbar(14:146,vertVolMean(14:146),vertVolErr2(14:146),vertVolErr1(14:146),'r.')
%errorbar(14:145,-dVdtMean2(14:145),dVdtErr1(14:145),dVdtErr2(14:145),'c.')
errorbar(14:145,gateFluxTot(14:145).'-dVdtMean2(14:145)+vfluxsPF(14:145).'+vertVolMean(14:145),errTot2b(14:145),errTot1b(14:145),'k.')
%errorbar(14:145,gateFluxTot(14:145).'-dVdtMean2(14:145)+vfluxsPF(14:145).'+vertVolMean(14:145),errTot2b(14:145),errTot1b(14:145),'m.')


figure;
% plot(14:146,gateFluxTot(14:146),'linewidth',2)
hold all
% plot(vfluxsPF,'linewidth',2)
% plot(14:146,vertVolMean(14:146),'linewidth',2)
plot(14:145,gateFluxTot(14:145).'+vfluxsPF(14:145).'+vertVolMean(14:145),'k--','linewidth',2)
plot(1:145,dVdtSet)
legend('gate+precip+crossIso','dVdt','dVdt','dVdt','dVdt','dVdt')%'gate','precip','cross-iso',
errorbar(14:145,gateFluxTot(14:145).'+vfluxsPF(14:145).'+vertVolMean(14:145),vertVolErr2(14:145)+errGateTot(14:145).',vertVolErr1(14:145)+errGateTot(14:145).','k.')
%% which days does the budget 'close'
totVflux(14:145)=gateFluxTot(14:145).'-dVdtMean2(14:145)+vfluxsPF(14:145).'+vertVolMean(14:145);
maxVflux(14:145)=totVflux(14:145)+abs(errTot1b(14:145)).';
minVflux(14:145)=totVflux(14:145)-abs(errTot2b(14:145)).';
figure; plot(totVflux); hold all; plot(maxVflux); plot(minVflux)
indexClosed=find(maxVflux>=0&minVflux<=0);
isClosed=maxVflux>=0&minVflux<=0;

%old, using errTot1,2
closedays=[14:15 17:19 22 27 29 32:33 35:37 41:43 47:48 50:52 57 59:60 64:66 71:72 75 79:81 83 85:88 95 98 102:103 106 108 111:112 116:117 119:120 129 132 135:136 138:141 144];
xclose=[13.5 16.5 21.5 26.5 28.5 31.5 34.5 40.5 46.5 49.5 56.5 58.5 63.5 70.5 74.5 78.5 82.5 84.5 94.5 97.5 101.5 105.5 107.5 110.5 115.5 118.5 128.5 131.5 134.5 137.5 143.5;...
        13.5 16.5 21.5 26.5 28.5 31.5 34.5 40.5 46.5 49.5 56.5 58.5 63.5 70.5 74.5 78.5 82.5 84.5 94.5 97.5 101.5 105.5 107.5 110.5 115.5 118.5 128.5 131.5 134.5 137.5 143.5;...
        15.5 19.5 22.5 27.5 29.5 33.5 37.5 43.5 48.5 52.5 57.5 60.5 66.5 72.5 75.5 81.5 83.5 88.5 95.5 98.5 103.5 106.5 108.5 112.5 117.5 120.5 129.5 132.5 136.5 141.5 144.5;
        15.5 19.5 22.5 27.5 29.5 33.5 37.5 43.5 48.5 52.5 57.5 60.5 66.5 72.5 75.5 81.5 83.5 88.5 95.5 98.5 103.5 106.5 108.5 112.5 117.5 120.5 129.5 132.5 136.5 141.5 144.5];
%yfill=repmat([-2e7;2e7; 2e7; -2e7],[1 13]);
yclose=repmat([-2e7;2e7;2e7;-2e7],[1 31]);
faceclose=[1:4;5:8;9:12;13:16;17:20;21:24;25:28;29:32;33:36;37:40;41:44;45:48;49:52;53:56;57:60;61:64;65:68;69:72;73:76;77:80;81:84;85:88;89:92;93:96;97:100;101:104;105:108;109:112;113:116;117:120;121:124];
figure; patch('faces',faceclose,'vertices',[xclose(:) yclose(:)],'facecolor','green','facealpha',0.3)
%% top-down vs bottom-up as source of errors, total lagrange volume budget

load('gateAdvection275indexed','gateFlux*')
gateFluxS263=gateFluxS263i;
%load('gateAdvectionS263','gateFlux*')
load('lagrangeWAGboundaryPortionFilledTB.mat','*dV*')
%index fix
dVdt(2:146)=dVdt(1:145); dVdtT(2:146)=dVdtT(1:145); dVdtB(2:146)=dVdtB(1:145);
errdVdt(2:146)=errdVdt(1:145);
load('vertVolIso275revIndexed','vertVol275*')
load('lagrangeWAGvolumeBudget.mat', 'vfluxsPF')
%index fix
vfluxsPF(2:147)=vfluxsPF(1:146);
clear gateFluxSign*

gateFlux1=0.5*gateFluxS263+0.5*gateFlux263S;
gateErr1=abs(gateFluxS263-gateFlux263S)/2;
gateFlux2=0.5*gateFlux263265+0.5*gateFlux265263;
gateErr2=abs(gateFlux263265-gateFlux265263)/2;
gateFlux3=0.5*gateFlux2652675+0.5*gateFlux2675265;
gateErr3=abs(gateFlux2652675-gateFlux2675265)/2;
gateFlux4=0.5*gateFlux267527+0.5*gateFlux272675;
gateErr4=abs(gateFlux267527-gateFlux272675)/2;
gateFlux5=0.5*gateFlux27275+0.5*gateFlux27527;
gateErr5=abs(gateFlux27275-gateFlux27527)/2;
gateFluxT=gateFluxS263+gateFlux263265+gateFlux2652675+gateFlux267527+gateFlux27275;
gateFluxB=gateFlux263S+gateFlux265263+gateFlux2675265+gateFlux272675+gateFlux27527;
gateErrTB=abs(gateFluxT-gateFluxB)/2;
gateErrSum=gateErr1+gateErr2+gateErr3+gateErr4+gateErr5;

vertVolSet=[vertVol275inWag27i  vertVol275inWag275i];%vertVol275inWag27v2
vertVolMean=mean(vertVolSet,2);
vertVolErr=range(vertVolSet,2)/2;
vertVolErrU=max(vertVolSet,[],2)-vertVolMean;
vertVolErrL=vertVolMean-min(vertVolSet,[],2);

errTot=vertVolErr(1:146).'+gateErrTB(1:146)+errdVdt(1:146).';
totFlux=vertVolMean(1:146).'+0.5*gateFluxT(1:146)+0.5*gateFluxB(1:146)+vfluxsPF(1:146)-dVdt.';
% load('lagrangeWAGvolumeBudget.mat', 'dVdtErr*')
% errBigU=vertVolErrU(1:145).'+gateErrSum(1:145)+dVdtErr2(1:145).';
% errBigL=vertVolErrL(1:145).'+gateErrSum(1:145)+dVdtErr1(1:145).';

%errBigUtb=vertVolErrU(1:145).'+gateErrSum(1:145)+errdVdt(1:145).';
%errBigLtb=vertVolErrL(1:145).'+gateErrSum(1:145)+errdVdt(1:145).';

figure;
plot(14:146,0.5*gateFluxT(14:146)+0.5*gateFluxB(14:146),'linewidth',2)
hold all
plot(vfluxsPF,'linewidth',2)
plot(14:146,vertVolMean(14:146),'linewidth',2)
plot(14:146,-dVdt(14:146),'c','linewidth',2)
plot(14:146,totFlux(14:146),'k--','linewidth',2)
legend('gate','precip','cross-iso','-dVdt','total')
errorbar(14:146,totFlux(14:146),errTot(14:146),'k.')
title('Volume budget, TB')

dVdtT(14)=0;
dVdtB(14)=0;
figure;
plot(14:146,gateFluxT(14:146),'linewidth',2)
hold all
plot(vfluxsPF,'linewidth',2)
plot(14:146,vertVol275inWag27i(14:146),'linewidth',2)
plot(14:146,dVdtT(14:146),'c','linewidth',2)
plot(14:146,gateFluxT(14:146)+vfluxsPF(14:146)+vertVol275inWag27i(14:146).','k--','linewidth',2)
legend('gate','precip','cross-iso','dVdt','gate+precip+crossiso')
errorbar(14:146,gateFluxT(14:146)+vfluxsPF(14:146)+vertVol275inWag27i(14:146).',vertVolErr(1:146).'+gateErrTB(1:146),'k.')
title('Volume budget, T')

figure;
plot(14:146,gateFluxB(14:146),'linewidth',2)
hold all
plot(vfluxsPF,'linewidth',2)
plot(14:146,vertVol275inWag275i(14:146),'linewidth',2)
plot(14:146,dVdtB(14:146),'c','linewidth',2)
plot(14:146,gateFluxB(14:146)+vfluxsPF(14:146)+vertVol275inWag275i(14:146).','k--','linewidth',2)
legend('gate','precip','cross-iso','-dVdt','total')
errorbar(14:146,gateFluxB(14:146)+vfluxsPF(14:146)+vertVol275inWag275i(14:146).',vertVolErr(14:146).'+gateErrTB(14:146),'k.')
title('Volume budget, B')


figure; plot(14:145,gateFluxT(14:145),'b','linewidth',2); hold all
plot(vfluxsPF,'g','linewidth',2)
plot(14:145,vertVol275inWag27(14:145),'r','linewidth',2)
plot(14:145,-dVdtT(14:145),'c','linewidth',2)
plot(14:145,-dVdtT(14:145).'+vertVol275inWag27(14:145).'+vfluxsPF(14:145)+gateFluxT(14:145),'k','linewidth',2)
plot(14:145,gateFluxB(14:145),'b--','linewidth',2); hold all
plot(14:145,vertVol275inWag275(14:145),'r--','linewidth',2)
plot(14:145,-dVdtB(14:145),'c--','linewidth',2)
plot(14:145,-dVdtB(14:145).'+vertVol275inWag275(14:145).'+vfluxsPF(14:145)+gateFluxB(14:145),'m--','linewidth',2)
legend('gate','precip','cross-iso','-dVdt','total')

maxVflux(14:146)=totFlux(14:146)+errTot(14:146);
minVflux(14:146)=totFlux(14:146)-errTot(14:146);
%indexClosed=find(maxVflux>=0&minVflux<=0);
isClosed=maxVflux>=0&minVflux<=0;
sum(isClosed)-13

%%
maxVflux(14:145)=totFlux(14:145)+errBigUtb(14:145);
minVflux(14:145)=totFlux(14:145)-errBigLtb(14:145);
indexClosed=find(maxVflux>=0&minVflux<=0);
isClosed=maxVflux>=0&minVflux<=0;
sum(isClosed)-13
%% how much is closed with range of dVdt from layers error option
%load('lagrangeWAGvolumeBudgetTB.mat')
load('lagrangeWAGvolumesLayers.mat','*dVdtL*','errdVdt');
filldays=[15 38:43 48 51:56 60:64 66:72 80:88 103 106:110 112 117:118 121 145:146];
errdVdtBig=max([errLdVdtL errUdVdtL errdVdt.'],[],2);
errBigU=vertVolErrU(1:145).'+gateErrSum(1:145)+errdVdtBig(1:145).';%errdVdt(1:145);
errBigL=vertVolErrL(1:145).'+gateErrSum(1:145)+errdVdtBig(1:145).';%errdVdt(1:145);
maxVflux(14:145)=totFlux(14:145)+errBigU(14:145);
minVflux(14:145)=totFlux(14:145)-errBigL(14:145);
indexClosed=find(maxVflux>=0&minVflux<=0);
isClosed=maxVflux>=0&minVflux<=0;
sum(isClosed)-13
sum(ismember(indexClosed,filldays))
sum(ismember(14:145,[indexClosed,filldays]))


filldays=[15 38:43 48 51:56 60:64 66:72 80:88 103 106:110 112 117:118 121 145:146];
xfill=[14.5 37.5 47.5 50.5 59.5 65.5 79.5 102.5 105.5 111.5 116.5 120.5 144.5; 14.5 37.5 47.5 50.5 59.5 65.5 79.5 102.5 105.5 111.5 116.5 120.5 144.5;...
       15.5 43.5 48.5 56.5 64.5 72.5 88.5 103.5 110.5 112.5 118.5 121.5 146.5; 15.5 43.5 48.5 56.5 64.5 72.5 88.5 103.5 110.5 112.5 118.5 121.5 146.5];
%index fix
filldays=filldays+1;
xfill=xfill+1;
   yfill=repmat([-4e7;4e7;4e7;-4e7],[1 13]);
facefill=[1:4;5:8;9:12;13:16;17:20;21:24;25:28;29:32;33:36;37:40;41:44;45:48;49:52];
xClose=[indexClosed(14:end)-0.5; indexClosed(14:end)-0.5; indexClosed(14:end)+0.5; indexClosed(14:end)+0.5];
yClose=repmat([-4e7;4e7;4e7;-4e7],[1 109]);
faceClose=reshape(1:4*109,4,[]);
figure; patch('faces',facefill,'vertices',[xfill(:) yfill(:)],'facecolor','red','facealpha',0.3)
hold all; patch('faces',faceClose,'vertices',[xClose(:) yClose(:)],'facecolor','green','facealpha',0.3)
%% 
volErr1=[abs(volWag(1:145)-volWagSnap(1:145)).'; abs(volWag(1:145)-volWagSnap2).'; abs(volWagSnap(1:145)-volWagSnap2).'];
volErr=max(volErr1,[],1);
volErr(146)=abs(volWag(146)-volWagSnap(146));
volDifErr(1:145)=sqrt(volErr(1:145).^2+volErr(2:146).^2);

figure; plot(15:145,fluxVolWag(15:145),'k','linewidth',2)
hold on; plot(15:145,fluxVolWag(15:145)+volDifErr(15:145).'./86400,'r--')
plot(15:145,fluxVolWag(15:145)-volDifErr(15:145).'./86400,'r--')

fvolErr1=[abs(fluxVolWag(1:144)-fluxVolWag2(1:144)).'; abs(fluxVolWag(1:144)-fluxVolWag22).'; abs(fluxVolWag2(1:144)-fluxVolWag22).'];
fvolErr=max(volErr1,[],1);
fvolErr(145)=abs(fluxVolWag(145)-fluxVolWag2(145));
figure; plot(15:145,fluxVolWag(15:145),'k','linewidth',2)
hold on; plot(15:145,fluxVolWag(15:145)+fvolErr(15:145).'./86400,'b--')
plot(15:145,fluxVolWag(15:145)-fvolErr(15:145).'./86400,'b--')


%%
filldays=[15 38:43 48 51:56 60:64 66:72 80:88 103 106:110 112 117:118 121 145:146];
xfill=[14.5 37.5 47.5 50.5 59.5 65.5 79.5 102.5 105.5 111.5 116.5 120.5 144.5; 14.5 37.5 47.5 50.5 59.5 65.5 79.5 102.5 105.5 111.5 116.5 120.5 144.5;...
       15.5 43.5 48.5 56.5 64.5 72.5 88.5 103.5 110.5 112.5 118.5 121.5 146.5; 15.5 43.5 48.5 56.5 64.5 72.5 88.5 103.5 110.5 112.5 118.5 121.5 146.5];
%yfill=repmat([-2e7;2e7; 2e7; -2e7],[1 13]);
yfill=repmat([0;4e12;4e12;0],[1 13]);
facefill=[1:4;5:8;9:12;13:16;17:20;21:24;25:28;29:32;33:36;37:40;41:44;45:48;49:52];
figure; patch('faces',facefill,'vertices',[xfill(:) yfill(:)],'facecolor','red','facealpha',0.3)
%figure; plot(gateFluxTot); hold all
%plot(vertVol275inWag27./2)
% %plot(vertVol275inWag27v2)
%plot(vfluxs)
%plot(-fluxVolWag);
%figure; 
hold all;
%plot(14:145,gateFluxTot(14:145).'+vertVol275inWag27(14:145))
plot(14:145,gateFluxTot(14:145).'+vertVol275inWag27(14:145)./2,'linewidth',2)
plot(14:145,fluxVolWag2(14:145),'linewidth',2)
plot(15:145,fluxVolWag(15:145)+volDifErr(15:145).'./86400,'b--')
plot(15:145,fluxVolWag(15:145)-volDifErr(15:145).'./86400,'b--')
%legend('gate','cross-\sigma=27.5','surface','-dVdt','total')
%hold all; plot(14:145,gateFluxTot(14:145).'+vertVol275inWag27v2(14:145));
%% find cross-isopycnal flux

%load('uvwNativeGridIsoDepth28.mat')
load('geometrySpinupSteady.mat','XC','YC','dInterface')
load('saTCtCpSigma162NF.mat', 'Sigma')
[~,~,nz,nt]=size(Sigma);
xmin=min(min(XC));
ymin=min(min(YC));
xcm=(XC-xmin).*111000.*cosd(YC);
ycm=(YC-ymin).*111000;
zcm=-0.5*dInterface(1:end-1)-0.5*dInterface(2:end);
dZ=diff(dInterface);
for k=1:150
    k
[gsy1,gsx1,gsz1]=gradient(Sigma(:,:,:,k),ycm,xcm,zcm(1:nz));
gsx(:,:,:,k)=gsx1;
gsy(:,:,:,k)=gsy1;
gsz(:,:,:,k)=gsz1;
end
clear Sigma
save('varyingSigmaGradientNF.mat','gs*','-v7.3')
% %isoDepthM=mean(isoDepth,3);
%%
 
 load('geometrySpinupSteady.mat','XC','YC','dInterface')
xmin=min(min(XC));
ymin=min(min(YC));
xcm=(XC-xmin).*111000.*cosd(YC);
ycm=(YC-ymin).*111000;
zcm=-0.5*dInterface(1:end-1)-0.5*dInterface(2:end);
dZ=diff(dInterface);
size(dZ)
load('varyingSigmaGradientNF.mat','gsx')
gsx=gsx(:,:,1:20,:);
%%
load('varyingSigmaGradientNF.mat','gsy')
gsy=gsy(:,:,1:20,:);
load('varyingSigmaGradientNF.mat','gsz')
gsz=gsz(:,:,1:20,:);
 load('uvwNativeGridIsoDepth28.mat','isoDepth')
for k=1:150
    k
for i=1:700
    for j=1:200
        gsxI(i,j,k)=interp1(zcm(1:20),squeeze(gsx(i,j,:,k)),-isoDepth(i,j,k));
        gsyI(i,j,k)=interp1(zcm(1:20),squeeze(gsy(i,j,:,k)),-isoDepth(i,j,k));
        gszI(i,j,k)=interp1(zcm(1:20),squeeze(gsz(i,j,:,k)),-isoDepth(i,j,k));
        index1=find(dInterface<=isoDepth(i,j,k),1,'last');
        dZi(i,j,k)=dZ(index1);
    end
end
end
clear gsx gsy gsz
save('sigma28sigmagradientNF.mat', '*I','dZi')
%I am projecting vec(u) onto grad(sigma)
%%
 %load('iso275depthNRsnap.mat')
 %isoDepthS=isoDepth;
% load('uvwNativeGridIsoDepth275NR.mat','*Iso','isoDepth')
 %load('lagrangeWAGboundary.mat')
% load('geometrySpinupSteady.mat','XC','YC','dInterface')
%load('sigma275sigmagradientNR.mat', '*I','dZi')
%load('geometrySpinupSteady','Angle*')
load('distancesAreas','RAC','dxg','dyg')
DXG=reshape(dxg,[700 200]);
DYG=reshape(dyg,[700 200]);
%I rotate uIso and vIso to eastward and northward to get correct
%product- orthogonal vectors in same directions
Urot=uIso(:,:,1:131).*repmat(AngleCS,[1 1 131]) - vIso(:,:,1:131).*repmat(AngleSN,[1 1 131]);  
Vrot=uIso(:,:,1:131).*repmat(AngleSN,[1 1 131]) + vIso(:,:,1:131).*repmat(AngleCS,[1 1 131]); 
clear uIso vIso
%dot product is magnitude of cross-isopycnal velocity
mag1=sqrt(gsxI.^2+gsyI.^2+gszI.^2);%magnitude of grad(sigma)
vecU=[reshape(Urot,1,[]);reshape(Vrot,1,[]);reshape(wIso(:,:,1:131),1,[])];
vecGI=[reshape(gsxI./mag1,1,[]);reshape(gsyI./mag1,1,[]);reshape(gszI./mag1,1,[])];
uDotGi=reshape(dot(vecU,vecGI),[700 200 131]);%magnitude of projected vec(u)

%east-north-up cross-iso vel components
udx=uDotGi.*gsxI./mag1;
udy=uDotGi.*gsyI./mag1;
udz=uDotGi.*gszI./mag1;
magCrossIsoVel=sqrt(udx.^2+udy.^2+udz.^2);
%native grid cross-iso vel components
udxN=udx.*repmat(AngleCS,[1 1 131])+udy.*repmat(AngleSN,[1 1 131]);
udyN=-udx.*repmat(AngleSN,[1 1 131])+udy.*repmat(AngleCS,[1 1 131]);
%uTest=Urot.*repmat(AngleCS,[1 1 131])+Vrot.*repmat(AngleSN,[1 1 131]);
%vTest=-Urot.*repmat(AngleSN,[1 1 131])+Vrot.*repmat(AngleCS,[1 1 131]);
% figure; quiver(XC(1:20:end,1:20:end),YC(1:20:end,1:20:end),uIso(1:20:end,1:20:end,20),vIso(1:20:end,1:20:end,20),'k');
% hold on;  quiver(XC(1:20:end,1:20:end),YC(1:20:end,1:20:end),uTest(1:20:end,1:20:end,20),vTest(1:20:end,1:20:end,20),'r');

%projecting onto grad(sigma) is outward direction, so negative is inward
vertVol275=squeeze(nansum(nansum((-udz.*repmat(RAC,[1 1 131])-udxN.*dZi.*repmat(DYG,[1 1 131])-udyN.*dZi.*repmat(DXG,[1 1 131])).*double(inWAG27(:,:,1:131))))); %
% vecU=[reshape(uIso(:,:,1:131),1,[]);reshape(vIso(:,:,1:131),1,[]);reshape(wIso(:,:,1:131),1,[])];
% vecGI=[reshape(gsxI,1,[]);reshape(gsyI,1,[]);reshape(gszI,1,[])];
% uDotGi=reshape(dot(vecU,vecGI),[700 200 131]);
 vertVol275v6=-squeeze(nansum(nansum((uDotGi.*repmat(RAC,[1 1 131])).*double(inWAG27(:,:,1:131)))));
 
 %% new 10/31/2017
 snapDiff=diff(isoDepthS(:,:,1:132),1,3);
 snapVel=snapDiff.*gszI./(mag1.*86400);
 surfaceZ=repmat(RAC,[1 1 131]);
 surfaceXN=dZi.*repmat(DYG,[1 1 131]);
 surfaceYN=dZi.*repmat(DXG,[1 1 131]);
 surfaceX=surfaceXN.*repmat(AngleCS,[1 1 131]) - surfaceYN.*repmat(AngleSN,[1 1 131]);  
 surfaceY=surfaceXN.*repmat(AngleSN,[1 1 131]) + surfaceYN.*repmat(AngleCS,[1 1 131]); 
 surfaces=reshape(dot([reshape(surfaceX,1,[]);reshape(surfaceY,1,[]);reshape(surfaceZ,1,[])],vecGI),[700 200 131]);
 vertVelCross=uDotGi-snapVel;
 vertVol=vertVelCross.*surfaces;
 vertVol275v7=squeeze(nansum(nansum(vertVol.*double(inWAG27(:,:,1:131)))));%mean -7e5
 %% distance isopycnal moves- not done yet
 %isoDepth must be from snapshots here
 deltaZ=diff(isoDepth,1,3);
 for k=2:149
 [dDdy1,dDdx1]=gradient(isoDepth(:,:,k),ycm,xcm);
 dDdy(:,:,k-1)=dDdy1;
 dDdx(:,:,k-1)=dDdx1;
 end
 deltaX=deltaZ(:,:,1:148)./dDdx(:,:,2:149);%this is eastward
 deltaY=deltaZ(:,:,1:148)./dDdy(:,:,2:149);
 distance=dot(vecGI,[reshape(deltaX,1,[]); reshape(deltaY,1,[]); reshape(deltaZ,1,[])]);
 uI=distance.*gsxI./(86400*mag1);%east
 vI=distance.*gsyI./(86400*mag1);%north
 wI=distance.*gszI./(86400*mag1);
 wCross=udz-wI;
 uXE=udx-uI;
 vXN=udy-vI;
uCross=uXE.*repmat(AngleCS,[1 1 131])+vXN.*repmat(AngleSN,[1 1 131]);
vCross=-uXE.*repmat(AngleSN,[1 1 131])+vXN.*repmat(AngleCS,[1 1 131]);
%% potential temperature and salinity contents
load('lagrangeWAGboundary.mat','inW*')
%load('varying148ts16levelsRho.mat','T');
load('fluxesT58day.mat','Ttend')

for i=2:58
    T1=T(:,:,:,i);
    Ttend1=Ttend(:,:,:,i);
    
    stayIn=(inWag3(:,:,:,i)==1)&(inWag3(:,:,:,i-1)==1);
    out1=(inWag3(:,:,:,i)==0)&(inWag3(:,:,:,i-1)==1);
    in1=(inWag3(:,:,:,i)==1)&(inWag3(:,:,:,i-1)==0);
    inNow=inWag3(:,:,:,i)==1;
    
    dTdt(i)=(sum(Ttend1(stayIn))+sum(T1(in1))-sum(T1(out1)))/86400;
    dTdt2(i)=squeeze(nansum(nansum(nansum(Ttend1(inNow)))))/86400;
end
clear Ttend T

load('varying148ts16levelsRho.mat','S');
load('fluxesS58day.mat','Stend')

for i=2:58
    S1=S(:,:,:,i);
    Stend1=Stend(:,:,:,i);
    
    stayIn=(inWag3(:,:,:,i)==1)&(inWag3(:,:,:,i-1)==1);
    out1=(inWag3(:,:,:,i)==0)&(inWag3(:,:,:,i-1)==1);
    in1=(inWag3(:,:,:,i)==1)&(inWag3(:,:,:,i-1)==0);
    inNow=inWag3(:,:,:,i)==1;
    
    dSdt(i)=(sum(Stend1(stayIn))+sum(S1(in1))-sum(S1(out1)))/86400;
    dSdt2(i)=squeeze(nansum(nansum(nansum(Stend1(inNow)))))/86400;
end
clear S Stend

save('lagrangeWAGcontentsTS58.mat','dTdt*','dSdt*')
%% T fluxes the Eulerian way, from diagnostics- load 
load('distancesAreas.mat','RAC','hFacC','dxc','dyc')
DXC=reshape(dxc,[700 200]); DYC=reshape(dyc,[700 200]);
load('geometrySpinupSteady.mat','dInterface')
nt=58;

dZ(1,1,1:46)=diff(dInterface);
cellVol=repmat(RAC,[1 1 46]).*hFacC.*repmat(dZ,[700 200 1]);
cellVol4=repmat(cellVol,[1 1 1 nt]);

load('fluxesT58day.mat')
load('lagrangeWAGboundary.mat','inWag3')
hold1=zeros([700 200 46 nt]);
[a1,a2,a3,a4]=size(inWag3);
hold1(1:a1,1:a2,1:a3,1:nt)=inWag3(:,:,:,1:nt);
inWag3=hold1; clear hold1
%% T fluxes the Eulerian way, from diagnostics
advTh=squeeze(nansum(nansum(nansum((AdvUt(1:end-1,1:end-1,:,:)-AdvUt(2:end,1:end-1,:,:)+AdvVt(1:end-1,1:end-1,:,:)...
    -AdvVt(1:end-1,2:end,:,:)).*inWag3(1:end-1,1:end-1,:,1:nt)./cellVol4(1:end-1,1:end-1,:,:)))));
advTz=squeeze(nansum(nansum(nansum(...
    (AdvWt(:,:,2:46,:)-AdvWt(:,:,1:45,:)).*inWag3(:,:,1:45,1:nt)./cellVol4(:,:,1:45,:)))));
difTz=squeeze(nansum(nansum(nansum((DifzEt(:,:,2:46,:)+DifIt(:,:,2:46,:)+KpT(:,:,2:46,:)-DifzEt(:,:,1:45,:)-DifIt(:,:,1:45,:)-KpT(:,:,1:45,:)).*inWag3(:,:,1:45,1:nt)./cellVol4(:,:,1:45,:)))));
difTh=squeeze(nansum(nansum(nansum((DifxEt(1:end-1,1:end-1,:,:)-DifxEt(2:end,1:end-1,:,:)+DifyEt(1:end-1,1:end-1,:,:)...
    -DifyEt(1:end-1,2:end,:,:)).*inWag3(1:end-1,1:end-1,:,1:nt)./cellVol4(1:end-1,1:end-1,:,:)))));
for i=1:nt
surfTtend3(i)=nansum(nansum(squeeze(inWag3(:,:,1,i)).*(-squeeze(WT(:,:,1,i)))./(dZ(1,1,1).*hFacC(:,:,1))));
end

load('varying148ts16levelsRho.mat', 'Rho')
Rho=squeeze(Rho(:,:,1,1:nt));
load('cp16levels148.mat', 'cp')
cp=squeeze(cp(:,:,1,1:nt));
surfQ=squeeze(nansum(nansum(Tflux.*inWag3(:,:,1,1:nt)./(dZ(1,1,1).*repmat(hFacC(:,:,1),[1 1 nt]).*Rho.*cp))));
clear Rho cp Tflux

load('varying148ts16levelsRho.mat', 'T');
T=T(:,:,1:15,1:nt);
load('dwadz148levels15.mat','dWdZ')
dWdZ=dWdZ(:,:,:,1:nt);
tdwdz=squeeze(nansum(nansum(nansum(T.*dWdZ.*inWag3(:,:,1:15,1:nt)))));
clear T dWdZ cellVol DXC DYC hFacC RAC

load('lagrangeWAGcontentsTS58.mat','dTdt*')
save('lagrangeWagTempBudget58.mat')
%% S fluxes, Eulerian way, load
clear
load('distancesAreas.mat','RAC','hFacC','dxc','dyc')
DXC=reshape(dxc,[700 200]); DYC=reshape(dyc,[700 200]);
load('geometrySpinupSteady.mat','dInterface')
nt=58;

dZ(1,1,1:46)=diff(dInterface);
cellVol=repmat(RAC,[1 1 46]).*hFacC.*repmat(dZ,[700 200 1]);
cellVol4=repmat(cellVol,[1 1 1 nt]);

load('fluxesS58day.mat')
load('lagrangeWAGboundary.mat','inWag3')
hold1=zeros([700 200 46 nt]);
[a1,a2,a3,a4]=size(inWag3);
hold1(1:a1,1:a2,1:a3,1:nt)=inWag3(:,:,:,1:nt);
inWag3=hold1; clear hold1

%% S fluxes from diagnostics
AdvSh=squeeze(nansum(nansum(nansum(...
    (AdvUs(1:end-1,1:end-1,:,:)-AdvUs(2:end,1:end-1,:,:)+AdvVs(1:end-1,1:end-1,:,:)...
    -AdvVs(1:end-1,2:end,:,:)).*inWag3(1:end-1,1:end-1,:,1:nt)./cellVol4(1:end-1,1:end-1,:,:)...
        ))));
AdvSz=squeeze(nansum(nansum(nansum(...
    (AdvWs(:,:,2:46,:)-AdvWs(:,:,1:45,:)).*inWag3(:,:,1:45,1:nt)./cellVol4(:,:,1:45,:)...
    ))));
difSz=squeeze(nansum(nansum(nansum((DifEs(:,:,2:46,:)+DifIs(:,:,2:46,:)+KpS(:,:,2:46,:)-DifEs(:,:,1:45,:)-DifIs(:,:,1:45,:)-KpS(:,:,1:45,:)).*inWag3(:,:,1:45,1:nt)./cellVol4(:,:,1:45,:)))));
difSh=squeeze(nansum(nansum(nansum((DifxEs(1:end-1,1:end-1,:,:)-DifxEs(2:end,1:end-1,:,:)+DifyEs(1:end-1,1:end-1,:,:)...
    -DifyEs(1:end-1,2:end,:,:)).*inWag3(1:end-1,1:end-1,:,1:nt)./cellVol4(1:end-1,1:end-1,:,:)))));
for i=1:nt
surfStend3(i)=nansum(nansum(inWag3(:,:,1,i).*(-squeeze(WS(:,:,1,i)))./(dZ(1,1,1).*hFacC(:,:,1))));
end
load('varying148ts16levelsRho.mat', 'Rho')
Rho=squeeze(Rho(:,:,1,1:nt));
surfF=squeeze(nansum(nansum(Sflux.*squeeze(inWag3(:,:,1,1:nt))./(dZ(1,1,1).*repmat(hFacC(:,:,1),[1 1 nt]).*Rho))));

load('varying148ts16levelsRho.mat', 'S');
S=S(:,:,1:15,1:nt);
load('dwadz148levels15.mat','dWdZ')
dWdZ=dWdZ(:,:,:,1:nt);
sdwdz=squeeze(nansum(nansum(nansum(S.*dWdZ.*inWag3(:,:,1:15,1:nt)))));
clear S dWdZ Sflux cellVol DXC DYC hFacC RAC

load('lagrangeWAGcontentsTS58.mat','dSdt*')
save('lagrangeWagSalBudget58.mat')
%% plots budget
load('lagrangeWAGcontentsTS58.mat','dSdt*')
load('lagrangeWagSalBudget58.mat','AdvS*','sdwdz','difS*','surfStend*','surfF')
figure; plot(-dSdt); hold all
plot(AdvSh+sdwdz); plot(AdvSz-sdwdz);
plot(difSh); plot(difSz); plot(surfStend3'+surfF)
plot(AdvSh(1:45)+AdvSz(1:45)+difSh(1:45)+difSz(1:45)+surfStend3(1:45)'+surfF(1:45))
legend('-dS/dt','AdvH','AdvV','difH','difV','surf','tot fluxes')

figure; plot(-dSdt2); hold all
plot(AdvSh+sdwdz); plot(AdvSz-sdwdz);
plot(difSh); plot(difSz); plot(surfStend3'+surfF)
plot(AdvSh(1:45)+AdvSz(1:45)+difSh(1:45)+difSz(1:45)+surfStend3(1:45)'+surfF(1:45))
legend('-dS/dt','AdvH','AdvV','difH','difV','surf','tot fluxes')

%% gate advective fluxes from lagrange calculations before, compared to the advective ones from diagnostics
load('gateAdvection26527.mat')
load('gateAdvection26265.mat')
load('gateAdvection27275.mat')
load('gateAdvectionSurface.mat')





