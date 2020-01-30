%% load
load('lobesByHandSigma265.mat')

load('lobes14dayIso265days1to50small.mat')

clear Sigma CT Angle* density* open* wagR np* notaBene nt2 nt3 perimeters areas
clear delZ dist* clusterorder ans *U *V *G C DRC ii j jc jcN jj kk k

%% plot, for a given day, the lobes on other layers and the manifolds on sigma 265

for i=19:25
    
    xlobeHS(xlobeHS==0)=nan;
ylobeHS(ylobeHS==0)=nan;
figure; plot(squeeze(xlobeHS(i-8,:,:)).',squeeze(ylobeHS(i-8,:,:)).','m','LineWidth',2)
xlobeH26(xlobeH26==0)=nan;
ylobeH26(ylobeH26==0)=nan;
hold on; plot(squeeze(xlobeH26(i,:,:)).',squeeze(ylobeH26(i,:,:)).','c')

plot(lontrFall(i,:),lattrFall(i,:),'r')
    hold on; plot(lontrBall(i,:),lattrBall(i,:),'b')
    plot(lonCoast,latCoast,'k')
    axis([-6.5 -2 35 37])
    title(num2str(i))
    
   %legend('surface','\sigma=26','\sigma=26.5 F','\sigma=26.5 B')
    
    
    
end
