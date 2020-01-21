%calculates ftles and arclengths

load('trajMITgcmNovDecA2week.mat');
tIntegrate=range(tmesh);

for intLoop=1:31
    intLoop
    iL=(intLoop-2)/2;
    %trajectories
    lattr=squeeze(lattrF(:,:,iL));
    lontr=squeeze(lontrF(:,:,iL));
    lattr2=squeeze(lattrB(:,:,iL));
    lontr2=squeeze(lontrB(:,:,iL));
    %convert to meters
    xtr=(lontr-xmin*ones(size(lontr))).*111000.*cosd(lattr); 
    ytr=(lattr-ymin*ones(size(lattr))).*111000;
    xtr2=(lontr2-xmin*ones(size(lontr2))).*111000.*cosd(lattr2); 
    ytr2=(lattr2-ymin*ones(size(lattr2))).*111000;
    
    
    xend=reshape(xtr(end,:),NY,NX); yend=reshape(ytr(end,:),NY,NX); 
     %xend=reshape(xtr(150,:),NY,NX); yend=reshape(ytr(150,:),NY,NX); 

    dx11=xend(2:end-1,3:end)-xend(2:end-1,1:end-2); 
    dx12=xend(3:end,2:end-1)-xend(1:end-2,2:end-1);
    dy21=yend(2:end-1,3:end)-yend(2:end-1,1:end-2);
    dy22=yend(3:end,2:end-1)-yend(1:end-2,2:end-1);
    dx0 =xinM(2:end-1,3:end)-xinM(2:end-1,1:end-2);
    dy0 =yinM(3:end,2:end-1)-yinM(1:end-2,2:end-1);

    a11=(dx11./dx0).^2.0+(dy21./dx0).*(dy21./dx0);	
    a12=(dx11./dx0).*(dx12./dy0)+(dy21./dx0).*(dy22./dy0);
    a21=(dx12./dy0).*(dx11./dx0)+(dy22./dy0).*(dy21./dx0);
    a22=(dx12./dy0).*(dx12./dy0)+(dy22./dy0).^2.0;
    lambda1=(a22+a11)/2.0+sqrt(((a22+a11)/2.0).^2.0-(a11.*a22-a12.*a21));
    lambda2=(a22+a11)/2.0-sqrt(((a22+a11)/2.0).^2.0-(a11.*a22-a12.*a21));
    ftleP(:,:,iL)=log(sqrt(max(lambda1,lambda2)))./abs(tIntegrate);
    
    xend2=reshape(xtr2(end,:),NY,NX); yend2=reshape(ytr2(end,:),NY,NX); 

    dx11=xend2(2:end-1,3:end)-xend2(2:end-1,1:end-2); 
    dx12=xend2(3:end,2:end-1)-xend2(1:end-2,2:end-1);
    dy21=yend2(2:end-1,3:end)-yend2(2:end-1,1:end-2);
    dy22=yend2(3:end,2:end-1)-yend2(1:end-2,2:end-1);
    dx0 =xinM(2:end-1,3:end)-xinM(2:end-1,1:end-2);
    dy0 =yinM(3:end,2:end-1)-yinM(1:end-2,2:end-1);

    a11=(dx11./dx0).^2.0+(dy21./dx0).*(dy21./dx0);	
    a12=(dx11./dx0).*(dx12./dy0)+(dy21./dx0).*(dy22./dy0);
    a21=(dx12./dy0).*(dx11./dx0)+(dy22./dy0).*(dy21./dx0);
    a22=(dx12./dy0).*(dx12./dy0)+(dy22./dy0).^2.0;
    lambda1=(a22+a11)/2.0+sqrt(((a22+a11)/2.0).^2.0-(a11.*a22-a12.*a21));
    lambda2=(a22+a11)/2.0-sqrt(((a22+a11)/2.0).^2.0-(a11.*a22-a12.*a21));
    ftleN(:,:,iL)=log(sqrt(max(lambda1,lambda2)))./abs(tIntegrate);
    
    dxdt=diff(xtr);
    dydt=diff(ytr);
    ds=sqrt(dxdt.^2 +dydt.^2);
    lengths(:,:,iL)=nansum(ds);
    
    dxdt2=diff(xtr2);
    dydt2=diff(ytr2);
    ds2=sqrt(dxdt2.^2 +dydt2.^2);
    lengths2(:,:,iL)=nansum(ds2);
end    
    fn='ftleADepth1ND2week.mat';
    save(fn,'-v7.3')

