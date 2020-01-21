%load('saTCtCpSigma162NR.mat', 'Sigma')
load('tsSigmaSnapshots162NF.mat','Sigma')
load('geometrySpinupSteady.mat')
%load('sigmaCT148.mat')
%size(CT)

size(Sigma)
%%
dBin=0.5.*dInterface(1:end-1)+0.5.*dInterface(2:end);
isopycs=[26.3 26.5 26.75 27 27.5 28 28.5 28.9 29 29.1];
isoStr=[263 265 2675 27 275 28 285 289 29 291];

%depths1=[2.5 3:dBin(end)];
for iiso=9:10
    isopyc=isopycs(iiso); %potential density
    nt=162;
    isoDepth=zeros([700 200 nt]);
%%    
    for i=1:nt
    i
        for j=1:700
            for k=1:200
                if d(j,k)>0
               depths1=[2.5 3:dBin(nLayers(j,k))];
               rho11=interp1(dBin,squeeze(Sigma(j,k,:,i)),depths1);

                di12=find(abs(rho11-isopyc)==min(abs(rho11-isopyc)),1,'first');       
                %di21=find(rho22<isopyc,1,'last');
                %di22=find(abs(rho22-isopyc)==min(abs(rho22-isopyc)),1,'first');
                if isempty(di12)
                    isoDepth(j,k,i)=0;
                else
                    isoDepth(j,k,i)=depths1(di12);
                end
                else
                    isoDepth(j,k,i)=0;
                end
            end
        end
    
    end
    
%figure; pcolor(XC,YC,mean(isoDepth,3)); shading 'flat'
%figure; pcolor(XC,YC,std(isoDepth,0,3)); shading 'flat'
fn=strcat('iso',num2str(isoStr(iiso)),'depthNFsnap.mat')
save(fn,'isoDepth')
end
