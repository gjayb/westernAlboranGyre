%load('isopycnalDepths2.mat')
load('uvwDailyNativeNF.mat')
u=U; clear U;
v=V; clear V;
w=W; clear W
load('iso29depthNF.mat')
addpath('/nobackup1/gbrett/mStuff/')
disp('load done')
%%
su=size(u)
sv=size(v)
sw=size(w)
%%
load('geometrySpinupSteady.mat','dInterface','d')
dBin=0.5.*dInterface(1:end-1)+0.5.*dInterface(2:end);
%depths1=[2.5 3:dBin(16)];
zw=dInterface(1:end-1);
    isopyc=29; %potential density
%isoDepth=iso29depth;
%clear iso275depth iso28depth iso26depth iso27depth iso265depth iso285depth iso29depth iso295depth
    nt=162;
    uIso=zeros([700 200 nt]); vIso=uIso; wIso=uIso;
%%    
    for i=1:nt
    i
        for j=1:700
            for k=1:200
                if isoDepth(j,k,i)>0 && isoDepth(j,k,i)<d(j,k) %because isodepth is based on top 16 layers only
                uIso(j,k,i)=interp1(dBin,squeeze(u(j,k,:,i)),isoDepth(j,k,i));
                vIso(j,k,i)=interp1(dBin,squeeze(v(j,k,:,i)),isoDepth(j,k,i));
                wIso(j,k,i)=interp1(zw,squeeze(w(j,k,:,i)),isoDepth(j,k,i));
                end
            end
        end
    
    end
    
fn='uvwNativeGridIsoDepth29.mat';
save(fn,'iso*','*Iso','-v7.3')
