load('isopycnalDepths2.mat')
%load('uvwDailyNativeGrid.mat')
addpath('/nobackup1/gbrett/mStuff/')
disp('load 1 done')
%%
% su=size(u)
% sv=size(v)
% sw=size(w)
 %%
dBin=0.5.*dInterface(1:end-1)+0.5.*dInterface(2:end);
%depths1=[2.5 3:dBin(16)];
zw=dInterface(1:end-1);
    isopyc=29; %potential density
isoDepth=iso29depth;
clear iso275depth iso28depth iso26depth iso27depth iso265depth iso285depth iso29depth iso295depth
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

%% make S, T, cp on isopycnals: loads
load('varying148ts16levelsRho.mat','S')
load('varying148ts16levelsRho.mat','T')
 %at the surface, in-situ t is potential temp, which is this T
load('varying148ts16levelsRho.mat','Rho')
disp('load 2 done')
%% make S, T, cp on isopycnals
%% make  cp 
isoDepth=iso26depth;
 nt=148;
 sIso=zeros([700 200 nt]); tIso=sIso; rhoIso=sIso; cpIso=sIso;
 
%  cp=zeros([700 200 16 nt]);
%  for i=1:16
% cp(:,:,i,:)=gsw_cp_t_exact(S(:,:,i,:),T(:,:,i,:),dInterface(i));
%  end
 
 disp('cp done')
%% make S, T, cp on isopycnals
dBin=0.5*(dInterface(1:end-1)+dInterface(2:end));
 for i=1:nt
    i
        for j=1:700
            for k=1:200
                if isoDepth(j,k,i)>0 && isoDepth(j,k,i)<d(j,k) %because isodepth is based on top 16 layers only
                sIso(j,k,i)=interp1(dBin(1:16),squeeze(S(j,k,:,i)),isoDepth(j,k,i));
                tIso(j,k,i)=interp1(dBin(1:16),squeeze(T(j,k,:,i)),isoDepth(j,k,i));
                rhoIso(j,k,i)=interp1(dBin(1:16),squeeze(Rho(j,k,:,i)),isoDepth(j,k,i));
                cpIso(j,k,i)=interp1(dBin(1:16),squeeze(cp(j,k,:,i)),isoDepth(j,k,i));
                end
            end
        end
    
 end

fn='tsrhocpNativeGridIsoDepth26.mat';
save(fn,'*Iso','-v7.3')

