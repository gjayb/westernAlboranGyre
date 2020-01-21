
addpath('/nobackup1/gbrett/mStuff/')
addpath('/nobackup1/gbrett/mStuff/library/') 
times=8640:8640:(162*8640)
nt=length(times)
days=1:162;
% PotTave=rdmds('Tave',times);
% PractSave=rdmds('Save',times);
% RhoAave=rdmds('RhoAAve',times);
% Rho0=rdmds('RhoRef');
% Phi=rdmds('Pave',times);

g=9.81;
load('geometrySpinupSteady.mat','*C','dInterface')
dBin=0.5.*dInterface(1:end-1)+0.5.*dInterface(2:end);

% load('rhoSSH148day.mat','Rho0')
% Pref=g*Rho0.*reshape(dBin,[1 1 46]);
% PressDbar=1e-4*(Phi.*repmat(Rho0,[700 200 1 nt])+repmat(Pref,[700 200 1 nt]));
% disp('save 1')
% save('tsPress162NF.mat','-v7.3')
% clear Pref RhoRef Phi RhoAave Rho0
% load('tsPress162NF.mat','PractSave','PressDbar','PotTave')
 
% for k=1:nt
%     k
%     for di=1:46
%         SA(:,:,di,k)=gsw_SA_from_SP(PractSave(:,:,di,k),PressDbar(:,:,di,k),XC,YC);
 %        CT(:,:,di,k)=gsw_CT_from_pt(SA(:,:,di,k),PotTave(:,:,di,k));
 %        tInsitu(:,:,di,k)=gsw_t_from_CT(SA(:,:,di,k),CT(:,:,di,k),PressDbar(:,:,di,k));
%         CP(:,:,di,k)=gsw_cp_t_exact(SA(:,:,di,k),tInsitu(:,:,di,k),PressDbar(:,:,di,k));
%         Sigma(:,:,di,k)=gsw_sigma0(SA(:,:,di,k),CT(:,:,di,k));
%         
%     end
% end
% clear PractSave PotTave PressDbar
% disp('save 2')
% save('saTCtCpSigma162NR.mat','-v7.3')
 
% clear SA CT tInsitu CP Sigma
% PotT=rdmds('T',times);
% PractS=rdmds('S',times);
% load('tsPress162NF.mat','PressDbar')
% for k=1:nt
%     k
%     for di=1:46
%         SA(:,:,di,k)=gsw_SA_from_SP(PractS(:,:,di,k),PressDbar(:,:,di,k),XC,YC);
%         CT(:,:,di,k)=gsw_CT_from_pt(SA(:,:,di,k),PotT(:,:,di,k));
%         Sigma(:,:,di,k)=gsw_sigma0(SA(:,:,di,k),CT(:,:,di,k));
%         
%     end
% end
 
% disp('save 3')
% save('tsSigmaSnapshots162NF.mat','-v7.3')
%% 
momHD=rdmds('momHDave',times);
%load('tsPress162NF.mat','PotTave');
%ThDiv=PotTave.*momHD;
%size(ThDiv)
%clear PotTave
%load('tsPress162NF.mat','PractSave');
%ShDiv=PractSave.*momHD;
%size(ShDiv)
%clear PractSave
%disp('save 4')
%save('divTS162NR.mat','ShDiv','ThDiv','-v7.3')
SSH=rdmds('SSHave',times)
save('sshMomHD.mat','-v7.3')
