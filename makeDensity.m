addpath('../mStuff')
addpath('../mStuff/gsw')

load('salinityDaily.mat'); S=S(:,:,:,1:148);
load('potentialTempDaily.mat'); T=T(:,:,:,1:148);
disp('load done')
size(S)
size(T)

for i=1:148
i
CT(:,:,:,i)=gsw_CT_from_pt(S(:,:,:,i),T(:,:,:,i));
end
disp('CT done')
size(CT)
clear T
for i=1:148
Sigma(:,:,:,i)=gsw_sigma0(S(:,:,:,i),CT(:,:,:,i));
end
disp('Sigma done')
size(Sigma)
clear S
save('sigmaCT148.mat','-v7.3')
disp('done')

