
%clear,clc
%close all
addpath('../mStuff')
load('geometrySpinupSteady.mat')
DRF = rdmds('DRF');
hFacW = rdmds('hFacW');
DYG = rdmds('DYG');

% I-index of the section whereby volume transport will be computed 
% (e.g. 88)  
I_index = 186;%250; 

Ny = 200; % Number of grid points in "y" direction
Nz = 46; % Number of vertical levels

% Matrix including the area of individual cells
Ar = zeros(Ny,Nz);

for j = 1:Ny
    for k = 1:Nz
        Ar(j,k) = DRF(k)*hFacW(I_index,j,k)*DYG(I_index,j);
    end
end

inflow=zeros([426 1]);
outflow=inflow;
interface1=inflow;
interface2=inflow;
umean=zeros([46 426]);
timeday=8640:8640:3680640;
for i=1:426
    i
    iter=timeday(i); % Iteration number of the model output
    % Reading u-velocity
    U = rdmds('Uave',iter);

    Usec = squeeze(U(I_index,:,:,:));
    Usec([1:10,194:200],:)=0;%correction for grid folding on itself
    Tr = Ar.*Usec;%repmat(Ar,[1 1 3]).*Usec;

    aux = Tr>0;
    inflow(i) = sum(sum(Tr.*aux))/1e6; % inflow in Sv
    aux = Tr<0;
    outflow(i) = sum(sum(Tr.*aux))/1e6; % outflow in Sv
    umean(:,i)=mean(Usec,1);
    interface1(i)=dInterface(find(umean(:,i)>0,1,'last')+1);
    interface2(i)=dInterface(find(umean(:,i)<0,1,'first'));
end

%figure; plot(abs(outflow)); title('Outflow in Sv'); xlabel('Days, 1=Nov 1, 2007')
%figure; plot(interface1); %hold all; plot(interface2); legend('interface1','interface2'); 
%xlabel('Days, 1=Nov 1, 2007'); title('Depth of Interface, below last inflow')
save('transportGibraltarDaily.mat')
 






