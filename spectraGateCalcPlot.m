%spectra of flux(es), point velocity, forcing
Fs=1/86400; %sampling frequency
Ts = 1/Fs;             % Sampling period       
LgF = 138;             % Length of signal, gate flux
LgV = 148;          %length of signal, velocities
tgF = (0:LgF-1)*Ts;        % Time vector
tgV = (0:LgV-1)*Ts;        % Time vector
fgF = Fs*(0:(LgF/2))/LgF; %frequency vector
fgV = Fs*(0:(LgV/2))/LgV; %frequency vector
%% flux
load('gateAdvection26527.mat')
load('gateAdvection26265.mat')
load('gateAdvection27275.mat')
load('gateAdvectionSurface.mat')

figure; %volume
plot(gateFluxS26./1e6); hold all;
plot(gateFlux26265./1e6); plot(gateFlux26527./1e6); plot(gateFlux27275./1e6)
plot((gateFlux26265+gateFlux26527+gateFlux27275+gateFluxS26(1:138))./1e6,'k','LineWidth',2)
title('Gate Volume Flux'); xlabel('simulation day'); ylabel('volume flux, Sverdrups')
legend('Surface to \sigma=26','\sigma=26 to \sigma=26.5','\sigma=26.5 to \sigma=27','\sigma=27 to \sigma=27.5','sum')

gateFluxT=gateFlux26265+gateFlux26527+gateFlux27275+gateFluxS26(1:138);
[PgS26,FgS26]=pwelch(gateFluxS26-nanmean(gateFluxS26),[],[],[],1/86400);
[Pg26265,Fg26265]=pwelch(gateFlux26265,[],[],[],1/86400);
[Pg26527,Fg26527]=pwelch(gateFlux26527,[],[],[],1/86400);
[Pg27275,Fg27275]=pwelch(gateFlux27275,[],[],[],1/86400);
[PgT,FgT]=pwelch(gateFluxT,[],[],[],1/86400);

fftGS26=fft(gateFluxS26); 
holdvar = abs(fftGS26/LgF);
fftGS26 = holdvar(1:LgF/2+1);
fftGS26(2:end-1) = 2*fftGS26(2:end-1);

fftG26265=fft(gateFlux26265); 
holdvar = abs(fftG26265/LgF);
fftG26265 = holdvar(1:LgF/2+1);
fftG26265(2:end-1) = 2*fftG26265(2:end-1);

fftG26527=fft(gateFlux26527); 
holdvar = abs(fftG26527/LgF);
fftG26527 = holdvar(1:LgF/2+1);
fftG26527(2:end-1) = 2*fftG26527(2:end-1);

fftG27275=fft(gateFlux27275); 
holdvar = abs(fftG27275/LgF);
fftG27275 = holdvar(1:LgF/2+1);
fftG27275(2:end-1) = 2*fftG27275(2:end-1);

fftGT=fft(gateFluxT); 
holdvar = abs(fftGT/LgF);
fftGT = holdvar(1:LgF/2+1);
fftGT(2:end-1) = 2*fftGT(2:end-1);

%% flux plot
figure
loglog(FgS26*86400,PgS26)
hold all
loglog(Fg26265*86400,Pg26265)
loglog(Fg26527*86400,Pg26527)
loglog(Fg27275*86400,Pg27275)
loglog(FgT*86400,PgT)
title('Gate Volume Flux Spectra'); xlabel('Frequency in 1/day'); ylabel('Power')
legend('Surface to \sigma=26','\sigma=26 to \sigma=26.5','\sigma=26.5 to \sigma=27','\sigma=27 to \sigma=27.5','sum')


figure
loglog(fgF*86400,fftGS26)
hold all
loglog(fgF*86400,fftG26265)
loglog(fgF*86400,fftG26527)
loglog(fgF*86400,fftG27275)
loglog(fgF*86400,fftGT)
title('Gate Volume Flux Spectra'); xlabel('Frequency in 1/day'); ylabel('Power')
legend('Surface to \sigma=26','\sigma=26 to \sigma=26.5','\sigma=26.5 to \sigma=27','\sigma=27 to \sigma=27.5','sum')
%% flux- fft sections
load('gateAdvection26527.mat')
load('gateAdvection26265.mat')
load('gateAdvection27275.mat')
load('gateAdvectionSurface.mat')

Fs=1/86400; %sampling frequency
Ts = 1/Fs;             % Sampling period       
LgS = 28;             % Length of signal, gate flux subsets
tgS = (0:LgS-1)*Ts;        % Time vector
fgS = Fs*(0:(LgS/2))/LgS; %frequency vector

fluxes=[gateFluxS26(1:140);0 gateFlux26265 0;0 gateFlux26527 0; 0 gateFlux27275 0;0 gateFluxT 0];
%fluxes=[gateFluxS26(1:130); gateFlux26265(1:130); gateFlux26527(1:130); gateFlux27275(1:130); gateFluxT(1:130)];

clear fftFluxes
for i=1:5 %depths
    for j=1:5 %time
        fftF=fft(fluxes(i,1+28*(j-1):28*j)); %1+26*(j-1):26*j%
        holdvar = abs(fftF/LgS);
        fftFluxes(i,j,:) = holdvar(1:LgS/2+1);
        fftFluxes(i,j,2:end-1) = 2*fftFluxes(i,j,2:end-1);
    end %for j, time
    
end %for i, depth

spectraMeanF=squeeze(mean(fftFluxes,2));
spectraMaxF=squeeze(max(fftFluxes,[],2));
spectraMinF=squeeze(min(fftFluxes,[],2));
%% plot flux spectra from fft sections with error bars
figure; 
errorbar(fgS*86400,spectraMeanF(1,:),spectraMeanF(1,:)-spectraMinF(1,:),spectraMaxF(1,:)-spectraMeanF(1,:))
hold all
for ii=2:5
    errorbar(fgS*86400,spectraMeanF(ii,:),spectraMeanF(ii,:)-spectraMinF(ii,:),spectraMaxF(ii,:)-spectraMeanF(ii,:))
end %ii
title('fft spectra with errors')
set(gca,'YScale','log')
set(gca,'XScale','log')
xlabel('frequency in 1/day')

%% flux- pwelch sections
load('gateAdvection26527.mat')
load('gateAdvection26265.mat')
load('gateAdvection27275.mat')
load('gateAdvectionSurface.mat')

Fs=1/86400; %sampling frequency
Ts = 1/Fs;             % Sampling period       
LgS = 28;             % Length of signal, gate flux subsets
tgS = (0:LgS-1)*Ts;        % Time vector
fgS = Fs*(0:(LgS/2))/LgS; %frequency vector

fluxes=[gateFluxS26(1:140);0 gateFlux26265 0;0 gateFlux26527 0; 0 gateFlux27275 0;0 gateFluxT 0];

for i=1:5 %depths
    for j=1:5 %time
        pwFluxes(i,j,:) = pwelch(fluxes(i,1+28*(j-1):28*j),[],[],fgS,Fs);
    end %for j, time
    
end %for i, depth

spectraMeanFp=squeeze(mean(pwFluxes,2));
spectraMaxFp=squeeze(max(pwFluxes,[],2));
spectraMinFp=squeeze(min(pwFluxes,[],2));
spectraStdFp=squeeze(std(pwFluxes,0,2));
%% plot flux spectra from pwelch sections with error bars
figure; 
errorbar(fgS*86400,spectraMeanFp(1,:),spectraStdFp(1,:))%spectraMeanFp(1,:)-spectraMinFp(1,:),spectraMaxFp(1,:)-spectraMeanFp(1,:))
hold all
for ii=2:5
    errorbar(fgS*86400,spectraMeanFp(ii,:),spectraStdFp(ii,:))%spectraMeanFp(ii,:)-spectraMinFp(ii,:),spectraMaxFp(ii,:)-spectraMeanFp(ii,:))
end %ii
title('Pwelch spectra with errors')

set(gca,'YScale','log')
set(gca,'XScale','log')
xlabel('frequency in 1/day')

%% flux in 2 sections
load('gateAdvection26527.mat')
load('gateAdvection26265.mat')
load('gateAdvection27275.mat')
load('gateAdvectionSurface.mat')

Fs=0.5/86400; %sampling frequency
Ts = 1/Fs;             % Sampling period       
LgS = 70;             % Length of signal, gate flux subsets
tgS = (0:LgS-1)*Ts;        % Time vector
fgS = Fs*(0:(LgS/2))/LgS; %frequency vector

gateFluxT=gateFlux26265+gateFlux26527+gateFlux27275+gateFluxS26(1:138);
fluxes=[gateFluxT(1:2:end-1) 0; 0 gateFluxT(2:2:end)];
%fluxes=[gateFluxS26(1:130); gateFlux26265(1:130); gateFlux26527(1:130); gateFlux27275(1:130); gateFluxT(1:130)];

clear fftFluxes2 pwFluxes2
for i=1:2 %depths
    %for j=1:5 %time
    pwFluxes2(i,:) = pwelch(fluxes(i,:),[],[],fgS,Fs);
        fftF=fft(fluxes(i,:)); 
        holdvar = abs(fftF/LgS);
        fftFluxes2(i,:) = holdvar(1:LgS/2+1);
        fftFluxes2(i,2:end-1) = 2*fftFluxes2(i,2:end-1);
    %end %for j, time
    
end %for i, depth

figure; loglog(fgS*86400,fftFluxes2(1,:)); hold all; loglog(fgS*86400,fftFluxes2(2,:));
figure; loglog(fgS*86400,pwFluxes2(1,:)); hold all; loglog(fgS*86400,pwFluxes2(2,:));

%% flux- decorrelation time, re-sectioning, proper spectra means with errorbars from std, not range

load('gateAdvection26527.mat')
load('gateAdvection26265.mat')
load('gateAdvection27275.mat')
load('gateAdvectionSurface.mat')

gateFluxT=gateFlux26265+gateFlux26527+gateFlux27275+gateFluxS26(1:138);

for i=1:130
    holdvar=corrcoef(gateFluxS26(1:139-i),gateFluxS26(i:138));
    autoCS(i)=holdvar(1,2);    
    holdvar=corrcoef(gateFlux26265(1:139-i),gateFlux26265(i:138));
    autoC26265(i)=holdvar(1,2);
    holdvar=corrcoef(gateFlux26527(1:139-i),gateFlux26527(i:138));
    autoC26527(i)=holdvar(1,2);
    holdvar=corrcoef(gateFlux27275(1:139-i),gateFlux27275(i:138));
    autoC27275(i)=holdvar(1,2);
    holdvar=corrcoef(gateFluxT(1:139-i),gateFluxT(i:138));
    autoCT(i)=holdvar(1,2);
end

figure; plot(0:129,autoCS);
hold all; plot(0:129,autoC26265);
plot(0:129,autoC26527); plot(0:129,autoC27275);
plot(0:129,autoCT);
legend('S26','26265','26527','27275','T')

%3 to 11 day decorrelation time

Fs=1/86400; %sampling frequency
Ts = 1/Fs;             % Sampling period       
LgS = 28;             % Length of signal, gate flux subsets
tgS = (0:LgS-1)*Ts;        % Time vector
fgS = Fs*(0:(LgS/2))/LgS; %frequency vector

fluxes=[gateFluxS26(1:138);gateFlux26265; gateFlux26527 ; gateFlux27275 ; gateFluxT];
%fluxes=[gateFluxS26(1:130); gateFlux26265(1:130); gateFlux26527(1:130); gateFlux27275(1:130); gateFluxT(1:130)];

clear fftFluxes
for i=[1 2 4 5] %depths (5 is sum), decorrelation time 3-4 days
    for j=1:23%11 %time
        fftF=fft(fluxes(i,5*j-4:23+5*j)); %11*j-10:17+11*j %1+26*(j-1):26*j%
        holdvar = abs(fftF/LgS);
        fftFluxes(i,j,:) = holdvar(1:LgS/2+1);
        fftFluxes(i,j,2:end-1) = 2*fftFluxes(i,j,2:end-1);
    end %for j, time
    
end %for i, depth

for i=3 %decorrelation time here is 11 days
    for j=1:11 %time
        fftF=fft(fluxes(i,11*j-10:17+11*j)); % %1+26*(j-1):26*j%
        holdvar = abs(fftF/LgS);
        fftFluxes(i,j,:) = holdvar(1:LgS/2+1);
        fftFluxes(i,j,2:end-1) = 2*fftFluxes(i,j,2:end-1);
    end %for j, time
    
end %for i, depth

fftFluxes(3,12:23,:)=nan;

spectraMeanF=squeeze(nanmean(fftFluxes,2));
%spectraMaxF=squeeze(max(fftFluxes,[],2));
%spectraMinF=squeeze(min(fftFluxes,[],2));
spectraStdF=squeeze(std(fftFluxes,0,2));
spectraStdF(3,:)=std(fftFluxes(3,1:11,:),0,2);

figure; 
errorbar(fgS*86400,spectraMeanF(1,:),spectraStdF(1,:))
hold all
for ii=2:5
    errorbar(fgS*86400,spectraMeanF(ii,:),spectraStdF(ii,:))
end %ii
title('fft spectra with errors')
set(gca,'YScale','log')
set(gca,'XScale','log')
xlabel('frequency in 1/day')
legend('Surface to \sigma=26','\sigma=26 to \sigma=26.5','\sigma=26.5 to \sigma=27','\sigma=27 to \sigma=27.5','sum')


for ii=1:5
    figure;
    loglog(fgS*86400,squeeze(fftFluxes(ii,:,:)));
    hold all
    errorbar(fgS*86400,spectraMeanF(ii,:),spectraStdF(ii,:),'linewidth',2)
set(gca,'YScale','log')
set(gca,'XScale','log')
xlabel('frequency in 1/day')
switch ii
    case 1 
        title('each spectra, gate flux surface to \sigma=26')
    case 2
        title('each spectra, gate flux \sigma=26 to \sigma=26.5')
    case 3
        title('each spectra, gate flux \sigma=26.5 to \sigma=27')
    case 4
        title('each spectra, gate flux \sigma=27 to \sigma=27.5')
    case 5
        title('each spectra, total gate flux')
end
end

%% velocity load
load('uva148levels14.mat', 'XC')
load('uva148levels14.mat', 'YC')
load('uva148levels14.mat', 'U')
u1=squeeze(U(331:332,160:161,10,:)); clear U
load('uva148levels14.mat', 'V')
v1=squeeze(V(331:332,160:161,10,:)); clear V
%% velocity calc
us=[squeeze(u1(1,1,:)).'; squeeze(u1(2,1,:)).'; squeeze(u1(1,2,:)).'; squeeze(u1(2,2,:)).'];
vs=[squeeze(v1(1,1,:)).'; squeeze(v1(2,1,:)).'; squeeze(v1(1,2,:)).'; squeeze(v1(2,2,:)).'];

for i=1:130
    holdvar=corrcoef(us(1,1:139-i),us(1,i:138));
    autoCU(1,i)=holdvar(1,2);    
    holdvar=corrcoef(us(2,1:139-i),us(2,i:138));
    autoCU(2,i)=holdvar(1,2);
    holdvar=corrcoef(us(3,1:139-i),us(3,i:138));
    autoCU(3,i)=holdvar(1,2);
    holdvar=corrcoef(us(4,1:139-i),us(4,i:138));
    autoCU(4,i)=holdvar(1,2);
    
    holdvar=corrcoef(vs(1,1:139-i),vs(1,i:138));
    autoCV(1,i)=holdvar(1,2);    
    holdvar=corrcoef(vs(2,1:139-i),vs(2,i:138));
    autoCV(2,i)=holdvar(1,2);
    holdvar=corrcoef(vs(3,1:139-i),vs(3,i:138));
    autoCV(3,i)=holdvar(1,2);
    holdvar=corrcoef(vs(4,1:139-i),vs(4,i:138));
    autoCV(4,i)=holdvar(1,2);
end

figure; plot(0:129,autoCU);
hold all; plot(0:129,autoCV);
legend('u','u','u','u','v','v','v','v')
grid on

%u decorrelation time is 15 days
%v decorrelation time is 9 days

Fs=1/86400; %sampling frequency
Ts = 1/Fs;             % Sampling period       
LgS = 28;             % Length of signal, gate flux subsets
tgS = (0:LgS-1)*Ts;        % Time vector
fgS = Fs*(0:(LgS/2))/LgS; %frequency vector
for i=1:4 %locations
    for j=1:14 %time
        if j<9
        pwUs(i,j,:) = pwelch(us(i,15*j-14:13+15*j),[],[],fgS,Fs);
        end
        pwVs(i,j,:) = pwelch(vs(i,9*j-8:19+9*j),[],[],fgS,Fs);
    end %for j, time
    
end %for i, location

spectraMeanUp=squeeze(mean(pwUs,2));
spectraStdUp=squeeze(std(pwUs,0,2));

spectraMeanVp=squeeze(mean(pwVs,2));
spectraStdVp=squeeze(std(pwVs,0,2));

 [PgU,FgU]=pwelch(u1,[],[],[],1/86400);
 [PgV,FgV]=pwelch(v1,[],[],[],1/86400);

%% velocity plots
figure; 
for i=1:2
    for j=1:2
        plot(1:148,squeeze(u1(i,j,:))); hold all; plot(1:148,squeeze(v1(i,j,:))); 
    end
end
xlabel('time (days)')
ylabel('speed'); legend('u','v','u','v','u','v','u','v'); title('Velocity at -4.4,36.4')

figure; loglog(FgU*86400,PgU); hold all; loglog(FgV*86400,PgV);
legend('u','v'); xlabel('Frequency in 1/day'); ylabel('Power')
title('Velocity Spectra at -4.4,36.4')
%% velocity plots
figure; 
errorbar(fgS*86400,spectraMeanUp(1,:),spectraStdUp(1,:))
hold all
for ii=2:4
    errorbar(fgS*86400,spectraMeanUp(ii,:),spectraStdUp(ii,:))
end %ii
title('Pwelch spectra with errors, U')
set(gca,'YScale','log')
set(gca,'XScale','log')
xlabel('frequency in 1/day')

figure; 
errorbar(fgS*86400,spectraMeanVp(1,:),spectraStdVp(1,:))
hold all
for ii=2:4
    errorbar(fgS*86400,spectraMeanVp(ii,:),spectraStdVp(ii,:))
end %ii
title('Pwelch spectra with errors, V')
set(gca,'YScale','log')
set(gca,'XScale','log')
xlabel('frequency in 1/day')
%% forcing: wind stress over time

load('C:\Users\JayB\Documents\MATLAB\MITgcm\addForceFull\windForEkman.mat')

tauN148=tauN(:,:,iNov12007:i148days);
tauE148=tauE(:,:,iNov12007:i148days);

tU=[squeeze(tauE148(16,22,:)).'; squeeze(tauE148(17,22,:)).'; squeeze(tauE148(16,23,:)).'; squeeze(tauE148(17,23,:)).'];
tV=[squeeze(tauN148(16,22,:)).'; squeeze(tauN148(17,22,:)).'; squeeze(tauN148(16,23,:)).'; squeeze(tauN148(17,23,:)).'];


for i=1:130
    holdvar=corrcoef(tU(1,1:139-i),tU(1,i:138));
    autoCtU(1,i)=holdvar(1,2);    
    holdvar=corrcoef(tU(2,1:139-i),tU(2,i:138));
    autoCtU(2,i)=holdvar(1,2);
    holdvar=corrcoef(tU(3,1:139-i),tU(3,i:138));
    autoCtU(3,i)=holdvar(1,2);
    holdvar=corrcoef(tU(4,1:139-i),tU(4,i:138));
    autoCtU(4,i)=holdvar(1,2);
    
    holdvar=corrcoef(tV(1,1:139-i),tV(1,i:138));
    autoCtV(1,i)=holdvar(1,2);    
    holdvar=corrcoef(tV(2,1:139-i),tV(2,i:138));
    autoCtV(2,i)=holdvar(1,2);
    holdvar=corrcoef(tV(3,1:139-i),tV(3,i:138));
    autoCtV(3,i)=holdvar(1,2);
    holdvar=corrcoef(tV(4,1:139-i),tV(4,i:138));
    autoCtV(4,i)=holdvar(1,2);
end

figure; plot(0:129,autoCtU);
hold all; plot(0:129,autoCtV);
legend('u','u','u','u','v','v','v','v')
grid on

%u decorrelation time is 12 timesteps, 3 days
%v decorrelation time is 9 timesteps, 2 days



Fs=4/86400; %sampling frequency
Ts = 1/Fs;             % Sampling period       
LgS = 112;             % Length of signal, gate flux subsets
tgS = (0:LgS-1)*Ts;        % Time vector
fgS = Fs*(0:(LgS/2))/LgS; %frequency vector
for i=1:4 %locations
    for j=1:41 %time
        pwTU(i,j,:) = pwelch(tU(i,12*j-11:100+12*j),[],[],fgS,Fs);
        pwTV(i,j,:) = pwelch(tV(i,9*j-8:103+9*j),[],[],fgS,Fs);
    end %for j, time
    
end %for i, location

spectraMeanTU=squeeze(mean(pwTU,2));
spectraStdTU=squeeze(std(pwTU,0,2));
%spectraMinTU=squeeze(min(pwTU,[],2));
spectraMeanTV=squeeze(mean(pwTV,2));
spectraStdTV=squeeze(std(pwTV,0,2));
%spectraMinTV=squeeze(min(pwTV,[],2));
 %[PgU,FgU]=pwelch(u1,[],[],[],1/86400);
 %[PgV,FgV]=pwelch(v1,[],[],[],1/86400);

 
 figure; plot(0:0.25:148,tU,'b'); hold on; plot(0:0.25:148,tV,'r'); legend('\tau E','\tau E','\tau E','\tau E','\tau N','\tau N','\tau N','\tau N')
 title('Timeseries of wind stresses','fontsize',16)
 xlabel('Simulation day','fontsize',16)
 ylabel('Velocity','fontsize',16)
 set(gca,'fontsize',14)
 
 
 figure; 
errorbar(fgS*86400,spectraMeanTU(1,:),spectraStdTU(1,:))
hold all
for ii=2:4
    errorbar(fgS*86400,spectraMeanTU(ii,:),spectraStdTU(ii,:))%,spectraMeanTU(ii,:)-spectraMinTU(ii,:),spectraMaxTU(ii,:)-spectraMeanTU(ii,:))
end %ii
title('Pwelch spectra with errors, eastward windtress','fontsize',16)
set(gca,'YScale','log')
set(gca,'XScale','log')
xlabel('frequency in 1/day','fontsize',16)
 set(gca,'fontsize',14)

figure; 
errorbar(fgS*86400,spectraMeanTV(1,:),spectraStdTV(1,:))%,spectraMeanTV(1,:)-spectraMinTV(1,:),spectraMaxTV(1,:)-spectraMeanTV(1,:))
hold all
for ii=2:4
    errorbar(fgS*86400,spectraMeanTV(ii,:),spectraStdTV(ii,:))%spectraMeanTV(ii,:)-spectraMinTV(ii,:),spectraMaxTV(ii,:)-spectraMeanTV(ii,:))
end %ii
title('Pwelch spectra with errors, northward windstress','fontsize',16)
set(gca,'YScale','log')
set(gca,'XScale','log')
xlabel('frequency in 1/day','fontsize',16)
 set(gca,'fontsize',14)
 
 
%  figure; hold on;
%  for i=1:4; for j=1:5
%      loglog(fgS*86400,squeeze(pwTU(i,j,:)))
%  end; end
%  set(gca,'YScale','log')
% set(gca,'XScale','log')

%% strait of gibraltar flow
load('transportGibraltarDaily.mat', 'inflow')

for i=1:130
    holdvar=corrcoef(inflow(1:149-i),inflow(i:148));
    %holdvar=corrcoef(inflow(1:end+1-i),inflow(i:end));
    autoCi(i)=holdvar(1,2);    
end

figure; plot(0:129,autoCi) %13 day decorrelation timescale

Fs=1/86400; %sampling frequency
Ts = 1/Fs;             % Sampling period       
LgS = 28;             % Length of signal, gate flux subsets
tgS = (0:LgS-1)*Ts;        % Time vector
fgS = Fs*(0:(LgS/2))/LgS; %frequency vector
for j=1:12 %time
    pwT(j,:) = pwelch(inflow(13*j-12:15+13*j),[],[],fgS,Fs);
end %for j, time


spectraMeanT=squeeze(mean(pwT,1));
spectraStdT=squeeze(std(pwT,0,1));
%spectraMinT=squeeze(min(pwT,[],1));

figure; errorbar(fgS*86400,spectraMeanT(1,:),spectraStdT(1,:))%MeanT(1,:)-spectraMinT(1,:),spectraMaxT(1,:)-spectraMeanT(1,:))
xlabel('frequency in 1/day','fontsize',16)
set(gca,'YScale','log')
set(gca,'XScale','log')
 set(gca,'fontsize',14)
 title('Gibraltar Inflow Spectrum','fontsize',16)
 
%% junk/notes

%pwelch example
[P1o,f1o]=pwelch(var1o);
figure
hold all
loglog(f1t,P1t)

%fft example
Fs = 1000;            % Sampling frequency                    
T = 1/Fs;             % Sampling period       
L = 1500;             % Length of signal
t = (0:L-1)*T;        % Time vector
Y = fft(X);
holdvar = abs(Y/L);
P1 = holdvar(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:(L/2))/L;
plot(f,P1) 
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')