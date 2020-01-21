times=8640:8640:(148*8640)%15*8640
Ua=rdmds('Uave',times);
U=mean(Ua,4); clear Ua
Va=rdmds('Vave',times);
V=mean(Va,4); clear Va
save('uvNativeAve148.mat')
%SSH=rdmds('SSHave',times);
%save('checkrunAFF.mat','-v7.3')