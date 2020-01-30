%load('allmanifoldsDays9to149Resolved.mat')

xhist=-6:0.05:-2;
yhist=35:0.05:36.9;

[xh,yh]=meshgrid(xhist,yhist);



isoStrs=[2635 265 2675 27 2725 275 2775 28 2825 285 2875 289 29 291];
for iiso=1:7
    isoStr=num2str(isoStrs(iiso))
    fnl=strcat('manifoldsSigma',isoStr,'Int8.mat');
    load(fnl);
[nt,~]=size(lontrFall);
heatF=zeros(size(xh));
heatB=heatF;
heat=heatF;
%%
for t=1:nt %for each time
    t
    XF=lontrFall(t,:);
    YF=lattrFall(t,:);
    XB=lontrBall(t,:);
    YB=lattrBall(t,:);
    for i=1:length(xhist)-1 %for each longitude bin
        i;
        XV=xhist([i i+1 i+1 i]);
        for j=1:length(yhist)-1 %for each latitude bin
            YV=yhist([j j j+1 j+1]);
            h1=sum(inpolygon(XF,YF,XV,YV));
            h2=sum(inpolygon(XB,YB,XV,YV));
            heatF(j,i) = heatF(j,i)+sign(h1);%add 1 if the manifold goes through the box
            heatB(j,i) = heatB(j,i)+sign(h2);
            heat(j,i) = heat(j,i)+sign(h1+h2);
        end
    end
end

pmanF=heatF./nt;
pmanB=heatB./nt;
pman=heat./nt;

fns=strcat('pmanSigma',isoStr,'Int8.mat');
save(fns,'heat*','pman*','xh','yh')

end %iiso

%% plots

% figure; pcolor(xh,yh,pmanF); shading 'flat'; colorbar
% title('Probability of Forward Manifold')
% hold on; plot(lonCoast,latCoast,'k')
% 
% figure; pcolor(xh,yh,pmanB); shading 'flat'; colorbar
% title('Probability of Backward Manifold')
% hold on; plot(lonCoast,latCoast,'k')
% 
% figure; pcolor(xh,yh,pman); shading 'flat'; colorbar
% title('Probability of Any Manifold')
% hold on; plot(lonCoast,latCoast,'k')

%% find area of the core

%figure; [c2,h2]=contour(xh,yh,pman,[0 0]);%[0.05 0.05],'k');
%surface 8day wide
%length1=c2(2,1);
%pman05x=c2(1,2:1+length1);
%pman05y=c2(2,2:1+length1);

%sigma=26 10-day
% pman0x=c2(1,198:end);
% pman0y=c2(2,198:end);

%surface narrow 8day
% pman05x=c2(1,1044:1044+204);
% pman05y=c2(2,1044:1044+204);
% figure; plot(pman05x,pman05y)

%surface 14day wide
% pman05x=c2(1,1276:1276+184);
% pman05y=c2(2,1276:1276+184);
% figure; plot(pman05x,pman05y)

%surface 8day wide
%[xi,yi,ii]=polyxpoly(pman05x,pman05y,lonCoast,latCoast);
%iintersec=find(xi<-3 & xi>-5.5 & yi<35.4);
%corex=cat(1,(pman05x(ii(iintersec(1),1):-1:ii(iintersec(2),1))).', lonCoast(ii(iintersec(2),2):-1:ii(iintersec(1),2))); %this is not automatic!
%corey=cat(1,(pman05y(ii(iintersec(1),1):-1:ii(iintersec(2),1))).', latCoast(ii(iintersec(2),2):-1:ii(iintersec(1),2)));

%surface narrow 8day, wide 14day
%  corex=pman05x;
%  corey=pman05y;
% 
% corexm=(corex-min(corex)*ones(size(corex))).*111000.*cosd(corey);
% coreym=(corey-min(corey)*ones(size(corey))).*111000;
% coreArea=polyarea(corexm,coreym);

%corexm=(corex-min(min(XC))*ones(size(corex))).*111000.*cosd(corey);
%coreym=(corey-min(min(YC))*ones(size(corey))).*111000;
%surface wide 8day zero probability core:
%figure; [c2,h2]=contour(xh,yh,pman,[0.0 0.0],'k');
%corex=c2(1,663:704);corey=c2(2,663:704);