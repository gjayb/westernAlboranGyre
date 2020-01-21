%gathers manifolds, plots all on one, saves out

addpath('../mStuff')
%load('edgesWAGeuler.mat')
%size(isedge)

day0=14
isopycStr='Iso265'%'Iso2775'
% dayN=155

fn=strcat('wagBoundaryday',num2str(day0),isopycStr,'Int14.mat');%'Sigma'
load(fn,'lontrF');
load(fn,'lonCoast');
load(fn,'latCoast');
[nt,np]=size(lontrF)
%lontrF=zeros(dayN-day0+1,np);
%lattrF=lontrF;
%lontrB=lattrF;
%lattrB=lattrF;

disp('start loop')
for i=day0:146%[19:25 58:64 103:106 113:117 128:133]
i
    fn=strcat('wagBoundaryday',num2str(i),isopycStr,'Int14.mat');%'Sigma'
    load(fn,'lontrF');
    load(fn,'lontrB');
    load(fn,'lattrF');
    load(fn,'lattrB');
[nt2,np2]=size(lontrF)
[nt3,np3]=size(lontrB)
    lontrFall(i,1:np2)=lontrF(end,:);
    lontrBall(i,1:np3)=lontrB(end,:);
    lattrBall(i,1:np3)=lattrB(end,:);
    lattrFall(i,1:np2)=lattrF(end,:);
end
disp('end loop')
sizeF=size(lontrFall)

%figure; plot(lonCoast,latCoast,'k')
%hold on; plot(lontrFall,lattrFall,'r')
%plot(lontrBall,lattrBall,'b')
%plot(XC(logical(isedge(:,:,1))),YC(logical(isedge(:,:,1))),'g.')
%axis([-5.5 -1 34 38])
%title('All manifolds day 15 to 155')
%savefig('allmanifolds15to155D14.fig')

%figure; plot(lonCoast,latCoast,'k')
%hold on; %plot(lontrF,lattrF,'ro')
%plot(lontrB,lattrB,'bo')
%plot(XC(logical(isedge(:,:,1))),YC(logical(isedge(:,:,1))),'g.')
%axis([-5.5 -1 34 38])
%title('All backward manifolds Nov Dec 2007')
%save2pdf('allmanifoldsB1nd.pdf')

%figure; plot(lonCoast,latCoast,'k')
%hold on; plot(lontrF,lattrF,'ro')
%%plot(lontrB,lattrB,'bo')
%plot(XC(logical(isedge(:,:,1))),YC(logical(isedge(:,:,1))),'g.')
%axis([-5.5 -1 34 38])
%title('All forward manifolds Nov Dec 2007')
%save2pdf('allmanifoldsF1nd.pdf')

%close all
%disp('done figures')

fn=strcat('manifolds',isopycStr,'Int14NF.mat');%'Sigma'
%fn='manifoldsIso2675Int8NF.mat'
save(fn,'-v7.3')
