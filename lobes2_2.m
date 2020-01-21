%find lobes between adjacent clusters
%polyxpoly for intersections of manifolds
%polyarea for area of the lobes found

%load('allmanifoldsDays9to149Resolved.mat')
addpath('/nobackup1/gbrett/mStuff')
load('geometrySpinupSteady.mat')
load('manifoldsIso265int14.mat')
[nt,~]=size(lontrFall)

lontrFall(lontrFall==0)=nan;
lontrBall(lontrBall==0)=nan;
lattrFall(lattrFall==0)=nan;
lattrBall(lattrBall==0)=nan;
%%
for i=114:120
    i
    [xi,yi,ii]=polyxpoly(lontrFall(i,:),lattrFall(i,:),lontrBall(i,:),lattrBall(i,:));
    np=length(xi)
    nint(i)=np;
    if nint(i)>0
    xint(i,1:np)=xi;
    yint(i,1:np)=yi;
    fiint(i,1:np)=ii(:,1);
    biint(i,1:np)=ii(:,2);
    end
  clear xi yi ii
end
disp('section 1 done')
%%
%load('lobes8day148DailyVSurfaceClean.mat')
plotting=false;
%clear lobAmany lobAfew nlobe compactness perimeters areas nClusters *Many *Few xlobe* ylobe*
xmin=min(min(lontrFall));
ymin=min(min(lattrFall));
%notaBene='skipping days 81 to 88'
for i=114:120%nt%[1:80 89:nt];%[1:74 76:nt];%1:nt
    i
    clear C nc
        if nint(i)>1 & nint(i)<2500
            if plotting
        figure; plot(lontrFall(i,:),lattrFall(i,:),'r')
        hold all; plot(lontrBall(i,:),lattrBall(i,:),'b')
        plot(lonCoast,latCoast,'k')
        plot(xint(i,1:nint(i)),yint(i,1:nint(i)),'go')
        axis([-6 -1 34.5 37.5])
        title(strcat('intersections day ',num2str(i+8)))
            end
            [C,nc]=clusterDist(xint(i,1:nint(i)),yint(i,1:nint(i)),0.1); 
            nc
		nClusters(i)=nc;
            
            if nc>1 & nc<13
                clear xc
                for j=1:nc
                    xc(j)=mean(xint(logical(C(j,:))));                    
                end %for j=1:nc
                disp('xc made')
                clusterorder=sortrows([ xc' (1:nc)']);
                pair=0;
                for jCluster=2:nc
                    clear areas perimeters xlobe* ylobe*
                    jc=clusterorder(jCluster,2);
                    jindices=find(C(jc,:)==1);
                    pair=pair+1;
                        jcN=clusterorder(jCluster-1,2);
                        kindices=find(C(jcN,:)==1);
                        intpair=0;
                            for jj=1:length(jindices)
                                            j=jindices(jj);
                              for kk=1:length(kindices)
                                                k=kindices(kk);
                                                intpair=intpair+1;
                                                %disp('checkpoint 0')
                                    if fiint(i,k)>fiint(i,j)
                                        if biint(i,k)>biint(i,j)
                                            xlobe=[lontrFall(i,fiint(i,j):fiint(i,k)) lontrBall(i,biint(i,k):-1:biint(i,j))];
                                            ylobe=[lattrFall(i,fiint(i,j):fiint(i,k)) lattrBall(i,biint(i,k):-1:biint(i,j))];
                                        else
                                            xlobe=[lontrFall(i,fiint(i,j):fiint(i,k)) lontrBall(i,biint(i,k):biint(i,j))];
                                            ylobe=[lattrFall(i,fiint(i,j):fiint(i,k)) lattrBall(i,biint(i,k):biint(i,j))];
                                        end
                                    elseif biint(i,k)>biint(i,j)
                                        xlobe=[lontrFall(i,fiint(i,j):-1:fiint(i,k)) lontrBall(i,biint(i,k):-1:biint(i,j))];
                                        ylobe=[lattrFall(i,fiint(i,j):-1:fiint(i,k)) lattrBall(i,biint(i,k):-1:biint(i,j))];
                                    else
                                        xlobe=[lontrFall(i,fiint(i,j):-1:fiint(i,k)) lontrBall(i,biint(i,k):biint(i,j))];
                                        ylobe=[lattrFall(i,fiint(i,j):-1:fiint(i,k)) lattrBall(i,biint(i,k):biint(i,j))];
                                    end
                                    %meters!
                                    xlobem=111000*cosd(ylobe).*(xlobe-xmin*ones(size(xlobe)));
                                    ylobem=111000*(ylobe-ymin*ones(size(ylobe)));
        		%clear xlobe ylobe            
	                xlobes(1:length(xlobe),intpair)=xlobem;
                                    ylobes(1:length(xlobe),intpair)=ylobem;
                                    if plotting
                                        figure; plot(xlobe,ylobe)
                                        title(strcat('possible lobe',num2str(pair),' intpair',num2str(intpair)))
                                    end

                                    areas(intpair)=polyarea(xlobem,ylobem);
                                    if areas(intpair)<1
                                        areas(intpair)=0;
                                    end
                                    fiints(1:2,intpair)=[fiint(i,j); fiint(i,k)];
                                    biints(1:2,intpair)=[biint(i,j); biint(i,k)];
                                    %disp('checkpoint 1')

                              end%kk
                            end%jj
                    xlobes(xlobes==0)=nan;
                    ylobes(ylobes==0)=nan;
                    perimeters=nansum(sqrt(diff(xlobes).^2 +diff(ylobes).^2));
                    perimeters(perimeters==0)=inf;
                    mostCompact=find(perimeters./areas==min(perimeters./areas),1,'first');
                    lobAfew(i,pair)=areas(mostCompact);
                    compactness(i,pair)=min(perimeters./areas);
                    fiintFew(i,pair,1:2)=fiints(:,mostCompact);
                    biintFew(i,pair,1:2)=biints(:,mostCompact);
                clear fiints biints 
                end%for jCluster=2:nc
                
            end%if nc>1
        else
		nClusters(i)=0;        
        end%if nint(i)>1
        
end
disp('checkpoint')
[length1,~]=size(compactness)
if length1<i
	compactness(i)=0;
	fiintFew(i,1,1:2)=0;
	biintFew(i,1,1:2)=0;
end

disp('section 2 done')
%%
%plot the lobes

clear xlobe ylobe inlobe %nj nLobesHere inside 

nj=length(fiintFew(1,:,1))
%nLobesHere=zeros([700 200]);
%inside=zeros([700 200 nt]);
%for i=113:131%nt %[1:74 76:nt]
%    i
%    if nClusters(i)>1
%        
%        
%        
%        for j=1:nj
%            
%            if fiintFew(i,j,1)>0 & compactness(i,j)<0.1
%            if fiintFew(i,j,2)>fiintFew(i,j,1)
%                    if biintFew(i,j,2)>biintFew(i,j,1)
%                        xlobe=[lontrFall(i,fiintFew(i,j,1):fiintFew(i,j,2)) lontrBall(i,biintFew(i,j,2):-1:biintFew(i,j,1))];
%                        ylobe=[lattrFall(i,fiintFew(i,j,1):fiintFew(i,j,2)) lattrBall(i,biintFew(i,j,2):-1:biintFew(i,j,1))];
%                    else
%                        xlobe=[lontrFall(i,fiintFew(i,j,1):fiintFew(i,j,2)) lontrBall(i,biintFew(i,j,2):biintFew(i,j,1))];
%                        ylobe=[lattrFall(i,fiintFew(i,j,1):fiintFew(i,j,2)) lattrBall(i,biintFew(i,j,2):biintFew(i,j,1))];
%                    end
%            elseif biintFew(i,j,2)>biintFew(i,j,1)
%                    xlobe=[lontrFall(i,fiintFew(i,j):-1:fiintFew(i,j,2)) lontrBall(i,biintFew(i,j,2):-1:biintFew(i,j,1))];
%                    ylobe=[lattrFall(i,fiintFew(i,j):-1:fiintFew(i,j,2)) lattrBall(i,biintFew(i,j,2):-1:biintFew(i,j,1))];
%            else
%                    xlobe=[lontrFall(i,fiintFew(i,j):-1:fiintFew(i,j,2)) lontrBall(i,biintFew(i,j,2):biintFew(i,j,1))];
%                    ylobe=[lattrFall(i,fiintFew(i,j):-1:fiintFew(i,j,2)) lattrBall(i,biintFew(i,j,2):biintFew(i,j,1))];
%            end
%            inside(:,:,i)=inside(:,:,i)+inpolygon(XC,YC,xlobe,ylobe);
%            end
%        end %for j
%        inlobe=(inside(:,:,i)>0);
        
%        figure;%(3); clf
%        plot(XC(inlobe),YC(inlobe),'g.')
%        hold all; plot(lontrFall(i,:),lattrFall(i,:),'r')
%        plot(lontrBall(i,:),lattrBall(i,:),'b')
%        plot(lonCoast,latCoast,'k')
%        plot(xint(i,1:nint(i)),yint(i,1:nint(i)),'ko')
%        axis([-6 -1 34.5 37.5])
%        title(strcat('Lobe Extent, Simulation Day ',num2str(i+8)))
        
%        nLobesHere=nLobesHere+inside(:,:,i);
%    end
%end
disp('section 3 done')
%%
%nLobesPlot=nLobesHere;
%nLobesPlot(nLobesPlot==0)=nan;
disp('saving')
save('lobes14dayIso265days114to120.mat','-v7.3')
disp('done')

% 
% figure; pcolor(XC,YC,nLobesPlot); hold on; 
% plot(lonCoast,latCoast,'k')
% colorbar; shading 'flat'
% axis([-6 -1 34.5 37.5])
% title('Lobe Number')
%   
% figure; pcolor(XC,YC,nLobesHere./nt); hold on; 
% plot(lonCoast,latCoast,'k')
% colorbar; shading 'flat'
% axis([-6 -1 34.5 37.5])
% title('Probability Map of a Lobe for a day')
       
