%polyxpoly for intersections of manifolds
%polyarea for area of the lobes found

load('allmanifoldsDays9to149Resolved.mat')

[nt,~]=size(lontrFall)

lontrFall(lontrFall==0)=nan;
lontrBall(lontrBall==0)=nan;
lattrFall(lattrFall==0)=nan;
lattrBall(lattrBall==0)=nan;

for i=1:nt
    i
    [xi,yi,ii]=polyxpoly(lontrFall(i,:),lattrFall(i,:),lontrBall(i,:),lattrBall(i,:));
    np=length(xi);
    nint(i)=np;
    if nint(i)>0
    xint(i,1:np)=xi;
    yint(i,1:np)=yi;
    fiint(i,1:np)=ii(:,1);
    biint(i,1:np)=ii(:,2);
    end
end
disp('section 1 done')
%%
%load('lobesSurface2.mat')
plotting=false;
%clear lobAmany lobAfew nlobe compactness nClusters *Many *Few
for i=[1:74 76:nt];%1:nt
    i
    clear C nc
        if nint(i)>1
            if plotting
        figure; plot(lontrFall(i,:),lattrFall(i,:),'r')
        hold all; plot(lontrBall(i,:),lattrBall(i,:),'b')
        plot(lonCoast,latCoast,'k')
        plot(xint(i,1:nint(i)),yint(i,1:nint(i)),'go')
        axis([-6 -1 34.5 37.5])
        title(strcat('intersections day ',num2str(i+8)))
            end
            %for each day, I want to locate lobes/filaments. I want only 1
            %from each pair of cluster of intersections, the most compact
            %one, aka smallest perimeter:area
            %first, make clusters of intersections
            %next, for each pair of clusters, find lobes between each 
            %possible pair of intersections, calc area and perimeter,
            %choose one to keep
            [C,nc]=clusterDist(xint(i,1:nint(i)),yint(i,1:nint(i)),0.1);
            nClusters(i)=nc;
            if nc>1
            pair=0;
            
            for jc=1:nc-1
                jindices=find(C(jc,:)==1);
                for kc=jc+1:nc
                    pair=pair+1;
                    kindices=find(C(kc,:)==1);
                    intpair=0;
                    clear xlobes ylobes areas perimeters fiints biints
                    for jj=1:sum(C(jc,:))
                        j=jindices(jj);
                        for kk=1:sum(C(kc,:))
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
                xlobes(1:length(xlobe),intpair)=xlobe;
                ylobes(1:length(xlobe),intpair)=ylobe;
                if plotting
                    figure; plot(xlobe,ylobe)
                    title(strcat('possible lobe',num2str(pair),' intpair',num2str(intpair)))
                end
                areas(intpair)=polyarea(xlobe,ylobe);
                fiints(1:2,intpair)=[fiint(i,j); fiint(i,k)];
                biints(1:2,intpair)=[biint(i,j); biint(i,k)];
                %disp('checkpoint 1')

                        end
                    end
                
                    xlobes(xlobes==0)=nan;
                    ylobes(ylobes==0)=nan;
                perimeters=nansum(sqrt(diff(xlobes).^2 +diff(ylobes).^2));
                
                %add this??!!%%%%
                perimeters(perimeters==0)=inf;
                %%%%
                
                %areas=polyarea(xlobes,ylobes);
                mostCompact=find(perimeters./areas==min(perimeters./areas),1,'first');
                
                if plotting
                figure; plot(xlobes(:,mostCompact),ylobes(:,mostCompact)); title(strcat('Day ',num2str(i+8),' lobe ',num2str(pair)))
                end
                %xl(i,j)=length(xlobe);
                lobAmany(i,pair)=areas(mostCompact);
                compactness(i,pair)=min(perimeters./areas);
                fiintMany(i,pair,1:2)=fiints(:,mostCompact);
                biintMany(i,pair,1:2)=biints(:,mostCompact);
                %j=k+1;
                %disp('checkpoint 2')
                
                end
            end
            
          
                temp=sortrows([compactness(i,:).' lobAmany(i,:).' fiintMany(i,:,1).' fiintMany(i,:,2).' biintMany(i,:,1).' biintMany(i,:,2).'],-1);
                [tempA,indicesT]=unique(temp(:,2));
                tempA=tempA(tempA>0);
                
                nlobe(i)=length(tempA);
                lobAfew(i,1:nlobe(i))=tempA.';
                
                nII=length(temp(indicesT,3));
                fiintFew(i,1:nII,1)=temp(indicesT,3);
                fiintFew(i,1:nII,2)=temp(indicesT,4);
                biintFew(i,1:nII,1)=temp(indicesT,5);
                biintFew(i,1:nII,2)=temp(indicesT,6);
                %disp('checkpoint 3')
            end
        end 
        
 end
    
%figure; pcolor(lobAfew); shading 'flat'; colorbar
disp('section 2 done')
%%
%plot the lobes

nj=length(fiintFew(1,:,1));
nLobesHere=zeros([700 200]);
inside=zeros([700 200 69]);
for i=[1:74 76:nt]
    i
    if nClusters(i)>1
        
        
        
        for j=1:nj
            if fiintFew(i,j,1)>0
            if fiintFew(i,j,2)>fiintFew(i,j,1)
                    if biintFew(i,j,2)>biintFew(i,j,1)
                        xlobe=[lontrFall(i,fiintFew(i,j,1):fiintFew(i,j,2)) lontrBall(i,biintFew(i,j,2):-1:biintFew(i,j,1))];
                        ylobe=[lattrFall(i,fiintFew(i,j,1):fiintFew(i,j,2)) lattrBall(i,biintFew(i,j,2):-1:biintFew(i,j,1))];
                    else
                        xlobe=[lontrFall(i,fiintFew(i,j,1):fiintFew(i,j,2)) lontrBall(i,biintFew(i,j,2):biintFew(i,j,1))];
                        ylobe=[lattrFall(i,fiintFew(i,j,1):fiintFew(i,j,2)) lattrBall(i,biintFew(i,j,2):biintFew(i,j,1))];
                    end
            elseif biintFew(i,j,2)>biintFew(i,j,1)
                    xlobe=[lontrFall(i,fiintFew(i,j):-1:fiintFew(i,j,2)) lontrBall(i,biintFew(i,j,2):-1:biintFew(i,j,1))];
                    ylobe=[lattrFall(i,fiintFew(i,j):-1:fiintFew(i,j,2)) lattrBall(i,biintFew(i,j,2):-1:biintFew(i,j,1))];
            else
                    xlobe=[lontrFall(i,fiintFew(i,j):-1:fiintFew(i,j,2)) lontrBall(i,biintFew(i,j,2):biintFew(i,j,1))];
                    ylobe=[lattrFall(i,fiintFew(i,j):-1:fiintFew(i,j,2)) lattrBall(i,biintFew(i,j,2):biintFew(i,j,1))];
            end
            inside(:,:,i)=inside(:,:,i)+inpolygon(XC,YC,xlobe,ylobe);
            end
        end
        inlobe=(inside(:,:,i)>0);
        
%         figure; plot(XC(inlobe),YC(inlobe),'g.')
%         hold all; plot(lontrFall(i,:),lattrFall(i,:),'r')
%         plot(lontrBall(i,:),lattrBall(i,:),'b')
%         plot(lonCoast,latCoast,'k')
%         plot(xint(i,1:nint(i)),yint(i,1:nint(i)),'ko')
%         axis([-6 -1 34.5 37.5])
%         title(strcat('Lobe Extent, Simulation Day ',num2str(i+8)))
        
        nLobesHere=nLobesHere+inside(:,:,i);
    end
end

%%
disp('saving')
save('lobes8day148DailyVSurface.mat','-v7.3')
disp('done')
% nLobesPlot=nLobesHere;
% nLobesPlot(nLobesPlot==0)=nan;
% figure; pcolor(XC,YC,nLobesPlot); hold on; 
% plot(lonCoast,latCoast,'k')
% colorbar; shading 'flat'
% axis([-6 -1 34.5 37.5])
% title('Lobe Number')
%   
% figure; pcolor(XC,YC,nLobesHere./69); hold on; 
% plot(lonCoast,latCoast,'k')
% colorbar; shading 'flat'
% axis([-6 -1 34.5 37.5])
% title('Probability Map of a Lobe for a day')
        

