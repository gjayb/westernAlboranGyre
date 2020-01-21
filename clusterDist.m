function [ C, nc ] = clusterDist( xin,yin,dist )
%CLUSTERDIST clusters input points (xin,yin) into groups within dist of
%each other
%   xin and yin must be vectors of the same length, and dist a positive scalar
%   returns C, cluster matrix size nc by length(xin), nc, scalar number of clusters

np=length(xin);
ingroup=zeros(1,np);
C=zeros(1,np);
row=0;

while sum(ingroup)<np
    row=row+1;
    p1=find(ingroup==0,1,'first');
    C(row,p1)=1;
    ingroup(p1)=1;
    outgroup=ones(1,np)-ingroup;
    %x=xin(outgroup); y=yin(outgroup);
    x1=xin(p1); y1=yin(p1);
    distances=sqrt((xin-x1).^2 + (yin-y1).^2);
    ptry=(distances<dist).*outgroup;
    C(row,ptry==1)=1;
    ingroup(ptry==1)=1;     
    %sum(ingroup)
end

C=C(1:row,:);
nc=row;

end

