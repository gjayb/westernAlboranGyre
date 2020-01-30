function [ integral1,vec1 ] = crossCurveFlux( xm,ym,u,v,xs,ys,field,depth )
%crossCurveFlux takes in xm,ym the grid in meters for velocity field u,v in
%the east-north directions, and xs,ys the closed curve around which to
%integrate, and field, the grid of values of the thing being fluxed, 
% and depth, the vertical extent on same grid as field. 
% the u,v are interpolated onto the midpoints of xs,ys, dotted
%with <dy,-dx>, and field*(vel perpendicular)*(segment length)*depth is
% integrated around the curve, giving integral1. 
%xm,ym are expected to be nx by ny
%u,v,field,depth are expected to be nx by ny by nt
%xs,ys are expected to be 1d vectors of the same length
n=length(xs);
n2=length(ys);
if n~=n2
    disp('error in size of xs,ys')
end
[nx,ny]=size(xm);
[nx2,ny2]=size(ym);
[nx3,ny3,nt1]=size(u);
[nx4,ny4,nt2]=size(v);
[nx5,ny5,nt3]=size(field);
[nx6,ny6,nt4]=size(depth);
if length(unique([nx nx2 nx3 nx4 nx5 nx6]))>1
    disp('error in dimension 1 length of xm,ym,u,v,field,depth')
end
if length(unique([ny ny2 ny3 ny4 ny5 ny6]))>1
    disp('error in dimension 2 length of xm,ym,u,v,field,depth')
end
if length(unique([nt1 nt2 nt3 nt4]))>1
    disp('error in dimension 3 length of u,v,field,depth')
end

%check if xs,ys 1==end
if xs(1)~=xs(n) || ys(1)~=ys(n)
    xs(n+1)=xs(1);
    ys(n+1)=ys(1);
end
%centerpoints of segments, dx,dy
xc=0.5*xs(1:end-1)+0.5*xs(2:end);
yc=0.5*ys(1:end-1)+0.5*ys(2:end);
dx=xs(2:end)-xs(1:end-1);
dy=ys(2:end)-ys(1:end-1);
n=length(xc);

%interpolate u,v, depth,field
uc=zeros(n,nt1);
vc=uc;
depths=uc;
fields=uc;
for i=1:nt1
uc(1:n,i)=griddata(xm,ym,u(:,:,i),xc,yc);
vc(1:n,i)=griddata(xm,ym,v(:,:,i),xc,yc);
depths(1:n,i)=griddata(xm,ym,depth(:,:,i),xc,yc);
fields(1:n,i)=griddata(xm,ym,field(:,:,i),xc,yc);
end
%vectors for dot product to get direction, integrate
integral1=zeros(nt1,1);
vecPerp=[dy(:).'; -dx(:).'];
vec1=zeros(n,nt1);
for i=1:nt1
    uvec=[uc(1:n,i).';vc(1:n,i).'];
    uds=dot(uvec,vecPerp);
    sizeU=size(uds)
    sizeD=size(depths(:,i))
    sizeF=size(fields(:,i))
    integral1(i)=nansum(uds.'.*depths(:,i).*fields(:,i));
    vec1(1:n,i)=uds.'.*depths(:,i).*fields(:,i);
end


end

