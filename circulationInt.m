function [ integral1, sign1,uds ] = circulationInt( xum,yum,u,xvm,yvm,v,xs,ys )
%circulationInt takes in xm,ym the grid in meters for vector field u,v in
%the east-north directions, and xs,ys the closed curve around which to
%integrate. the u,v are interpolated onto the midpoints of xs,ys, dotted
%with <dx,dy>, and integrated around the curve, giving integral1. sign1 is
%the sum of the sign of the cross product of each segment on xs,ys with the
%next, meaning it will be positive for counterclockwise closed loops and
%negative for clockwise closed loops.
%xm,ym, should be nx by ny
%u,v should be nx by ny by nt
%xs ys should be 1d vectors
%NB: updating for multi-time calcs
n=length(xs);
n2=length(ys);
if n~=n2
    disp('error in size of xs,ys')
end
[nx,ny]=size(xum);
[nx2,ny2]=size(yum);

[nxa,nya]=size(xvm);
[nx2a,ny2a]=size(yvm);

[nx3,ny3,nt1]=size(u);
[nx4,ny4,nt2]=size(v);

if length(unique([nx nxa nx2a nx2 nx3 nx4]))>1
    disp('error in dimension 1 length of xu,xv,yu,yv,u,v')
end
if length(unique([ny nya ny2a ny2 ny3 ny4]))>1
    disp('error in dimension 2 length of xu,xv,yu,yv,u,v')
end
if length(unique([nt1 nt2]))>1
    disp('error in dimension 3 length of u,v')
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
%create sign1
crosses=cross([dx(1:end-1).' dy(1:end-1).' zeros(size(dx(1:end-1))).'],[dx(2:end).' dy(2:end).' zeros(size(dx(2:end))).']);
sign1=sum(sign(crosses(3,:)));
%interpolate u,v
uc=zeros(n,nt1);
vc=uc;
for i=1:nt1
uc(1:n,i)=griddata(xum,yum,u(:,:,i),xc,yc,'cubic');
vc(1:n,i)=griddata(xvm,yvm,v(:,:,i),xc,yc,'cubic');
end
%vectors for dot product to get direction, integrate
dsvec=[dx(:).';dy(:).'];
integral1=zeros(nt1,1);
for i=1:nt1
    uvec=[uc(1:n,i).';vc(1:n,i).'];
    uds=dot(uvec,dsvec);
    integral1(i)=nansum(uds);
end
end

