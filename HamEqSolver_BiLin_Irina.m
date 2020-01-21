function dzdt = HamEqSolver_BiLin(S,z,U,V,xarray,yarray,tarray)
% S

DXvel= xarray(2)-xarray(1); 
Xmin= min(xarray); 
DYvel= yarray(2)-yarray(1); 
Ymin= min(yarray); 
DTvel=tarray(2)-tarray(1); 

n=length(z);
dzdt=zeros(n,1);

zz1=z(1:2:end); 
zz2=z(2:2:end); 

vx=[]; vy=[];

% %%%%%%%%%%%% time-dependent
Tmin=min(tarray);
int1=floor((S-Tmin)/DTvel)+1; %int1=1;
int2=int1+1;

for i=1:n/2
    Y=[zz1(i) zz2(i)];
    inx1=floor((Y(1)-Xmin)/DXvel)+1; %! +1 if starting with x(1)=0
    iny1=floor((Y(2)-Ymin)/DYvel)+1;
    inx2=inx1+1;
    iny2=iny1+1;

    if (inx1<1)||(inx2>length(xarray))||(iny1<1)||(iny2>length(yarray))||isnan(inx1)||isnan(inx2)||isnan(iny1)||isnan(iny2)
        u0=0; v0=0;
        %print('out of bounds point')
    else

%CCC at t1
%CC linear interp at y1
	U1_t1=((U(iny1,inx2,int1)-U(iny1,inx1,int1))*Y(1)+...
     U(iny1,inx1,int1)*(Xmin+DXvel*(inx2-1))-...
     U(iny1,inx2,int1)*(Xmin+DXvel*(inx1-1)))/DXvel;

	V1_t1=((V(iny1,inx2,int1)-V(iny1,inx1,int1))*Y(1)+...
     V(iny1,inx1,int1)*(Xmin+DXvel*(inx2-1))-...
     V(iny1,inx2,int1)*(Xmin+DXvel*(inx1-1)))/DXvel;

%CC linear interp at y2
	U2_t1=((U(iny2,inx2,int1)-U(iny2,inx1,int1))*Y(1)+...
     U(iny2,inx1,int1)*(Xmin+DXvel*(inx2-1))-...
     U(iny2,inx2,int1)*(Xmin+DXvel*(inx1-1)))/DXvel;

	V2_t1=((V(iny2,inx2,int1)-V(iny2,inx1,int1))*Y(1)+...
     V(iny2,inx1,int1)*(Xmin+DXvel*(inx2-1))-...
     V(iny2,inx2,int1)*(Xmin+DXvel*(inx1-1)))/DXvel;

%CC linear interp at x
	u0_t1=((U2_t1-U1_t1)*Y(2)+U1_t1*(Ymin+DYvel*(iny2-1))-...
     U2_t1*(Ymin+DYvel*(iny1-1)))/DYvel;
	v0_t1=((V2_t1-V1_t1)*Y(2)+V1_t1*(Ymin+DYvel*(iny2-1))-...
     V2_t1*(Ymin+DYvel*(iny1-1)))/DYvel;
%CCC at t2
%CC linear interp at y1
	U1_t2=((U(iny1,inx2,int2)-U(iny1,inx1,int2))*Y(1)+...
     U(iny1,inx1,int2)*(Xmin+DXvel*(inx2-1))-...
     U(iny1,inx2,int2)*(Xmin+DXvel*(inx1-1)))/DXvel;

	V1_t2=((V(iny1,inx2,int2)-V(iny1,inx1,int2))*Y(1)+...
     V(iny1,inx1,int2)*(Xmin+DXvel*(inx2-1))-...
     V(iny1,inx2,int2)*(Xmin+DXvel*(inx1-1)))/DXvel;

%CC linear interp at y2
	U2_t2=((U(iny2,inx2,int2)-U(iny2,inx1,int2))*Y(1)+...
     U(iny2,inx1,int2)*(Xmin+DXvel*(inx2-1))-...
     U(iny2,inx2,int2)*(Xmin+DXvel*(inx1-1)))/DXvel;

	V2_t2=((V(iny2,inx2,int2)-V(iny2,inx1,int2))*Y(1)+...
     V(iny2,inx1,int2)*(Xmin+DXvel*(inx2-1))-...
     V(iny2,inx2,int2)*(Xmin+DXvel*(inx1-1)))/DXvel;

%CC linear interp at x
	u0_t2=((U2_t2-U1_t2)*Y(2)+U1_t2*(Ymin+DYvel*(iny2-1))-...
     U2_t2*(Ymin+DYvel*(iny1-1)))/DYvel;
	v0_t2=((V2_t2-V1_t2)*Y(2)+V1_t2*(Ymin+DYvel*(iny2-1))-...
     V2_t2*(Ymin+DYvel*(iny1-1)))/DYvel; 
 
%CCC interpolation in t between t1 and t2
	u0=((u0_t2-u0_t1)*S+u0_t1*(Tmin+DTvel*(int2-1))-u0_t2*(Tmin+DTvel*(int1-1)))/DTvel;
	v0=((v0_t2-v0_t1)*S+v0_t1*(Tmin+DTvel*(int2-1))-v0_t2*(Tmin+DTvel*(int1-1)))/DTvel;
    end
    
     vx=[vx u0];
     vy=[vy v0];

end

dzdt=[vx'; vy'];

dzdt=reshape(dzdt,[n/2 2])';
dzdt=dzdt(:);