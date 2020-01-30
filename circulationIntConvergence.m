x=0:30;
y=0:30;
[xg,yg]=meshgrid(x,y);
a=2; b=-2;
u=a*yg;
v=b*xg;
r=10; theta1=0:0.1:2*pi; theta2=0:0.05:2*pi; theta3=0:0.01:2*pi;
xs1=11+r*cos(theta1); ys1=11+r*sin(theta1);
xs2=11+r*cos(theta2); ys2=11+r*sin(theta2);
xs3=11+r*cos(theta3); ys3=11+r*sin(theta3);

vorticityC=b-a;
vorticityInt=vorticityC*pi*r^2;

ds1=sqrt(diff(xs1).^2+diff(ys1).^2);
ds2=sqrt(diff(xs2).^2+diff(ys2).^2);
ds3=sqrt(diff(xs3).^2+diff(ys3).^2);

[circInt1,sign1,uds1]=circulationInt(xg,yg,u,xg,yg,v,xs1,ys1);
[circInt2,sign2,uds2]=circulationInt(xg,yg,u,xg,yg,v,xs2,ys2);
[circInt3,sign3,uds3]=circulationInt(xg,yg,u,xg,yg,v,xs3,ys3);