%% %%%%%% Manifolds %%%%%%%%
options=odeset('RelTol',10^(-6),'AbsTol',10^(-6));

RR=.25*10^4; 
xc=[2;6].*10^4; % for forward int from 35dt
% xc=[12.5;6].*10^4; % for forward int from 106dt   

theta1=(0:2*pi/100:2*pi); 

t1=linspace(dt*35,dt*71,10);
% t1=linspace(dt*106,dt*71,10);

x0=xc(1)+RR*cos(theta1); y0=xc(2)+RR*sin(theta1); 
x0=[x0 x0(1)]; y0=[y0 y0(1)];
    
% t1=linspace((it-1)*T1/4+3.*T1, (it-1)*T1/4, 10); % times for trajectory integration 
% x0=xc(1)+RR*cos(theta1); y0=xc(2)+RR*sin(theta1); 
% x0=[x0 x0(1)]; y0=[y0 y0(1)];

for i=1:length(t1)-1 % different integration time
        
    z0=[x0;y0];     
    tmesh1=[t1(i) (t1(i+1)+t1(i))/2 t1(i+1)];

        [~,zz]=ode45(@HamEqSolver_BiLin_Irina,tmesh1,z0(:),options,u,v,x,y,t);

        xunst=zz(:,1:2:end-1); 
        yunst=zz(:,2:2:end);
             
        %%%%% Resampling 1 %%%%%
        if i<length(t1)-1 %if not the last integration
           
            dlmin=RR/3; %10^(-3);

            xunst=[zz(end,1:2:end-1) zz(end,1)]; %end of previous integration, make a loop
            yunst=[zz(end,2:2:end) zz(end,2)];

            tunst=(1:1:length(xunst)); t1unst=(1:0.005:length(xunst));
            x1unst=interp1(tunst,xunst,t1unst); y1unst=interp1(tunst,yunst,t1unst); %interpolate to 200x as many points

            dxunst=diff(x1unst); dyunst=diff(y1unst); dlunst=sqrt(dxunst.^2+dyunst.^2); %distances between points

            nunst=[];

            ii=1; 
            while ii<(length(dlunst)-1) %while index ii is less than the number of points
                if dlunst(ii)>dlmin
                    nunst=[nunst ii]; ii=ii+1; %if this distance is bigger than the minimum, add ii to list nunst
                else 
                    jj=1;
                    while ((dlunst(ii)+dlunst(ii+jj)<dlmin)&((ii+jj)<length(dlunst)))
                          dlunst(ii)=dlunst(ii)+dlunst(ii+jj); jj=jj+1;  %otherwise, find how many points forward you can move until it gets too far
                     end
                    nunst=[nunst ii+jj]; ii=ii+jj; %add the next index, where the distance just goes over dlmin, to nunst
                end
                %ii
            end
            x0=x1unst(nunst); y0=y1unst(nunst); %set input to next integration to all the interpolated points that are needed to have separations no larger than dlmin
            %%%%%%%% end of resampling 1 %%%%%%%

        else
                 x0=zz(end,1:2:end-1);
                 y0=zz(end,2:2:end);
            
        end

                  
      i  
end
%%% modulus 
x0=x(1)+x0-floor((x0-x(1))./(x(end)-x(1)))*(x(end)-x(1));
y0=y(1)+y0-floor((y0-y(1))./(y(end)-y(1)))*(y(end)-y(1));
%%% end of modulus
figure(11); subplot(223); hold on; plot(x0,y0,'.b','MarkerSize',2); % axis([-1.5 1.5 -1 1]); title([num2str((it-1)/4) 'T']) 

