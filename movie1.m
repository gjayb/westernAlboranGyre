% movie

%coriolis
%fc=2*7.2921e-5.*sind(YC);
%
v = VideoWriter('salinity35mAddforce.avi');
v.FrameRate=3;
open(v)
f1=figure('renderer','zbuffer');

for i=1:426
    clf
    pcolor(XC,YC,Sa(:,:,6,i)); shading 'flat'; caxis([36 38]); colorbar
%     pcolor(XC,YC,vorta(:,:,i)./fc); shading 'flat';caxis([-1e-4 1e-4]); colorbar
%     colormap(b2r(-1.5,1.5)); hold all;
%     contour(XC,YC,inWag2(:,:,1),[1 1],'m','linewidth',2)
     hold all; 
     plot(lonCoast,latCoast,'k','linewidth',2)
     %hold all
    %plot(xm(i,:),ym(i,:),'linewidth',2)
    title(num2str(i))
    %legend('coast','boundary')
    %axis([2e5 8e5 2e5 5e5])
    axis([-6 0 35 37.5])
        writeVideo(v,getframe(gcf))
        writeVideo(v,getframe(gcf))
    %plot(xtr(i,:),ytr(i,:),'--','linewidth',2)
    %legend('coast','boundary','yesterdays boundary advected')
    %axis([2e5 8e5 2e5 5e5])
%     subplot(3,1,1)
%     plot(vfluxs,'linewidth',2)
%     hold on
%     plot(i,vfluxs(i),'ro')
%     subplot(3,1,2)
%     plot(surfF,'linewidth',2)
%     hold on
%     plot(i,surfF(i),'ro')
%     subplot(3,1,3)
%     pcolor(lonEg.',latEg.',ep148(:,:,2*i-1)); shading 'flat'
%     colormap(b2r(-2e-8,5e-8))
%     colorbar
%     %caxis([0 1])
%     hold on
%     plot(lonCoast,latCoast,'k','linewidth',2)
%     axis([-6.5 -0.75 34.5 37.5])
%      title(num2str(i))
%      writeVideo(v,getframe(gcf))
%     writeVideo(v,getframe(gcf))
%     
%     clf
%     subplot(3,1,1)
%     plot(vfluxs,'linewidth',2)
%     hold on
%     plot(i,0.5*vfluxs(i)+0.5*vfluxs(i+1),'ro')
%     subplot(3,1,2)
%     plot(surfF,'linewidth',2)
%     hold on
%     plot(i,0.5*surfF(i)+0.5*surfF(i+1),'ro')
%     subplot(3,1,3)
%     pcolor(lonEg.',latEg.',ep148(:,:,2*i)); shading 'flat'
%     colormap(b2r(-2e-8,5e-8))
%     colorbar
%     %caxis([0 1])
%     hold on
%     plot(lonCoast,latCoast,'k','linewidth',2)
%     axis([-6.5 -0.75 34.5 37.5])
%      title(num2str(i))
%      writeVideo(v,getframe(gcf))
%     writeVideo(v,getframe(gcf))
    %pcolor(XC(:,75),-dBin(1:20),squeeze(inWag3(:,75,:,i)).'); shading 'flat'
    %axis([-6 -2 -220 0])
% contourf(XC,YC,SSHa(:,:,i),-1:0.05:1)
% hold on; contour(XC,YC,inWagH(:,:,1),[1 1],'m--','linewidth',2)
% %plot(center1(1),center1(2),'m*')
% contour(XC,YC,inWAGS(:,:,i),[1 1],'b','linewidth',2)
% %plot(centerSsh(1,i),centerSsh(2,i),'bo')
% vel1=velEastward(:,:,i); vel2=velNorthward(:,:,i);
% index2=abs(vel1)>0; index3=find(index2==1);
% x1=downsample(XC(index3),5); y1=downsample(YC(index3),5); v1=downsample(vel1(index3),5); v2=downsample(vel2(index3),5);
% quiver(x1,y1,v1,v2,'k','linewidth',2)
% colormap(b2r(-0.25, 0.25))
% %     plot(lonCoast,latCoast,'k')
% %     hold on
% %     plot(lontrFall(i,:),lattrFall(i,:),'r')
% %     plot(lontrBall(i,:),lattrBall(i,:),'b')
%     %pcolor(XC,YC,SSHa(:,:,i)); shading 'flat'; caxis([-0.25 0.25])
%     %hold on
%     %quiver(XC(1:5:end,1:5:end),YC(1:5:end,1:5:end),Urot(1:5:end,1:5:end,i),Vrot(1:5:end,1:5:end,i));
     
%     legend('SSH','euler WAG','lagrange WAG','euler WAG edge velocity')
%     writeVideo(v,getframe(gcf))
% %     plot(lontrFallB(i,:),lattrFallB(i,:),'m')
% %     plot(lontrBallB(i,:),lattrBallB(i,:),'c')
    
    writeVideo(v,getframe(gcf))
    %writeVideo(v,getframe(gcf))
    %writeVideo(v,getframe(gcf))
end
close(v)