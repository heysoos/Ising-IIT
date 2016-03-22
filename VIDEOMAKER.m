figure(1)
filename = 'testgif.gif';
for n = 1:length(temp)
      imagesc(Corr_DTI(:,:,n))
      xlabel(['T = ', num2str(temp(n))])
      title('Ising Correlations')
      colorbar
      caxis([-1, 1])
      colormap(linspecer)
      drawnow
      frame = getframe(1);
      im = frame2im(frame);
      [imind,cm] = rgb2ind(im,256);
      if n == 1;
          imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
      else
          imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime', 0);
      end
end

vobj=VideoWriter('MyMovieFile', 'Motion JPEG AVI');
vobj.FrameRate=10;
vobj.Quality=75
open(vobj);
for i=500:550
    imagesc(Corr_DTI(:,:,i))
    xlabel(['T = ', num2str(temp(i))])
    title('Ising Correlations')
    colorbar
    caxis([-1, 1])
    colormap(linspecer)
%     drawnow
    F=getframe(gcf);
    writeVideo(vobj, F);
    cla(gca)
end
close(vobj)