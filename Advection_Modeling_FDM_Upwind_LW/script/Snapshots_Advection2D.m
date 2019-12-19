dirname = 'data/Advection2D_Snapshots/';
files = dir(strcat([dirname, 'advection2D_*.txt']));
n = length(files);
Fig = figure('units','normalized','position',[0.1,0.1,0.5,0.75]);
filename = 'advection2D.gif';
for k = 1: n
    t = strsplit(files(k).name, '_');
    t = t(2);
    d = load(strcat([dirname,files(k).name]));
    cmap = load('data\brain.txt');
    colormap(cmap);
    imagesc(d);
    colorbar();
    caxis([0, 0.5]);
    xlabel('X [m]');
    ylabel('Y [m]')
    title(strcat([strtrim(t), ' s']));
    set(gca, 'yDir', 'normal', 'fontsize', 15, 'fontweight', 'bold');
    drawnow;
    
    frame = getframe(Fig); 
    im = frame2im(frame); 
    [imind,cm] = rgb2ind(im,256);
    if k == 1
        imwrite(imind,cm,filename,'gif','WriteMode','overwrite', ...
            'Loopcount',inf);
    else
        imwrite(imind,cm,filename,'gif','WriteMode','append', ...
            'DelayTime',0.01);
   end
end