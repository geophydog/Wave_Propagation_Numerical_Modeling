n = 9;
figure(1);
for i = 1: n
    subplot(3, 3, i)
    colormap(bone);
    imagesc(x, z, reshape(u(floor(i*10+40), :), [nz, nz])');
    colorbar();
    hold on;
    scatter(sx*dx, sz*dz, 100, 'filled', 'o', 'markerfacecolor', 'g', ...
        'markeredgecolor', 'r', 'linewidth', 1.5);
    hold on;
    xlabel('X [m]');
    ylabel('Z [m]');
    strtime = strcat( [num2str((i*10+40)*dt), ' s']);
    title(strtime);
    set(gca, 'fontsize', 12, 'fontweight', 'bold', 'ydir', 'normal');
end