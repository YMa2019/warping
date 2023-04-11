function plot_warping(t, warping, fig)
    figure (fig);clf;
    lsize = 16; % Label fontsize
    nsize = 18; % Axis fontsize
    plot(t, warping,'linewidth', 1);% plot the simulated warping function
    axis equal;
    ylim([0,1]);
    xlim([0,1]);
    xticks([0 0.2 0.4 0.6 0.8 1]);
    set(gca, 'Fontsize', nsize,'linewidth', 1.5)
    set(gcf,'paperpositionmode','auto');
    set(gcf,'windowstyle','normal');
    set(gca,'LooseInset',max(get(gca,'TightInset'), 0.02))
    set(gca,'fontweight','normal')
    opts.Colors     = get(groot,'defaultAxesColorOrder');
    opts.saveFolder = 'img/';
    opts.width      = 12;
    opts.height     = 10;
    opts.fontType   = 'Times';
end 