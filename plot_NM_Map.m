function plot_NM_Map(data,clim,xlab,ylab,titleval,Colormap,CmapL)
    imAlpha = ones(size(data));
    imAlpha(isnan(data)) = 1;
    imagesc(data','AlphaData',imAlpha',clim);
    set(gca, 'YDir', 'normal');
    colormap(Colormap);
    xlabel(xlab);
    ylabel(ylab);
    title(titleval);
    set(gca,'FontSize',16);
    if ~isempty(CmapL)
        lcolorbar(CmapL);
    else
        c = colorbar;
        c.Ticks = linspace(clim(1),clim(2),5);
    end
    hold on
%     for i = 1:size(data,2)
%         xline(i+0.5);
%     end
%     for i = 1:size(data,1)
%         yline(i+0.5);
%     end
    plot(1:100,1:100,'--g','LineWidth',1);
    for i = 2:5
        plot(1:100,i * (1:100),'--r','LineWidth',1);
        plot(1:100,1/i * (1:100),'--r','LineWidth',1);
    end
    hold off
    set(gca,'TickLength',[0 0]);
    set(gca,'FontName','Ariel');
end