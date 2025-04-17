function NMC_plot(NMCoh,xlab,ylab,Xlim,Ylim)
    imagesc(NMCoh', [0 1]);
    colormap hot;
    set(gca, 'YDir', 'normal');   
    xlim(Xlim);
    ylim(Ylim);
    xlabel(xlab,'Fontsize',12);
    ylabel(ylab,'Fontsize',12);
end