function prettyPlot(gca)

box off;
set(gca,'tickdir','out');
set(gca,'TickLength', [.005 .005])
ax = gca;
ax.YColor = [0 0 0];
ax.XColor = [0 0 0];
set(0,'DefaultAxesFontName','Arial')