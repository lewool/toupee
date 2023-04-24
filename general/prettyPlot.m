function prettyPlot(gca, tl)

box off;
set(gca,'tickdir','out');
if nargin < 2
    set(gca,'TickLength', [.01 .01])
else
    set(gca,'TickLength', [tl tl])
end
ax = gca;
ax.YColor = [0 0 0];
ax.XColor = [0 0 0];
set(0,'DefaultAxesFontName','Arial')