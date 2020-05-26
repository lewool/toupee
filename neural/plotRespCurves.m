function plotRespCurves(eventWindow, meanResp, semResp, respColors,ls)

for m = 1:size(meanResp,1)
    plotMean = plot(eventWindow,(meanResp(m,:)),'linestyle',ls,'linewidth',2,'color',respColors(m,:));
    hold on;
    plotSEM = fill([eventWindow';flipud(eventWindow')],[(meanResp(m,:)-semResp(m,:))';flipud((meanResp(m,:)+semResp(m,:))')],respColors(m,:), 'LineStyle', 'none');
    alpha(0.2);
end