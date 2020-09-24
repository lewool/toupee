function plotPSTHs(eventWindow, meanResp, semResp, contrastColors,ls)

for m = 1:size(meanResp,1)
    meanResp(m,:) = smooth(meanResp(m,:),1);
    plotMean = plot(eventWindow,(meanResp(m,:)),'linestyle',ls,'linewidth',2,'color',contrastColors(m,:));
    hold on;
    plotSEM = fill([eventWindow';flipud(eventWindow')],[(meanResp(m,:)-semResp(m,:))';flipud((meanResp(m,:)+semResp(m,:))')],contrastColors(m,:), 'LineStyle', 'none');
    alpha(0.2);
end