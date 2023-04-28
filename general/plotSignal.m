function plotSignal(x, signal, CIU, CIL, colors,ls)

for m = 1:size(signal,1)
    signal(m,:) = smooth(signal(m,:),1);
    plotSignal = plot(x,(signal(m,:)),'linestyle',ls,'linewidth',1,'color',colors(m,:));
    hold on;
    plotConfidence = fill([x';flipud(x')],[(CIL(m,:))';flipud((CIU(m,:))')],colors(m,:), 'LineStyle', 'none');
    alpha(0.2);
end