function violinPlotter(xCategories, yValues,randRange,bandwidth,width,color,plotV)

for ix = 1:length(xCategories)
    xRand = randi([-10 10],1,length(yValues{ix}))*(randRange*.1);
    plotScatter = scatter(xCategories(ix)+xRand,yValues{ix},10,'go');
    set(plotScatter,'MarkerFaceColor',color,'MarkerEdgeColor','none','MarkerFaceAlpha',.3)
    hold on
    
    meanY = nanmean(yValues{ix});
    semY = nanstd(yValues{ix})/sqrt(length(yValues{ix}));
    
    plotMean = errorbar(xCategories(ix),meanY,semY,'o',...
        'MarkerFaceColor',color,...
        'LineStyle','none',...
        'color',color,...
        'MarkerEdgeColor','none',...
        'capsize',0);
    
     if ~isempty(yValues{ix}) & strcmp(plotV,'on')
        [density,value] = ksdensity(yValues{ix},'Bandwidth', bandwidth);
        w = width*length(yValues{ix});
        pViolin = fill([density*w+xCategories(ix) -density(end:-1:1)*w+xCategories(ix)], [value value(end:-1:1)], color,...
            'LineStyle', '-',...
            'EdgeColor',color,...
            'EdgeAlpha',0.2,...
            'FaceAlpha',0.2);
        uistack(pViolin,'bottom');
    else
    end
end