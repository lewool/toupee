oriCond = [];
uniqueOris = unique(expInfo(2).block.events.oriValues);
for o = 1:length(uniqueOris)
    oriCond{o} = find(expInfo(2).block.events.oriValues(1:length(expInfo(2).block.events.endTrialValues)) == uniqueOris(o));
end

[baselineResps, stimResps] = getStimResps(neuralData(2).eta);
iCell = 2491;

for o = 1:length(uniqueOris)
    OSI(o) = mean(stimResps(oriCond{o},iCell) - neuralData(2).eta.alignedResps{1}(oriCond{o},20,iCell));
end

figure;
plot(uniqueOris, OSI,'ko','LineStyle','-','MarkerFaceColor','k')
set(gca,'tickdir','out')
xlim([-10 370])
xticks([0 45 90 135 180 225 270 315 360])
% set(gca, 'XTickLabels', {'0', '45', '90', '135', '180', '225', '270', '315', '360'})
box off
xlabel('grating direction')
ylabel('response (stim – baseline)')

figure;
set(gcf,'position', [280 1000 1540 130]);
for o = 1:length(uniqueOris)
    subplot(1,9,o)
    maxy(o) = max((nanmean(neuralData(2).eta.alignedResps{1}(oriCond{o},:,iCell),1)  - mean(baselineResps(oriCond{o},iCell))));
    plot(eventWindow,(nanmean(neuralData(2).eta.alignedResps{1}(oriCond{o},:,iCell),1) - mean(neuralData(2).eta.alignedResps{1}(oriCond{o},20,iCell))),'LineWidth',1,'Color','k');
end

for o = 1:length(uniqueOris)
    subplot(1,9,o)
    line([0 0],[-max(maxy) max(maxy)],'LineStyle','--','Color',[.5 .5 .5])
    ylim([-max(maxy) max(maxy)]);
    xlim([-0.5 2]);
    box off
    set(gca,'tickdir','out')
    title(num2str(uniqueOris(o)))
end