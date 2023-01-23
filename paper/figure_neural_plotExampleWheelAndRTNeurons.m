vels = abs(behavioralData.wheelMoves.epochs(5).peakVel);
[baselineResps, stimResps, pmovResps, movResps, rewResps, preCueResps] = getEpochResps(neuralData.eta);
RTs = behavioralData.wheelMoves.epochs(5).onsetTimes - behavioralData.eventTimes(1).daqTime;
trials = find((RTs < prctile(RTs,97.5)) & (RTs > prctile(RTs,2.5)));

%%

v{1} = find(vels < prctile(vels,20));
v{2} = find((vels >= prctile(vels,20)) & (vels < prctile(vels,40)));
v{3} = find((vels >= prctile(vels,40)) & (vels < prctile(vels,60)));
v{4} = find((vels >= prctile(vels,60)) & (vels < prctile(vels,80)));
v{5} = find(vels >= prctile(vels,80));

figure;
set(gcf,'position',[1250 1220 1220 406])

cellno = 27;
lineColors = colormapThruGray([.75 0 .75],[0 .75 0],6,.5);
subplot(1,2,1)
for pv = 1:5
    plot(neuralData.eta.eventWindow,nanmean(neuralData.eta.alignedResps{2}(v{pv},:,cellno),1),'color',lineColors(pv,:),'linewidth',1)
    hold on
end
line([0 0],[0 .6],'Color',[.5 .5 .5],'LineStyle','--')

prettyPlot(gca)
legend({'20th %ile' '40th' '60th' '80th' '100th'})
legend('boxoff')
xlabel('Time from movement (s)')
ylabel('Response')

subplot(1,2,2)
hold on
for pv = 1:length(v)
    scatter(vels(v{pv}), movResps(v{pv},cellno),12,'MarkerFaceColor',lineColors(pv,:),'MarkerEdgeColor','none')
end
prettyPlot(gca)
xlabel('Peak wheel velocity (mm/s)')
ylabel('Mean movement activity')

%%
    
r{1} = find((RTs >= prctile(RTs(trials),0)) & (RTs < prctile(RTs(trials),20)));
r{2} = find((RTs >= prctile(RTs(trials),20)) & (RTs < prctile(RTs(trials),40)));
r{3} = find((RTs >= prctile(RTs(trials),40)) & (RTs < prctile(RTs(trials),60)));
r{4} = find((RTs >= prctile(RTs(trials),60)) & (RTs < prctile(RTs(trials),80)));
r{5} = find((RTs >= prctile(RTs(trials),80)) & (RTs < prctile(RTs(trials),100)));

figure;
set(gcf,'position',[1250 1220 1220 406])
cellno = 894;
lineColors = colormapThruGray([1 0 .75],[0 .75 1],6,.5);

subplot(1,2,1)
for pv = 1:5
    plot(neuralData.eta.eventWindow,(nanmean(neuralData.eta.alignedResps{1}(r{pv},:,cellno),1)),'color',lineColors(pv,:),'linewidth',1)
    hold on
end
line([0 0],[0 .3],'Color',[.5 .5 .5],'LineStyle','--')
prettyPlot(gca)
legend({'20th %ile' '40th' '60th' '80th' '100th'})
legend('boxoff')
xlabel('Time from stim on (s)')
ylabel('Response')

subplot(1,2,2)
hold on
for pr = 1:length(r)
    scatter(RTs(r{pr}), baselineResps(r{pr},cellno),12,'MarkerFaceColor',lineColors(pr,:),'MarkerEdgeColor','none')
end
prettyPlot(gca)
xlabel('Response times (s)')
ylabel('Mean baseline activity')

figure;
plot(xcorr(baselineResps(:,cellno),20,'normalized'));
xlim([1 41])
xticks([1 11 21 31 41])
set(gca, 'XTickLabels', {'-20' '-10' '0' '+10' '+20'})
xlabel('Lags')
ylabel('Correlation')
prettyPlot(gca)
