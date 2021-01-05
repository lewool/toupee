expInfo = initExpInfo({{'LEW031'}},{{'2020-02-03',1,[1]}});
[expInfo, neuralData, behavioralData] = processExperiment(expInfo);
trialTypes = getTrialTypes(expInfo, behavioralData, 'late');
alignedResps = neuralData.eta.alignedResps;
eventWindow = neuralData.eta.eventWindow;


%%
subpop = 'advanceMov';
plotCells = getWhichCells(subpop,neuralData);
pct = round((length(plotCells)/size(alignedResps{1},3))*100);

cp = zeros(length(plotCells),length(eventWindow));
ifSig = zeros(length(plotCells),length(eventWindow));

for t = 1:length(eventWindow)
    for iCell = 1:length(plotCells)
        [cp(iCell,t), ~, ifSig(iCell,t)] = getCP(plotCells(iCell), squeeze(alignedResps{1}(:,t,:)), trialTypes.intVar.all.contrast_direction);
    end
    disp(strcat(num2str(round((t/length(eventWindow))*100)),'% complete'))
end

pcp = sum(ifSig)/length(plotCells);
figure;
hold on;
set(gcf,'position',[465   252   617   407]);
plot(eventWindow,pcp,'LineWidth',2)
line([0 0],[0 .6],'LineStyle',':','Color',[.5 .5 .5])
line([0.8 0.8],[0 .6],'LineStyle',':','Color',[.5 .5 .5])
xlim([-1 2])
xlabel('Time from stimulus onset (s)')
ylabel('Proportion of cells with significant CP')
title(strcat({expInfo.expDate},{' '},{expInfo.mouseName},{': '},{subpop},{' cells ('}, num2str(pct),{'%)'}))