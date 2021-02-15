figure;
hold on;
set(gcf,'position',[465   252   617   407]);
line([0 0],[0 1],'LineStyle',':','Color',[.5 .5 .5])
line([0.8 0.8],[0 1],'LineStyle',':','Color',[.5 .5 .5])
line([-.5 1.5],[0.05 0.05],'LineStyle',':','Color',[.5 .5 .5])
xlim([-.5 1.5])
ylim([0 .8]);
yticks([0 .2 .4 .6 .8 1])
set(gca,'tickdir','out')
xlabel('Time – stimOn (s)')
ylabel('Proportion of cells with significant CP')
%%

subpop = 'contraMov';

for iX = 1:2
    trialTypes = getTrialTypes(expInfo(iX), behavioralData(iX), 'late');
    alignedResps = neuralData(iX).eta.alignedResps;
    eventWindow = neuralData(iX).eta.eventWindow;
    if strcmp(subpop,'contraMov') && expInfo(iX).hemisphere > 0
        plotCells = getWhichCells('leftMov',neuralData(iX));
    elseif strcmp(subpop,'contraMov') && expInfo(iX).hemisphere < 0
        plotCells = getWhichCells('rightMov',neuralData(iX));
    elseif strcmp(subpop,'contraStim') && expInfo(iX).hemisphere > 0
        plotCells = getWhichCells('leftStim',neuralData(iX));
    elseif strcmp(subpop,'contraStim') && expInfo(iX).hemisphere < 0
        plotCells = getWhichCells('rightStim',neuralData(iX));
    elseif strcmp(subpop,'ipsiMov') && expInfo(iX).hemisphere > 0
        plotCells = getWhichCells('rightMov',neuralData(iX));
    elseif strcmp(subpop,'ipsiMov') && expInfo(iX).hemisphere < 0
        plotCells = getWhichCells('leftMov',neuralData(iX));
    elseif strcmp(subpop,'ipsiStim') && expInfo(iX).hemisphere > 0
        plotCells = getWhichCells('rightStim',neuralData(iX));
    elseif strcmp(subpop,'ipsiStim') && expInfo(iX).hemisphere < 0
        plotCells = getWhichCells('leftStim',neuralData(iX));
    elseif strcmp(subpop,'all')
        plotCells = getWhichCells('all',neuralData(iX));
    end
    pct = round((length(plotCells)/size(alignedResps{1},3))*100);

    cp = zeros(length(plotCells),length(eventWindow));
    ifSig = zeros(length(plotCells),length(eventWindow));

    for t = 16:36
        for iCell = 1:length(plotCells)
            [cp(iCell,t), ~, ifSig(iCell,t)] = getCP(plotCells(iCell), squeeze(alignedResps{1}(:,t,:)), trialTypes.intVar.all.contrast_direction);
        end
        disp(strcat('session',num2str(iX),', t = ',num2str(eventWindow(t)),'s'))
    end
end

pcp = sum(ifSig)/length(plotCells);
plot(eventWindow,pcp,'LineWidth',2)
hold on;

% title(strcat({expInfo.expDate},{' '},{expInfo.mouseName},{': '},{subpop},{' cells ('}, num2str(pct),{'%)'}))