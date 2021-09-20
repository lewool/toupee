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
subpop = 'all';

for iX = 1:length(expInfo)
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

    
    for iCell = 1:length(plotCells)
        for t = 5:36
            [cp(iCell,t), ~, ifSig(iCell,t)] = getCP(plotCells(iCell), squeeze(alignedResps{2}(:,t,:)), trialTypes.intVar.all.contrast_direction);
        end
%         [cp(iCell,t), ~, ifSig(iCell,t)] = getCP(plotCells(iCell), squeeze(nanmean(alignedResps{1}(:,22:26,:),2)), trialTypes.intVar.all.contrast_direction);
%         disp(strcat('session',num2str(iX),', t = ',num2str(eventWindow(t)),'s'))
        disp(strcat({'session '},num2str(iX),{', '},num2str(100*iCell/length(plotCells)),'% complete'))
    end
%     
%     figure;
% hold on;
% set(gcf,'position',[465   252   617   407]);
% line([0 0],[0 1],'LineStyle',':','Color',[.5 .5 .5])
% line([0.8 0.8],[0 1],'LineStyle',':','Color',[.5 .5 .5])
% line([-.5 1.5],[0.05 0.05],'LineStyle',':','Color',[.5 .5 .5])
% xlim([-.5 1.5])
% ylim([0 .8]);
% yticks([0 .2 .4 .6 .8 1])
% set(gca,'tickdir','out')
% xlabel('Time – stimOn (s)')
% ylabel('Proportion of cells with significant CP')

pcp = sum(ifSig)/length(plotCells);
pcp_all(iX,:) = pcp;
% plot(eventWindow,pcp,'LineWidth',2)
% title(strcat({expInfo(iX).expDate},{' '},{expInfo(iX).mouseName},{': '},{subpop},{' cells ('}, num2str(pct),{'%)'}))
end

% pcp = sum(ifSig)/length(plotCells);
% plot(eventWindow,pcp,'LineWidth',2)
% hold on;

% title(strcat({expInfo.expDate},{' '},{expInfo.mouseName},{': '},{subpop},{' cells ('}, num2str(pct),{'%)'}))