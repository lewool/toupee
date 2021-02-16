for iX = 1:14
    %% discard trials with early movements
contrasts = getUniqueContrasts(expInfo(iX));
[~, whichTrials] = selectCondition(expInfo(iX), contrasts, behavioralData(iX), initTrialConditions('movementTime','late'));
trialTypes = getTrialTypes(expInfo(iX), behavioralData(iX), 'late');
% [~, prevRightTrials] = selectCondition(expInfo(iX), contrasts, behavioralData(iX), initTrialConditions('movementTime','late','pastMovementDir','ccw'));
% [~, prevLeftTrials] = selectCondition(expInfo(iX), contrasts(iX), behavioralData, initTrialConditions('movementTime','late','pastMovementDir','cw'));

%% compute action vs stimulus coding for each cell
[baselineResps, stimResps, pmovResps, movResps, rewResps] = getEpochResps(neuralData(iX).eta);

action{iX} = nanmean((movResps(whichTrials,:)-baselineResps(whichTrials,:)),1);
stim{iX} = nanmean((stimResps(whichTrials,:)-baselineResps(whichTrials,:)),1);
relMovResps = movResps - baselineResps;
relStimResps = stimResps - baselineResps;

    if expInfo(iX).hemisphere < 0
        
        dir{iX} = ... 
            nanmean(relMovResps(trialTypes.singleVar.direction{1},:)) - ...
            nanmean(relMovResps(trialTypes.singleVar.direction{2},:));

        side{iX} = ...
            nanmean(relStimResps(trialTypes.singleVar.side{1},:)) - ...
            nanmean(relStimResps(trialTypes.singleVar.side{3},:));
    else
        
        dir{iX} = ... 
            nanmean(relMovResps(trialTypes.singleVar.direction{2},:)) - ...
            nanmean(relMovResps(trialTypes.singleVar.direction{1},:));

        side{iX} = ...
            nanmean(relStimResps(trialTypes.singleVar.side{3},:)) - ...
            nanmean(relStimResps(trialTypes.singleVar.side{1},:));
    end
end
%%
figure(1);
actionIndex = cat(2,action{:});
stimIndex = cat(2,stim{:});
minLim = min([actionIndex, stimIndex]);
maxLim = max([actionIndex, stimIndex]);
line([0 0],[minLim maxLim],'LineStyle',':','Color',[.5 .5 .5])
hold on;
line([minLim maxLim],[0 0],'LineStyle',':','Color',[.5 .5 .5])
line([minLim maxLim],[minLim maxLim],'LineStyle',':','Color',[.5 .5 .5])
scatter(stimIndex,actionIndex,2,...
    'MarkerEdgeColor','none','MarkerFaceColor','k','MarkerFaceAlpha', 0.7)
% scatplot(stimIndex,actionIndex,[],[],[],[],[],3);
xlim([minLim maxLim]);
ylim([minLim maxLim]);
axis square
box off
set(gca,'tickdir','out')
xlabel('stimulus response')
ylabel('movement response')

%%

dirIndex = cat(2,dir{:});
sideIndex = cat(2,side{:});

figure(2);
minLim = min([dirIndex, sideIndex]);
maxLim = max([dirIndex, sideIndex]);
line([0 0],[minLim maxLim],'LineStyle',':','Color',[.5 .5 .5])
hold on;
line([minLim maxLim],[0 0],'LineStyle',':','Color',[.5 .5 .5])
line([minLim maxLim],[minLim maxLim],'LineStyle',':','Color',[.5 .5 .5])
scatter(sideIndex,dirIndex,2,...
    'MarkerEdgeColor','none','MarkerFaceColor','k','MarkerFaceAlpha', .7);
% scatplot(sideIndex,dirIndex,[],[],[],[],[],3);
xlim([minLim maxLim]);
ylim([minLim maxLim]);
axis square
box off
set(gca,'tickdir','out')
xlabel('stimulus side preference')
ylabel('choice direction preference')

%% 

prevDirIndex = nanmean(baselineResps(prevRightTrials,:)) - nanmean(baselineResps(prevLeftTrials,:));
baselineStimIndex = nanmean(baselineResps(trialTypes.singleVar.side{3}(1:end-1),:)) - nanmean(baselineResps(trialTypes.singleVar.side{1}(1:end-1),:));

figure(1);
minLim = min([prevDirIndex, dirIndex]);
maxLim = max([prevDirIndex, dirIndex]);
line([0 0],[minLim maxLim],'LineStyle',':','Color',[.5 .5 .5])
hold on;
line([minLim maxLim],[0 0],'LineStyle',':','Color',[.5 .5 .5])
line([minLim maxLim],[minLim maxLim],'LineStyle',':','Color',[.5 .5 .5])
scatter(dirIndex,prevDirIndex,'k.')
xlim([minLim maxLim]);
ylim([minLim maxLim]);
axis square
box off
set(gca,'tickdir','out')
xlabel('choice direction preference (t)')
ylabel('choice direction preference (t–1)')


%%

whichNeuron = 623;
whichNeuron = neuralData(1).stats.bfcH(:,2) == 1 ;
% whichNeuron = 1:length(neuralData(1).stats.bfcH);
clear sLmL sLmR
for c = 1:length(trialTypes.intVar.cb2D.contrast_direction)
sLmL(c,:) = nanmean(stimResps(trialTypes.intVar.cb2D.contrast_direction{c,1},whichNeuron)-baselineResps(trialTypes.intVar.cb2D.contrast_direction{c,1},whichNeuron),1);
sLmR(c,:) = nanmean(stimResps(trialTypes.intVar.cb2D.contrast_direction{c,2},whichNeuron)-baselineResps(trialTypes.intVar.cb2D.contrast_direction{c,2},whichNeuron),1);
end



figure;
for c = 1:9
    subplot(1,9,c)
    scatter(sLmL(c,:),sLmR(c,:),'k.')
    hold on
    line([-.1 .25],[-.1 .25],'LineStyle',':','Color',[.5 .5 .5])
    xlim([-.05 .25]);
    ylim([-.05 .25]);
    axis square
    box off
end