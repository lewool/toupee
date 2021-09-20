for iX = 1:5
    %% discard trials with early movements
contrasts = getUniqueContrasts(expInfo(iX));
[~, whichTrials] = selectCondition(expInfo(iX), contrasts, behavioralData(iX), initTrialConditions('movementTime','late'));
trialTypes = getTrialTypes(expInfo(iX), behavioralData(iX), 'late');
% [~, prevRightTrials] = selectCondition(expInfo(iX), contrasts, behavioralData(iX), initTrialConditions('movementTime','late','pastMovementDir','ccw'));
% [~, prevLeftTrials] = selectCondition(expInfo(iX), contrasts(iX), behavioralData, initTrialConditions('movementTime','late','pastMovementDir','cw'));

[~, earlyTrials] = selectCondition(expInfo(iX), getUniqueContrasts(expInfo(iX)), behavioralData(iX), ...
    initTrialConditions('preStimMovement','quiescent','movementTime','early','specificRTs',[.1 inf]));

    [~, lateTrials] = selectCondition(expInfo(iX), getUniqueContrasts(expInfo(iX)), behavioralData(iX), ...
    initTrialConditions('preStimMovement','quiescent','movementTime','late','specificRTs',[.1 inf]));

[correct_highSide, correct_lowSide, incorrect] = getHighLoRewardTrials(expInfo(iX), behavioralData(iX));

%% compute action vs stimulus coding for each cell
[baselineResps, stimResps, pmovResps, movResps, rewResps] = getEpochResps(neuralData(iX).eta);

action{iX} = nanmean((movResps(whichTrials,:)-baselineResps(whichTrials,:)),1);
stim{iX} = nanmean((stimResps(whichTrials,:)-baselineResps(whichTrials,:)),1);
rew{iX} = nanmean((rewResps(whichTrials,:)-baselineResps(whichTrials,:)),1);
relMovResps = movResps - baselineResps;
relStimResps = stimResps - baselineResps;
relRewResps = rewResps - baselineResps;

    if expInfo(iX).hemisphere < 0
        
        dir{iX} = ... 
            nanmean(relMovResps(trialTypes.singleVar.direction{1},:)) - ...
            nanmean(relMovResps(trialTypes.singleVar.direction{2},:));

        side{iX} = ...
            nanmean(relStimResps(trialTypes.singleVar.side{1},:)) - ...
            nanmean(relStimResps(trialTypes.singleVar.side{3},:));
        
        feed{iX} = ...
            nanmean(relRewResps(trialTypes.singleVar.outcome{1},:)) - ...
            nanmean(relRewResps(trialTypes.singleVar.outcome{2},:));
        
        timing{iX} = ...
            nanmean(relMovResps(earlyTrials,:)) - ...
            nanmean(relMovResps(lateTrials,:));
        
        hilo{iX} = ...
            nanmean(relRewResps(correct_highSide,:)) - ...
            nanmean(relRewResps(correct_lowSide,:));
    else
        
        dir{iX} = ... 
            nanmean(relMovResps(trialTypes.singleVar.direction{2},:)) - ...
            nanmean(relMovResps(trialTypes.singleVar.direction{1},:));

        side{iX} = ...
            nanmean(relStimResps(trialTypes.singleVar.side{3},:)) - ...
            nanmean(relStimResps(trialTypes.singleVar.side{1},:));
        
        feed{iX} = ...
            nanmean(relRewResps(trialTypes.singleVar.outcome{1},:)) - ...
            nanmean(relRewResps(trialTypes.singleVar.outcome{2},:));
        
        timing{iX} = ...
            nanmean(relMovResps(earlyTrials,:)) - ...
            nanmean(relMovResps(lateTrials,:));
        
        hilo{iX} = ...
            nanmean(relRewResps(correct_highSide,:)) - ...
            nanmean(relRewResps(correct_lowSide,:));
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
% scatter(stimIndex,actionIndex,2,...
%     'MarkerEdgeColor','none','MarkerFaceColor','k','MarkerFaceAlpha', 0.7)
scatplot(stimIndex,actionIndex,[],[],[],[],[],4);
xlim([minLim maxLim]);
ylim([minLim maxLim]);
axis square
box off
set(gca,'tickdir','out')
xlabel('stimulus epoch')
ylabel('movement epoch')

figure(10);
actionIndex = cat(2,action{:});
rewIndex = cat(2,rew{:});
minLim = min([actionIndex, rewIndex]);
maxLim = max([actionIndex, rewIndex]);
line([0 0],[minLim maxLim],'LineStyle',':','Color',[.5 .5 .5])
hold on;
line([minLim maxLim],[0 0],'LineStyle',':','Color',[.5 .5 .5])
line([minLim maxLim],[minLim maxLim],'LineStyle',':','Color',[.5 .5 .5])
% scatter(rewIndex,actionIndex,2,...
%     'MarkerEdgeColor','none','MarkerFaceColor','k','MarkerFaceAlpha', 0.7)
scatplot(rewIndex,actionIndex,[],[],[],[],[],4);
xlim([minLim maxLim]);
ylim([minLim maxLim]);
axis square
box off
set(gca,'tickdir','out')
xlabel('reward epoch')
ylabel('movement epoch')


%%

dirIndex = cat(2,dir{:});
sideIndex = cat(2,side{:});
rewIdx = cat(2,feed{:});
timeIdx = cat(2,timing{:});
hiloIdx = cat(2,hilo{:});

figure(2);
minLim = min([dirIndex, sideIndex]);
maxLim = max([dirIndex, sideIndex]);
line([0 0],[minLim maxLim],'LineStyle',':','Color',[.5 .5 .5])
hold on;
line([minLim maxLim],[0 0],'LineStyle',':','Color',[.5 .5 .5])
% line([minLim maxLim],[minLim maxLim],'LineStyle',':','Color',[.5 .5 .5])
% scatter(sideIndex,dirIndex,2,...
%     'MarkerEdgeColor','none','MarkerFaceColor','k','MarkerFaceAlpha', .7);
scatplot(sideIndex,dirIndex,[],[],[],[],[],4);
xlim([minLim maxLim]);
ylim([minLim maxLim]);
axis square
box off
set(gca,'tickdir','out')
xlabel('stimulus side selectivity')
ylabel('choice selectivity')

figure(20);
minLim = min([dirIndex, rewIdx]);
maxLim = max([dirIndex, rewIdx]);
line([0 0],[minLim maxLim],'LineStyle',':','Color',[.5 .5 .5])
hold on;
line([minLim maxLim],[0 0],'LineStyle',':','Color',[.5 .5 .5])
% line([minLim maxLim],[minLim maxLim],'LineStyle',':','Color',[.5 .5 .5])
% scatter(sideIndex,dirIndex,2,...
%     'MarkerEdgeColor','none','MarkerFaceColor','k','MarkerFaceAlpha', .7);
scatplot(rewIdx,dirIndex,[],[],[],[],[],4);
xlim([minLim maxLim]);
ylim([minLim maxLim]);
axis square
box off
set(gca,'tickdir','out')
xlabel('reward selectivity')
ylabel('choice selectivity')

figure(200);
minLim = min([dirIndex, timeIdx]);
maxLim = max([dirIndex, timeIdx]);
line([0 0],[minLim maxLim],'LineStyle',':','Color',[.5 .5 .5])
hold on;
line([minLim maxLim],[0 0],'LineStyle',':','Color',[.5 .5 .5])
% line([minLim maxLim],[minLim maxLim],'LineStyle',':','Color',[.5 .5 .5])
% scatter(sideIndex,dirIndex,2,...
%     'MarkerEdgeColor','none','MarkerFaceColor','k','MarkerFaceAlpha', .7);
scatplot(timeIdx,dirIndex,[],[],[],[],[],4);
xlim([minLim maxLim]);
ylim([minLim maxLim]);
axis square
box off
set(gca,'tickdir','out')
xlabel('timing selectivity')
ylabel('choice selectivity')

figure(2000);
minLim = min([rewIdx, hiloIdx]);
maxLim = max([rewIdx, hiloIdx]);
line([0 0],[minLim maxLim],'LineStyle',':','Color',[.5 .5 .5])
hold on;
line([minLim maxLim],[0 0],'LineStyle',':','Color',[.5 .5 .5])
% line([minLim maxLim],[minLim maxLim],'LineStyle',':','Color',[.5 .5 .5])
% scatter(sideIndex,dirIndex,2,...
%     'MarkerEdgeColor','none','MarkerFaceColor','k','MarkerFaceAlpha', .7);
scatplot(dirIndex,rewIdx,[],[],[],[],[],4);
xlim([minLim maxLim]);
ylim([minLim maxLim]);
axis square
box off
set(gca,'tickdir','out')
xlabel('reward selectivity (yes/no)')
ylabel('reward selectivity (hi/lo)')


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

whichNeuron = 16;
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