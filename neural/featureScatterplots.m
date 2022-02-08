function [fit, correlation] = featureScatterplots(expInfo, behavioralData,neuralData, hemisphere)

expInfo.hemisphere = hemisphere;

    %% discard trials with early movements
contrasts = getUniqueContrasts(expInfo);
[~, whichTrials] = selectCondition(expInfo, contrasts, behavioralData, initTrialConditions('movementTime','late'));
trialTypes = getTrialTypes(expInfo, behavioralData, 'late');
% [~, prevRightTrials] = selectCondition(expInfo(iX), contrasts, behavioralData(iX), initTrialConditions('movementTime','late','pastMovementDir','ccw'));
% [~, prevLeftTrials] = selectCondition(expInfo(iX), contrasts(iX), behavioralData, initTrialConditions('movementTime','late','pastMovementDir','cw'));

[~, earlyTrials] = selectCondition(expInfo, getUniqueContrasts(expInfo), behavioralData, ...
    initTrialConditions('preStimMovement','quiescent','movementTime','early','specificRTs',[.1 inf]));

    [~, lateTrials] = selectCondition(expInfo, getUniqueContrasts(expInfo), behavioralData, ...
    initTrialConditions('preStimMovement','quiescent','movementTime','late','specificRTs',[.1 inf]));

[correct_highSide, correct_lowSide, incorrect] = getHighLoRewardTrials(expInfo, behavioralData);

% whichCells = find(neuralData.stats.bfcH(:,5) | neuralData.stats.bfcH(:,6));
% whichCells = find(neuralData.stats.bfcH(:,5));
whichCells = 1:length(neuralData.stats.bfcH(:,5));


%% compute action vs stimulus coding for each cell
[baselineResps, stimResps, pmovResps, movResps, rewResps] = getEpochResps(neuralData.eta);

action = nanmean((movResps(whichTrials,:)));%-baselineResps(whichTrials,:)),1);
stim = nanmean((stimResps(whichTrials,:)));%-baselineResps(whichTrials,:)),1);
rew = nanmean((rewResps(whichTrials,:)));%-baselineResps(whichTrials,:)),1);
relMovResps = movResps;
relStimResps = stimResps;
relRewResps = rewResps;
        
if expInfo.hemisphere < 0

    dir = ... 
        (nanmean(relMovResps(trialTypes.intVar.all.contrast_direction{1},whichCells)) - ...
        nanmean(relMovResps(trialTypes.intVar.all.contrast_direction{2},whichCells))) ./ ...
        (nanmean(relMovResps(trialTypes.intVar.all.contrast_direction{1},whichCells)) + ...
        nanmean(relMovResps(trialTypes.intVar.all.contrast_direction{2},whichCells)));

    side = ...
        (nanmean(relStimResps(trialTypes.singleVar.side{1},whichCells)) - ...
        nanmean(relStimResps(trialTypes.singleVar.side{3},whichCells))) ./ ...
        (nanmean(relStimResps(trialTypes.singleVar.side{1},whichCells)) + ...
        nanmean(relStimResps(trialTypes.singleVar.side{3},whichCells)));

    feed = ...
        (nanmean(relRewResps(trialTypes.singleVar.outcome{1},whichCells)) - ...
        nanmean(relRewResps(trialTypes.singleVar.outcome{2},whichCells))) ./ ...
        (nanmean(relRewResps(trialTypes.singleVar.outcome{1},whichCells)) - ...
        nanmean(relRewResps(trialTypes.singleVar.outcome{2},whichCells)));

    timing = ...
        (nanmean(relMovResps(earlyTrials,whichCells)) - ...
        nanmean(relMovResps(lateTrials,whichCells)));

    hilo = ...
        (nanmean(relRewResps(correct_highSide,whichCells)) - ...
        nanmean(relRewResps(correct_lowSide,whichCells)));
else

    dir = ... 
        (nanmean(relMovResps(trialTypes.intVar.all.contrast_direction{2},whichCells)) - ...
        nanmean(relMovResps(trialTypes.intVar.all.contrast_direction{1},whichCells))) ./ ...
        (nanmean(relMovResps(trialTypes.intVar.all.contrast_direction{2},whichCells)) + ...
        nanmean(relMovResps(trialTypes.intVar.all.contrast_direction{1},whichCells)));

    side = ...
        (nanmean(relStimResps(trialTypes.singleVar.side{3},whichCells)) - ...
        nanmean(relStimResps(trialTypes.singleVar.side{1},whichCells))) ./ ...
        (nanmean(relStimResps(trialTypes.singleVar.side{3},whichCells)) + ...
        nanmean(relStimResps(trialTypes.singleVar.side{1},whichCells)));

    feed = ...
        (nanmean(relRewResps(trialTypes.singleVar.outcome{1},whichCells)) - ...
        nanmean(relRewResps(trialTypes.singleVar.outcome{2},whichCells))) ./ ...
        (nanmean(relRewResps(trialTypes.singleVar.outcome{1},whichCells)) + ...
        nanmean(relRewResps(trialTypes.singleVar.outcome{2},whichCells)));

    timing = ...
        (nanmean(relMovResps(earlyTrials,whichCells)) - ...
        nanmean(relMovResps(lateTrials,whichCells)));

    hilo = ...
        (nanmean(relRewResps(correct_highSide,whichCells)) - ...
        nanmean(relRewResps(correct_lowSide,whichCells)));
end

%%
meanColor = [1 0 0];
ms = 4;
mm = 50;
% close all;
dirIndex = dir;
sideIndex = side;
rewIdx = feed;
timeIdx = timing;
hiloIdx = hilo;
[~,sortDirIdx] = sort(dirIndex);
[~,sortRewIdx] = sort(rewIdx);

fig = figure(1);
set(gcf,'position', [1000 1050 948 290])
subplot(1,3,1)
minLim = min([dirIndex, sideIndex]);
maxLim = max([dirIndex, sideIndex]);
% line([0 0],[minLim maxLim],'LineStyle',':','Color',[.5 .5 .5])
hold on;
% line([minLim maxLim],[0 0],'LineStyle',':','Color',[.5 .5 .5])
% line([minLim maxLim],[minLim maxLim],'LineStyle',':','Color',[.5 .5 .5])
scatter(dirIndex,sideIndex,ms,...
    'MarkerEdgeColor','none','MarkerFaceColor','k','MarkerFaceAlpha', .7);
% scatplot(sideIndex,dirIndex,[],[],[],[],[],4);
% plot(movmean(dirIndex(sortDirIdx),1),movmean(sideIndex(sortDirIdx),mm),'Color',meanColor,'LineWidth',2);
xlim([minLim maxLim]);
ylim([minLim maxLim]);
axis square
box off
set(gca,'tickdir','out')
ylabel('stim selectivity (ipsi - contra)')
xlabel('choice selectivity (ipsi - contra)')

subplot(1,3,2)
minLim = min([dirIndex, rewIdx]);
maxLim = max([dirIndex, rewIdx]);
% line([0 0],[minLim maxLim],'LineStyle',':','Color',[.5 .5 .5])
hold on;
% line([minLim maxLim],[0 0],'LineStyle',':','Color',[.5 .5 .5])
% line([minLim maxLim],[minLim maxLim],'LineStyle',':','Color',[.5 .5 .5])
scatter(dirIndex,rewIdx,ms,...
    'MarkerEdgeColor','none','MarkerFaceColor','k','MarkerFaceAlpha', .7);
% scatplot(rewIdx,dirIndex,[],[],[],[],[],4);
% plot(movmean(dirIndex(sortDirIdx),1),movmean(rewIdx(sortDirIdx),mm),'Color',meanColor,'LineWidth',2);
xlim([minLim maxLim]);
ylim([minLim maxLim]);
axis square
box off
set(gca,'tickdir','out')
ylabel('reward selectivity (correct - error)')
xlabel('choice selectivity (ipsi - contra)')

% figure(200);
% minLim = min([dirIndex, timeIdx]);
% maxLim = max([dirIndex, timeIdx]);
% line([0 0],[minLim maxLim],'LineStyle',':','Color',[.5 .5 .5])
% hold on;
% line([minLim maxLim],[0 0],'LineStyle',':','Color',[.5 .5 .5])
% % line([minLim maxLim],[minLim maxLim],'LineStyle',':','Color',[.5 .5 .5])
% scatter(dirIndex,timeIdx,ms,...
%     'MarkerEdgeColor','none','MarkerFaceColor','k','MarkerFaceAlpha', .7);
% % scatplot(timeIdx,dirIndex,[],[],[],[],[],4);
% plot(movmean(dirIndex(sortDirIdx),1),movmean(timeIdx(sortDirIdx),mm),'Color',meanColor,'LineWidth',2);
% xlim([minLim maxLim]);
% ylim([minLim maxLim]);
% axis square
% box off
% set(gca,'tickdir','out')
% ylabel('timing selectivity')
% xlabel('choice selectivity')

subplot(1,3,3)
minLim = min([sideIndex, rewIdx]);
maxLim = max([sideIndex, rewIdx]);
% line([0 0],[minLim maxLim],'LineStyle',':','Color',[.5 .5 .5])
hold on;
% line([minLim maxLim],[0 0],'LineStyle',':','Color',[.5 .5 .5])
% line([minLim maxLim],[minLim maxLim],'LineStyle',':','Color',[.5 .5 .5])
scatter(sideIndex,rewIdx,ms,...
    'MarkerEdgeColor','none','MarkerFaceColor','k','MarkerFaceAlpha', .7);
% scatplot(hiloIdx,dirIndex,[],[],[],[],[],4);
% plot(movmean(dirIndex(sortDirIdx),1),movmean(hiloIdx(sortDirIdx),mm),'Color',meanColor,'LineWidth',2);
xlim([minLim maxLim]);
ylim([minLim maxLim]);
axis square
box off
set(gca,'tickdir','out')
ylabel('reward selectivity (correct - error)')
xlabel('side selectivity (ipsi - contra)')

fit.dir_side = polyfit(dirIndex,sideIndex,1);
fit.dir_reward = polyfit(dirIndex,rewIdx,1);
fit.side_reward = polyfit(sideIndex,rewIdx,1);

correlation.dir_side = corr(dirIndex',sideIndex');
correlation.dir_reward = corr(dirIndex',rewIdx');
correlation.side_reward = corr(sideIndex',rewIdx');


%%

[expRef, ~] = data.constructExpRef(expInfo.mouseName,expInfo.expDate,expInfo.expNum);

printfig(fig, strcat('feature scatterplots zs_',expRef))