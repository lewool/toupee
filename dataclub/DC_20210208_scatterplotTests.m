clear action stim rew dir side feed timing hilo
for iX = 1
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

whichCells = find(neuralData.stats.bfcH(:,5) | neuralData.stats.bfcH(:,6));
%% compute action vs stimulus coding for each cell
[baselineResps, stimResps, pmovResps, movResps, rewResps] = getEpochResps(neuralData(iX).eta);

action{iX} = nanmean((movResps(whichTrials,whichCells)));%-baselineResps(whichTrials,:)),1);
stim{iX} = nanmean((stimResps(whichTrials,whichCells)));%-baselineResps(whichTrials,:)),1);
rew{iX} = nanmean((rewResps(whichTrials,whichCells)));%-baselineResps(whichTrials,:)),1);
relMovResps = movResps;
relStimResps = stimResps;
relRewResps = rewResps;
        
    if expInfo(iX).hemisphere < 0
        
        dir{iX} = ... 
            (nanmean(relMovResps(trialTypes.singleVar.direction{1},whichCells)) - ...
            nanmean(relMovResps(trialTypes.singleVar.direction{2},whichCells)));

        side{iX} = ...
            (nanmean(relStimResps(trialTypes.singleVar.side{1},whichCells)) - ...
            nanmean(relStimResps(trialTypes.singleVar.side{3},whichCells)));
        
        feed{iX} = ...
            (nanmean(relRewResps(trialTypes.singleVar.outcome{1},whichCells)) - ...
            nanmean(relRewResps(trialTypes.singleVar.outcome{2},whichCells)));
        
        timing{iX} = ...
            (nanmean(relMovResps(earlyTrials,whichCells)) - ...
            nanmean(relMovResps(lateTrials,whichCells)));
        
        hilo{iX} = ...
            (nanmean(relRewResps(correct_highSide,whichCells)) - ...
            nanmean(relRewResps(correct_lowSide,whichCells)));
    else
        
        dir{iX} = ... 
            (nanmean(relMovResps(trialTypes.singleVar.direction{2},whichCells)) - ...
            nanmean(relMovResps(trialTypes.singleVar.direction{1},whichCells)));

        side{iX} = ...
            (nanmean(relStimResps(trialTypes.singleVar.side{3},whichCells)) - ...
            nanmean(relStimResps(trialTypes.singleVar.side{1},whichCells)));
        
        feed{iX} = ...
            (nanmean(relRewResps(trialTypes.singleVar.outcome{1},whichCells)) - ...
            nanmean(relRewResps(trialTypes.singleVar.outcome{2},whichCells)));
        
        timing{iX} = ...
            (nanmean(relMovResps(earlyTrials,whichCells)) - ...
            nanmean(relMovResps(lateTrials,whichCells)));
        
        hilo{iX} = ...
            (nanmean(relRewResps(correct_highSide,whichCells)) - ...
            nanmean(relRewResps(correct_lowSide,whichCells)));
    end
end
%%
% figure(1);
% actionIndex = cat(2,action{:});
% stimIndex = cat(2,stim{:});
% minLim = min([actionIndex, stimIndex]);
% maxLim = max([actionIndex, stimIndex]);
% line([0 0],[minLim maxLim],'LineStyle',':','Color',[.5 .5 .5])
% hold on;
% line([minLim maxLim],[0 0],'LineStyle',':','Color',[.5 .5 .5])
% line([minLim maxLim],[minLim maxLim],'LineStyle',':','Color',[.5 .5 .5])
% % scatter(stimIndex,actionIndex,2,...
% %     'MarkerEdgeColor','none','MarkerFaceColor','k','MarkerFaceAlpha', 0.7)
% scatplot(stimIndex,actionIndex,[],[],[],[],[],4);
% xlim([minLim maxLim]);
% ylim([minLim maxLim]);
% axis square
% box off
% set(gca,'tickdir','out')
% xlabel('stimulus epoch')
% ylabel('movement epoch')
% 
% figure(10);
% actionIndex = cat(2,action{:});
% rewIndex = cat(2,rew{:});
% minLim = min([actionIndex, rewIndex]);
% maxLim = max([actionIndex, rewIndex]);
% line([0 0],[minLim maxLim],'LineStyle',':','Color',[.5 .5 .5])
% hold on;
% line([minLim maxLim],[0 0],'LineStyle',':','Color',[.5 .5 .5])
% line([minLim maxLim],[minLim maxLim],'LineStyle',':','Color',[.5 .5 .5])
% % scatter(rewIndex,actionIndex,2,...
% %     'MarkerEdgeColor','none','MarkerFaceColor','k','MarkerFaceAlpha', 0.7)
% scatplot(rewIndex,actionIndex,[],[],[],[],[],4);
% xlim([minLim maxLim]);
% ylim([minLim maxLim]);
% axis square
% box off
% set(gca,'tickdir','out')
% xlabel('reward epoch')
% ylabel('movement epoch')


%%
meanColor = [1 0 0];
ms = 4;
mm = 50;
% close all;
dirIndex = cat(2,dir{:});
sideIndex = cat(2,side{:});
rewIdx = cat(2,feed{:});
timeIdx = cat(2,timing{:});
hiloIdx = cat(2,hilo{:});
[~,sortDirIdx] = sort(dirIndex);
[~,sortRewIdx] = sort(rewIdx);

figure;
set(gcf,'position', [1000 1050 948 290])
subplot(1,3,1)
minLim = min([dirIndex, sideIndex]);
maxLim = max([dirIndex, sideIndex]);
line([0 0],[minLim maxLim],'LineStyle',':','Color',[.5 .5 .5])
hold on;
line([minLim maxLim],[0 0],'LineStyle',':','Color',[.5 .5 .5])
% line([minLim maxLim],[minLim maxLim],'LineStyle',':','Color',[.5 .5 .5])
scatter(dirIndex,sideIndex,ms,...
    'MarkerEdgeColor','none','MarkerFaceColor','k','MarkerFaceAlpha', .7);
% scatplot(sideIndex,dirIndex,[],[],[],[],[],4);
plot(movmean(dirIndex(sortDirIdx),1),movmean(sideIndex(sortDirIdx),mm),'Color',meanColor,'LineWidth',2);
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
line([0 0],[minLim maxLim],'LineStyle',':','Color',[.5 .5 .5])
hold on;
line([minLim maxLim],[0 0],'LineStyle',':','Color',[.5 .5 .5])
% line([minLim maxLim],[minLim maxLim],'LineStyle',':','Color',[.5 .5 .5])
scatter(dirIndex,rewIdx,ms,...
    'MarkerEdgeColor','none','MarkerFaceColor','k','MarkerFaceAlpha', .7);
% scatplot(rewIdx,dirIndex,[],[],[],[],[],4);
plot(movmean(dirIndex(sortDirIdx),1),movmean(rewIdx(sortDirIdx),mm),'Color',meanColor,'LineWidth',2);
xlim([minLim maxLim]);
ylim([minLim maxLim]);
axis square
box off
set(gca,'tickdir','out')
ylabel('reward selectivity (hit - miss)')
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
minLim = min([dirIndex, hiloIdx]);
maxLim = max([dirIndex, hiloIdx]);
line([0 0],[minLim maxLim],'LineStyle',':','Color',[.5 .5 .5])
hold on;
line([minLim maxLim],[0 0],'LineStyle',':','Color',[.5 .5 .5])
% line([minLim maxLim],[minLim maxLim],'LineStyle',':','Color',[.5 .5 .5])
scatter(dirIndex,hiloIdx,ms,...
    'MarkerEdgeColor','none','MarkerFaceColor','k','MarkerFaceAlpha', .7);
% scatplot(hiloIdx,dirIndex,[],[],[],[],[],4);
plot(movmean(dirIndex(sortDirIdx),1),movmean(hiloIdx(sortDirIdx),mm),'Color',meanColor,'LineWidth',2);
xlim([minLim maxLim]);
ylim([minLim maxLim]);
axis square
box off
set(gca,'tickdir','out')
ylabel('reward selectivity (hi - lo)')
xlabel('choice selectivity (ipsi - contra)')

% figure(2000);
% minLim = min([rewIdx, hiloIdx]);
% maxLim = max([rewIdx, hiloIdx]);
% line([0 0],[minLim maxLim],'LineStyle',':','Color',[.5 .5 .5])
% hold on;
% line([minLim maxLim],[0 0],'LineStyle',':','Color',[.5 .5 .5])
% % line([minLim maxLim],[minLim maxLim],'LineStyle',':','Color',[.5 .5 .5])
% scatter(hiloIdx,rewIdx,ms,...
%     'MarkerEdgeColor','none','MarkerFaceColor','k','MarkerFaceAlpha', .7);
% % scatplot(hiloIdx,rewIdx,[],[],[],[],[],4);
% plot(movmean(hiloIdx(sortRewIdx),mm),movmean(rewIdx(sortRewIdx),1),'Color',meanColor,'LineWidth',2)
% xlim([minLim maxLim]);
% ylim([minLim maxLim]);
% axis square
% box off
% set(gca,'tickdir','out')
% xlabel('reward selectivity (hi/lo)')
% ylabel('reward selectivity (yes/no)')
% 

%% 
% 
% if expInfo(iX).hemisphere < 0
%         
%         dir{iX} = ... 
%             (nanmean(relMovResps(trialTypes.singleVar.direction{1},:)) - ...
%             nanmean(relMovResps(trialTypes.singleVar.direction{2},:))) ./ ...
%             (nanmean(relMovResps(trialTypes.singleVar.direction{1},:)) + ...
%             nanmean(relMovResps(trialTypes.singleVar.direction{2},:)));
% 
%         side{iX} = ...
%             (nanmean(relStimResps(trialTypes.singleVar.side{1},:)) - ...
%             nanmean(relStimResps(trialTypes.singleVar.side{3},:))) ./ ...
%             (nanmean(relStimResps(trialTypes.singleVar.side{1},:)) + ...
%             nanmean(relStimResps(trialTypes.singleVar.side{3},:)));
%         
%         feed{iX} = ...
%             (nanmean(relRewResps(trialTypes.singleVar.outcome{1},:)) - ...
%             nanmean(relRewResps(trialTypes.singleVar.outcome{2},:))) ./ ...
%             (nanmean(relRewResps(trialTypes.singleVar.outcome{1},:)) + ...
%             nanmean(relRewResps(trialTypes.singleVar.outcome{2},:)));
%         
%         timing{iX} = ...
%             (nanmean(relMovResps(earlyTrials,:)) - ...
%             nanmean(relMovResps(lateTrials,:))) ./ ...
%             (nanmean(relMovResps(earlyTrials,:)) + ...
%             nanmean(relMovResps(lateTrials,:)));
%         
%         hilo{iX} = ...
%             (nanmean(relRewResps(correct_highSide,:)) - ...
%             nanmean(relRewResps(correct_lowSide,:))) ./ ...
%             (nanmean(relRewResps(correct_highSide,:)) + ...
%             nanmean(relRewResps(correct_lowSide,:)));
%     else
%         
%         dir{iX} = ... 
%             (nanmean(relMovResps(trialTypes.singleVar.direction{2},:)) - ...
%             nanmean(relMovResps(trialTypes.singleVar.direction{1},:))) ./ ...
%             (nanmean(relMovResps(trialTypes.singleVar.direction{2},:)) + ...
%             nanmean(relMovResps(trialTypes.singleVar.direction{1},:)));
% 
%         side{iX} = ...
%             (nanmean(relStimResps(trialTypes.singleVar.side{3},:)) - ...
%             nanmean(relStimResps(trialTypes.singleVar.side{1},:))) ./ ...
%             (nanmean(relStimResps(trialTypes.singleVar.side{3},:)) + ...
%             nanmean(relStimResps(trialTypes.singleVar.side{1},:)));
%         
%         feed{iX} = ...
%             (nanmean(relRewResps(trialTypes.singleVar.outcome{1},:)) - ...
%             nanmean(relRewResps(trialTypes.singleVar.outcome{2},:))) ./ ...
%             (nanmean(relRewResps(trialTypes.singleVar.outcome{1},:)) + ...
%             nanmean(relRewResps(trialTypes.singleVar.outcome{2},:)));
%         
%         timing{iX} = ...
%             (nanmean(relMovResps(earlyTrials,:)) - ...
%             nanmean(relMovResps(lateTrials,:))) ./ ...
%             (nanmean(relMovResps(earlyTrials,:)) + ...
%             nanmean(relMovResps(lateTrials,:)));
%         
%         hilo{iX} = ...
%             (nanmean(relRewResps(correct_highSide,:)) - ...
%             nanmean(relRewResps(correct_lowSide,:))) ./ ...
%             (nanmean(relRewResps(correct_highSide,:)) + ...
%             nanmean(relRewResps(correct_lowSide,:)));
% end
    %%
% 
% prevDirIndex = nanmean(baselineResps(prevRightTrials,:)) - nanmean(baselineResps(prevLeftTrials,:));
% baselineStimIndex = nanmean(baselineResps(trialTypes.singleVar.side{3}(1:end-1),:)) - nanmean(baselineResps(trialTypes.singleVar.side{1}(1:end-1),:));
% 
% figure(1);
% minLim = min([prevDirIndex, dirIndex]);
% maxLim = max([prevDirIndex, dirIndex]);
% line([0 0],[minLim maxLim],'LineStyle',':','Color',[.5 .5 .5])
% hold on;
% line([minLim maxLim],[0 0],'LineStyle',':','Color',[.5 .5 .5])
% line([minLim maxLim],[minLim maxLim],'LineStyle',':','Color',[.5 .5 .5])
% scatter(dirIndex,prevDirIndex,'k.')
% xlim([minLim maxLim]);
% ylim([minLim maxLim]);
% axis square
% box off
% set(gca,'tickdir','out')
% xlabel('choice direction preference (t)')
% ylabel('choice direction preference (t–1)')
% 
% 
% %%
% 
% whichNeuron = 16;
% whichNeuron = neuralData(1).stats.bfcH(:,2) == 1 ;
% % whichNeuron = 1:length(neuralData(1).stats.bfcH);
% clear sLmL sLmR
% for c = 1:length(trialTypes.intVar.cb2D.contrast_direction)
% sLmL(c,:) = nanmean(stimResps(trialTypes.intVar.cb2D.contrast_direction{c,1},whichNeuron)-baselineResps(trialTypes.intVar.cb2D.contrast_direction{c,1},whichNeuron),1);
% sLmR(c,:) = nanmean(stimResps(trialTypes.intVar.cb2D.contrast_direction{c,2},whichNeuron)-baselineResps(trialTypes.intVar.cb2D.contrast_direction{c,2},whichNeuron),1);
% end
% 
% 
% 
% figure;
% for c = 1:9
%     subplot(1,9,c)
%     scatter(sLmL(c,:),sLmR(c,:),'k.')
%     hold on
%     line([-.1 .25],[-.1 .25],'LineStyle',':','Color',[.5 .5 .5])
%     xlim([-.05 .25]);
%     ylim([-.05 .25]);
%     axis square
%     box off
% end