function plotWheelPseudorasters(expInfo, behavioralData, neuralData, moveTime)

stimOnTimes = behavioralData.eventTimes(1).daqTime;
feedbackTimes = behavioralData.eventTimes(5).daqTime;
choiceDir = expInfo.block.events.responseValues;
contrasts = expInfo.block.events.contrastValues(1:end-1);
uniqueContrasts = getUniqueContrasts(expInfo);

feedbackLatencies = behavioralData.eventTimes(5).daqTime - behavioralData.eventTimes(1).daqTime;
firstMoveLatencies = behavioralData.eventTimes(5).daqTime - behavioralData.wheelMoves.epochs(5).onsetTimes;

plotWindow = -.5:0.001:3.5;
eventIdx = find(plotWindow == 0);
fromStart = plotWindow(1)*1000;
toEnd = plotWindow(end)*1000;
%
for t = 1:length(behavioralData.wheelMoves.traces.time)
    traceTime = behavioralData.wheelMoves.traces.time{1,t}  - behavioralData.eventTimes(1).daqTime(t);
    tracePos = behavioralData.wheelMoves.traces.pos{1,t};
    velTime = traceTime;
    vel = [0 movmean(diff(tracePos),100)];
    zeroIdx = find(velTime == 0);
    startIdx = zeroIdx + fromStart;
    endIdx = zeroIdx + toEnd;
    
    try
        truncTime = velTime(startIdx:endIdx);
        truncVel = vel(startIdx:endIdx);
    catch
        try
            padding = endIdx - length(traceTime);
            truncTime = [velTime(startIdx:end), nan(1,padding)];
            truncVel = [vel(startIdx:end), nan(1,padding)];
        catch
            try
                padding = -startIdx+1;
                truncTime = [nan(1,padding), velTime(1:endIdx)];
                truncVel = [nan(1,padding),vel(1:endIdx)];
            catch
                s_pad = -startIdx+1;
                e_pad = endIdx - length(traceTime);
                truncTime = [nan(1,s_pad),velTime(1:end),nan(1,e_pad)];
                truncVel = [nan(1,s_pad),vel(1:end),nan(1,e_pad)];
            end
        end
    end
    
    bigVel(t,:) = truncVel;
end

trialConditions{1} = initTrialConditions('movementTime','early');
trialConditions{2} = initTrialConditions('movementTime','late','specificRTs',[0 3]);
trialConditions{3} = initTrialConditions('specificRTs',[3 Inf]);

if strcmp(moveTime,'early')
    T = 1;
elseif strcmp(moveTime,'late')
    T = 2;
end

%%
fig = figure; 
set(gcf,'position', [1250 440 760 1180])
for c = 1:length(uniqueContrasts)

    [whichTrials] = selectCondition(expInfo, uniqueContrasts(c), behavioralData, trialConditions{T});
    leftVels = bigVel(choiceDir == -1 & whichTrials,:);
    rightVels = bigVel(choiceDir == 1 & whichTrials,:);
    if T == 1
        [~, leftSortIdx] = sort(firstMoveLatencies(choiceDir == -1 & whichTrials));
        [~, rightSortIdx] = sort(firstMoveLatencies(choiceDir == 1 & whichTrials));
    else
        [~, leftSortIdx] = sort(feedbackLatencies(choiceDir == -1 & whichTrials));
        [~, rightSortIdx] = sort(feedbackLatencies(choiceDir == 1 & whichTrials));
    end
    
    meanTrajectoryLeft(c,:) = nanmean(-leftVels(leftSortIdx,:),1);
    meanTrajectoryRight(c,:) = nanmean(-rightVels(rightSortIdx,:),1);
    

    subplot(length(uniqueContrasts),2,2*c -1);
    imagesc(plotWindow, 1:length(leftSortIdx), -leftVels(leftSortIdx,:),'AlphaData',~isnan(leftVels));
    hold on;
    line([0 0],[0 100],'LineStyle','--','Color',[.5 .5 .5])
    caxis([-.25 .25]);
    xlim([-.5 2])
    ylabel(num2str(uniqueContrasts(c)))
    set(gca,'tickdir','out')
    box off
    if c < length(uniqueContrasts)
        set(gca, 'XTickLabels', {''});
    end
    if c == 1
        title('chose left')
    end
    
    subplot(length(uniqueContrasts),2,2*c);
    
    imagesc(plotWindow, 1:length(rightSortIdx), -rightVels(rightSortIdx,:),'AlphaData',~isnan(rightVels));
    hold on;
    line([0 0],[0 100],'LineStyle','--','Color',[.5 .5 .5])
    colormap(BlueWhiteRed);
    caxis([-.25 .25]);
    xlim([-.5 2])
    set(gca,'tickdir','out')
    box off
    if c == 1
        title('chose right')
    end
    if c < length(uniqueContrasts)
        set(gca, 'XTickLabels', {''});
    end
end

if expInfo.hemisphere < 0
    tmp = meanTrajectoryRight;
    meanTrajectoryRight = -flipud(meanTrajectoryLeft);
    meanTrajectoryLeft = -flipud(tmp);
end

% close(fig)
%%

zeroIdx = find(uniqueContrasts == 0);
walkup = length(uniqueContrasts) - zeroIdx;
walkback = zeroIdx - 1;
% allColors = [.25 0 0;.5 0 0 ;1 0 0;.8 .45 .45;.75 .75 .75;.55 .55 .55;.35 .35 .35;.15 .15 .15;0 0 0];
allColors = [0 0 .5; 0 0 1; 0 .4 1; .6 .8 1; .75 .75 .75; .8 .45 .45; 1 0 0; .5 0 0; .25 0 0];

zeroGray = find(allColors(:,1) == .75);
colors = allColors(zeroGray-walkback:zeroGray + walkup,:);

figure;
for c = 1:length(uniqueContrasts)
    if c < 5
        subplot(2,2,1)
        plot(plotWindow, -meanTrajectoryLeft(c,:),'LineWidth',2,'Color', colors(c,:));
        line([0 0],[-.2 .2],'LineStyle','--','Color',[.5 .5 .5])
    elseif c > 5
        subplot(2,2,4)
        plot(plotWindow, -meanTrajectoryLeft(c,:),'LineWidth',2,'Color', colors(c,:));
        line([0 0],[-.2 .2],'LineStyle','--','Color',[.5 .5 .5])
    end
    hold on;
    xlim([-.5 2]);
    ylim([-.01 .15]);
    set(gca,'tickdir','out')
    box off
end
subplot(2,2,1);
% title('chose left')
for c = 1:length(uniqueContrasts)
    if c < 5
        subplot(2,2,2)
        plot(plotWindow, meanTrajectoryRight(c,:),'LineWidth',2,'Color', colors(c,:));
        line([0 0],[-.2 .2],'LineStyle','--','Color',[.5 .5 .5])
    elseif c > 5
        subplot(2,2,3)
        plot(plotWindow, meanTrajectoryRight(c,:),'LineWidth',2,'Color', colors(c,:));
        line([0 0],[-.2 .2],'LineStyle','--','Color',[.5 .5 .5])
    end
    
    hold on;
    xlim([-.5 2]);
    ylim([-.01 .15]);
    set(gca,'tickdir','out')
    box off
end
subplot(2,2,2);
% title('chose right')

%%
fig = figure; 
set(gcf,'position', [1250 440 760 1180])
bins = 0:.1:3;
maxy = 25;
for c = 1:length(uniqueContrasts)
    
    [whichTrials] = selectCondition(expInfo, uniqueContrasts(c), behavioralData, trialConditions{T});
    lRTs = behavioralData.wheelMoves.epochs(5).onsetTimes(choiceDir == -1 & whichTrials == 1) - ...
        behavioralData.eventTimes(1).daqTime(choiceDir == -1 & whichTrials == 1);
    rRTs = behavioralData.wheelMoves.epochs(5).onsetTimes(choiceDir == 1 & whichTrials == 1) - ...
        behavioralData.eventTimes(1).daqTime(choiceDir == 1 & whichTrials == 1);
    medianLRTs = nanmedian(lRTs);
    medianRRTs = nanmedian(rRTs);
    subplot(length(uniqueContrasts),2,2*c -1);
    histogram(lRTs, bins)
    hold on;
    line([0 0],[0 100],'LineStyle','--','Color',[.5 .5 .5])
    line([medianLRTs medianLRTs],[0 100],'LineStyle','-','LineWidth',2,'Color',[.25 .25 .25])
    ylim([0 maxy])
    xlim([-.5 2])
    ylabel(num2str(uniqueContrasts(c)))
    set(gca,'tickdir','out')
    box off
    if c < length(uniqueContrasts)
        set(gca, 'XTickLabels', {''});
    end
    if c == 1
        title('chose left')
    end
    subplot(length(uniqueContrasts),2,2*c);
    histogram(rRTs, bins)
    hold on;
    line([0 0],[0 100],'LineStyle','--','Color',[.5 .5 .5])
    line([medianRRTs medianRRTs],[0 100],'LineStyle','-','LineWidth',2,'Color',[.25 .25 .25])
    ylim([0 maxy])
    xlim([-.5 2])
    set(gca,'tickdir','out')
    box off
    if c < length(uniqueContrasts)
        set(gca, 'XTickLabels', {''});
    end
    if c == 1
        title('chose right')
    end
end

% close(fig)
%%

whichCells = find(neuralData.stats.bfcH(:,strcmp(neuralData.stats.labels,'rightMov')) > 0);

for c = 1:length(uniqueContrasts)
    [whichTrials] = selectCondition(expInfo, uniqueContrasts(c), behavioralData, trialConditions{T});
    popResp = squeeze(nanmean(neuralData.eta.alignedResps{1}((choiceDir == 1 & whichTrials == 1),:,whichCells),1));
    meanPopResp(c,:) = nanmean(popResp,2);
end

trimmedIdx = find(plotWindow >= -.5 & plotWindow <=2);
trimWindow = plotWindow(trimmedIdx);
trimTraj = meanTrajectoryRight(:,trimmedIdx);
trimmedNeurIdx = find(neuralData.eta.eventWindow >= -.5 & neuralData.eta.eventWindow <=2);
trimNeurWindow = neuralData.eta.eventWindow(trimmedNeurIdx);
trimNeur = meanPopResp(:,trimmedNeurIdx);

for c = 1:length(uniqueContrasts)
    meanPopResp_int(c,:) = interp1(trimNeurWindow,trimNeur(c,:),trimWindow,'linear');
end

