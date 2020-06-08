%% initialize experiment details

% expInfo = initExpInfo({{'LEW032'}},{{'2020-02-14',3,[3]}});
expInfo = initExpInfo('LEW031');

%% load data

[expInfo, neuralData, behavioralData] = processExperiment(expInfo, 'matched');
[neuralData] = alignResps(expInfo, neuralData, behavioralData);
[neuralData] = getSignificantActivity(expInfo, behavioralData, neuralData);
combinedNeuralData = combineNeuralData(expInfo, behavioralData, neuralData,'matched');

%% extract some variables to make this code still work 

alignedResps = combinedNeuralData.matched.eta.alignedResps;
eventWindow = combinedNeuralData.matched.eta.eventWindow;
bfcH = combinedNeuralData.matched.stats.bfcH;
pLabels = combinedNeuralData.matched.stats.labels;

%% select cells with the properties you want

plotCells = find(bfcH(:,strcmp(pLabels,'mov')) > 0);
plotCells = 1:size(alignedResps{1},3);

whichETA = 1;
Fs = 0.1;

%% set up trial conditions to compare

contrasts = getUniqueContrasts(expInfo);

%set up trial conditions for hi-L and hi-R blocks
trialConditions{1} = initTrialConditions('highRewardSide','left','responseType','correct','movementTime','late');
trialConditions{2} = initTrialConditions('highRewardSide','right','responseType','correct','movementTime','late');
trialConditions{3} = initTrialConditions('responseType','correct','movementTime','late');
trialLabels{1} = 'highLeft';
trialLabels{2} = 'highRight';
trialLabels{3} = 'all';

contrastConditions{1} = contrasts(contrasts<0);
contrastConditions{2} = contrasts(contrasts>0);
contrastConditions{3} = contrasts(contrasts~=0);
contrastLabels{1} = 'stimLeft';
contrastLabels{2} = 'stimRight';
contrastLabels{3} = 'all';

testTrials = 1:2:size(alignedResps{whichETA},1);
trainTrials = 2:2:size(alignedResps{whichETA},1);

d = 1;
for c = 1:length(contrastConditions)
    for t = 1:length(trialConditions)
        [~, condIdx{d}.all] = selectCondition(expInfo, contrastConditions{c}, behavioralData, trialConditions{t});
        condIdx{d}.test = intersect(testTrials,condIdx{d}.all);
        condIdx{d}.train = intersect(trainTrials,condIdx{d}.all);
        l = sort({contrastLabels{c},trialLabels{t}});
        labels{d} = strcat(l{1},l{2});
        d = d+1;
    end
end

%% report the mean response of each cell to your event, per trial, at a given timepoint
%output is a matrix of size trials x cells

fig = figure;
set(gcf, 'Position', [544  254  1257 716]);
ax1Lim = [0 0.7];
ax2Lim = [-0.5 0.5 0 15];
axLim = [-0.25 0.25];

%initialize some plot values
colors = [0.1 0.7 0.1; 1 .6 0; 0 .4 1; 1 0 0];
mSize = 12;

for timeIdx = 1:length(eventWindow)
    
    %responses of all cells for timepoint timeIdx
    onsetResps = squeeze(alignedResps{whichETA}(:,timeIdx,plotCells));

    % check for cells with nans and remove
    someNaNs = find(isnan(onsetResps));
    badCells = unique(ceil(someNaNs/size(onsetResps,1)));
    onsetResps(:,badCells) = [];
    
    %% compute dot products

    % compute means from 'train' trails
    meanHighLeft = nanmean(onsetResps(condIdx{strcmp(labels,'allhighLeft')}.train,:),1);
    meanHighRight = nanmean(onsetResps(condIdx{strcmp(labels,'allhighRight')}.train,:),1);
    meanStimLeft = nanmean(onsetResps(condIdx{strcmp(labels,'allstimLeft')}.train,:),1);
    meanStimRight = nanmean(onsetResps(condIdx{strcmp(labels,'allstimRight')}.train,:),1);


    % compute dot products between each 'test' trial and mean 'train' trials
    % (split per condition type)

    condType = {'allhighLeft' 'allhighRight' 'allstimLeft' 'allstimRight'};

    fullNorm = 1; % 1 if just angle, 0 if angle and relative magnitude
    for c = 1:length(condType)
        for t = 1:length(condIdx{strcmp(labels,condType{c})}.test)
            iTrial = condIdx{strcmp(labels,condType{c})}.test(t);
            if fullNorm == 1
                tNorm = norm(onsetResps(iTrial,:));
            else
                tNorm = 1;
            end
            dotProd_highSide{c,1}(t,:) = dot(onsetResps(iTrial,:),meanHighLeft)/(tNorm*norm(meanHighLeft));
            dotProd_highSide{c,2}(t,:) = dot(onsetResps(iTrial,:),meanHighRight)/(tNorm*norm(meanHighRight));
            dotProd_stimSide{c,1}(t,:) = dot(onsetResps(iTrial,:),meanStimLeft)/(tNorm*norm(meanStimLeft));
            dotProd_stimSide{c,2}(t,:) = dot(onsetResps(iTrial,:),meanStimRight)/(tNorm*norm(meanStimRight));
        end
    end

    % compute diagonal histograms
    for d = 1:4
        [density_high(:,d),value_high(:,d)] = ksdensity(dotProd_highSide{d,1}-dotProd_highSide{d,2},'Bandwidth', max(ax1Lim)/20);
        [density_stim(:,d),value_stim(:,d)] = ksdensity(dotProd_stimSide{d,1}-dotProd_stimSide{d,2},'Bandwidth', max(ax1Lim)/20);
    end
    
    %compute (hiL - hiR) and (stimL - stimR) means from 'train' trials
    meanResps_reward = mean(onsetResps(condIdx{strcmp(labels,'allhighLeft')}.train,:),1) - ...
                            mean(onsetResps(condIdx{strcmp(labels,'allhighRight')}.train,:),1);
    meanResps_stim = mean(onsetResps(condIdx{strcmp(labels,'allstimLeft')}.train,:),1) - ...
                            mean(onsetResps(condIdx{strcmp(labels,'allstimRight')}.train,:),1);
                        
    % compute dot products between each 'test' trial and mean 'train' trials
    % (split per condition type)
    projDot_reward = dot(onsetResps, repmat(meanResps_reward,size(onsetResps,1),1),2)./(vectorNorm(onsetResps,2).*vectorNorm(repmat(meanResps_reward,size(onsetResps,1),1),2));
    projDot_stim = dot(onsetResps, repmat(meanResps_stim,size(onsetResps,1),1),2)./(vectorNorm(onsetResps,2).*vectorNorm(repmat(meanResps_stim,size(onsetResps,1),1),2));
    
    % compute the mean dot product per trial type
    trajectoryRewards_mean(timeIdx,1,1) = mean(projDot_stim(condIdx{strcmp(labels,'allhighLeft')}.test));
    trajectoryRewards_mean(timeIdx,2,1) = mean(projDot_reward(condIdx{strcmp(labels,'allhighLeft')}.test));
    trajectoryRewards_mean(timeIdx,1,2) = mean(projDot_stim(condIdx{strcmp(labels,'allhighRight')}.test));
    trajectoryRewards_mean(timeIdx,2,2) = mean(projDot_reward(condIdx{strcmp(labels,'allhighRight')}.test));
    
    trajectoryRewards_sem(timeIdx,1,1) = std(projDot_stim(condIdx{strcmp(labels,'allhighLeft')}.test))/sqrt(length(projDot_reward));
    trajectoryRewards_sem(timeIdx,2,1) = std(projDot_reward(condIdx{strcmp(labels,'allhighLeft')}.test))/sqrt(length(projDot_reward));
    trajectoryRewards_sem(timeIdx,1,2) = std(projDot_stim(condIdx{strcmp(labels,'allhighRight')}.test))/sqrt(length(projDot_reward));
    trajectoryRewards_sem(timeIdx,2,2) = std(projDot_reward(condIdx{strcmp(labels,'allhighRight')}.test))/sqrt(length(projDot_reward));

    trajectoryStims_mean(timeIdx,1,1) = mean(projDot_stim(condIdx{strcmp(labels,'allhighLeft')}.test));
    trajectoryStims_mean(timeIdx,2,1) = mean(projDot_reward(condIdx{strcmp(labels,'allstimLeft')}.test));
    trajectoryStims_mean(timeIdx,1,2) = mean(projDot_stim(condIdx{strcmp(labels,'allstimRight')}.test));
    trajectoryStims_mean(timeIdx,2,2) = mean(projDot_reward(condIdx{strcmp(labels,'allstimRight')}.test));

    trajectoryStims_sem(timeIdx,1,1) = std(projDot_stim(condIdx{strcmp(labels,'allhighLeft')}.test))/sqrt(length(projDot_stim));
    trajectoryStims_sem(timeIdx,2,1) = std(projDot_reward(condIdx{strcmp(labels,'allstimLeft')}.test))/sqrt(length(projDot_stim));
    trajectoryStims_sem(timeIdx,1,2) = std(projDot_stim(condIdx{strcmp(labels,'allstimRight')}.test))/sqrt(length(projDot_stim));
    trajectoryStims_sem(timeIdx,2,2) = std(projDot_reward(condIdx{strcmp(labels,'allstimRight')}.test))/sqrt(length(projDot_stim));
    
  %% plot 
  for s = 1:6
      subplot(2,3,s)
      cla;
  end
    % configure the subplots before adding any data
    for s = [1 4]
        subplot(2,3,s)
        ax1 = gca;
        hold on;
        axis square;
        plotLine = line([ax1Lim(1) ax1Lim(2)],[ax1Lim(1) ax1Lim(2)]);
        set(plotLine,'Color', [.5 .5 .5],'LineStyle','--')
        axis([ax1Lim(1) ax1Lim(2) ax1Lim(1) ax1Lim(2)]);
        hold on;
        ax1.TickDir = 'out';
        if s == 1 
            title('trials split by reward block');
        elseif s == 3
            title('trials split by stimulus side');
        end
    end

    %plot all trials compared to high-L and high-R mean activity
    %(sorted by reward or stim)
    subplot(2,3,1)
    hold on
    plotHighL = scatter(dotProd_highSide{1,1},dotProd_highSide{1,2},mSize,'go');
    set(plotHighL,'MarkerEdgeColor',colors(1,:), 'Marker','o','MarkerFaceColor',colors(1,:),'MarkerFaceAlpha',.4)
    plotHighR = scatter(dotProd_highSide{2,1},dotProd_highSide{2,2},mSize,'ro');
    set(plotHighR,'MarkerEdgeColor',colors(2,:), 'Marker','o','MarkerFaceColor',colors(2,:),'MarkerFaceAlpha',.4)
    ylabel('correlation to mean high-right activity')
    xlabel('correlation to mean high-left activity')

    %plot all trials compared to stim-L and stim-R mean activity
    %(sorted by reward or stim)
    subplot(2,3,4)
    plotStimL = scatter(dotProd_stimSide{3,1},dotProd_stimSide{3,2},mSize,'go');
    set(plotStimL,'MarkerEdgeColor',colors(3,:), 'Marker','o','MarkerFaceColor',colors(3,:),'MarkerFaceAlpha',.4)
    plotStimR = scatter(dotProd_stimSide{4,1},dotProd_stimSide{4,2},mSize,'ro');
    set(plotStimR,'MarkerEdgeColor',colors(4,:), 'Marker','o','MarkerFaceColor',colors(4,:),'MarkerFaceAlpha',.4)       
    ylabel('correlation to mean stim-right activity')
    xlabel('correlation to mean stim-left activity')

      %plot density compared to reward
    subplot(2,3,2)
    hold on;
    pHighL = plot(value_high(:,1),density_high(:,1));
    set(pHighL, 'LineWidth', 2, 'Color', colors(1,:));
    pHighL.Color(4) = 0.3;
    pHighR = plot(value_high(:,2),density_high(:,2));
    set(pHighR, 'LineWidth', 2, 'Color', colors(2,:));
    pHighR.Color(4) = 0.3;
    title('trials split by reward block');
    axis([ax2Lim(1) ax2Lim(2) ax2Lim(3) ax2Lim(4)]);
    xlabel('corr_{hiLeft} - corr_{hiRight}')
    ylabel('number of trials')

    %plot density compared to stimulus
    subplot(2,3,5)
    hold on;
    pHighL = plot(value_stim(:,3),density_stim(:,3));
    set(pHighL, 'LineWidth', 2, 'Color', colors(3,:));
    pHighL.Color(4) = 0.3;
    pHighR = plot(value_stim(:,4),density_stim(:,4));
    set(pHighR, 'LineWidth', 2, 'Color', colors(4,:));
    pHighR.Color(4) = 0.3;
    title('trials split by stimulus side');
    axis([ax2Lim(1) ax2Lim(2) ax2Lim(3) ax2Lim(4)]);
    xlabel('corr_{stimLeft} - corr_{stimRight}')
    ylabel('number of trials')
    
    subplot(2,3,3);
    hold on;
    p3 = plot(trajectoryRewards_mean(1:timeIdx,1,1),trajectoryRewards_mean(1:timeIdx,2,1),'Color',colors(1,:),'LineWidth',1);
    p3.Color(4) = 0.5;
    p4 = plot(trajectoryRewards_mean(1:timeIdx,1,2),trajectoryRewards_mean(1:timeIdx,2,2),'Color',colors(2,:),'LineWidth',1);
    p4.Color(4) = 0.5;
    if timeIdx >= 21
        plot(trajectoryRewards_mean(21,1,1),trajectoryRewards_mean(21,2,1),'ko','MarkerFaceColor',colors(1,:));
        plot(trajectoryRewards_mean(21,1,2),trajectoryRewards_mean(21,2,2),'ko','MarkerFaceColor',colors(2,:));
    end
    xlim(axLim);
    ylim(axLim);
    axis square
    ylabel('highL - highR')
    % legend({'left high','right high'})
    xlabel('stimL - stimR')
    
    subplot(2,3,6)
    hold on;
    p1 = plot(trajectoryStims_mean(1:timeIdx,1,1),trajectoryStims_mean(1:timeIdx,2,1),'Color',colors(3,:),'LineWidth',1);
    p1.Color(4) = 0.5;
    p2 = plot(trajectoryStims_mean(1:timeIdx,1,2),trajectoryStims_mean(1:timeIdx,2,2),'Color',colors(4,:),'LineWidth',1);
    p2.Color(4) = 0.5;
    if timeIdx >= 21
        plot(trajectoryStims_mean(21,1,1),trajectoryStims_mean(21,2,1),'ko','MarkerFaceColor',colors(3,:));
        plot(trajectoryStims_mean(21,1,2),trajectoryStims_mean(21,2,2),'ko','MarkerFaceColor',colors(4,:));
    end
    xlim(axLim);
    ylim(axLim);
    axis square
    ylabel('highL - highR')
    % legend({'left stim','right stim'})
    xlabel('stimL - stimR')
    pause
%     printfig(gcf,strcat('trajectory',num2str(timeIdx)))
    
end

%