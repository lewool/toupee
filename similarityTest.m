%% overview

% this script plots ETAs aligned to any event you want, and compares them
% between two conditions of your choosing (high/low reward, t-1 left or 
% right, correct/incorrect, etc.) Should work with most lewieWorld_B2AFC
% expDefs

% 2018-11-20: LEW created (plotCRF_stimulusHistory.m)
% 2019-05-20: Changed name; modified to more easily change condition comparison
% and adapt to varied contrast conditions across subjects 
% 2019-07-02: moved to git repository 'toupee'; extracted cell type chooser
% to separate function call
% 2020-01-14 streamlining


%%
clear all
mouseList = {{'LEW007'},{'LEW008'},{'LEW008'},{'LEW008'},{'LEW031'},{'LEW032'}};
expList = {{'2018-10-09',1,[1]},{'2019-01-29',1,[1]},{'2019-02-07',1,[1]},{'2019-02-12',1,[1]},{'2020-02-03',1,[1]},{'2020-02-03',1,[1]}}; 

mouseList = ...
    {{'LEW005'},...
    {'LEW006'},...
    {'LEW006'},...
    {'LEW013'},...
    {'LEW013'},...
    {'LEW015'},...
    {'LEW015'}};

expList = ...
    {{'2018-06-10',2,[2 3]},...
    {'2018-06-14',1,[1 2]},...
    {'2018-06-15',1,[1 2]},...
    {'2019-03-26',1,1},...
    {'2019-03-27',1,1},...
    {'2019-03-21',1,1},...
    {'2019-04-12',1,1}};

mouseList = {{'LEW031'},{'LEW032'}};
expList = {{'2020-02-03',1,[1]},{'2020-02-03',1,[1]}};
whichETA = 2;
allHistos_high = {};
allHistos_stim = {};

%%
for m = 1:length(mouseList)
    
clearvars -except m mouseList expList allHistos_high allHistos_stim whichETA
disp(m);

%% load experiment details

expInfo = initExpInfo(mouseList(m),expList(m));

%% load data

expInfo = data.loadExpData(expInfo);

%% get event timings and wheel trajectories

[eventTimes, wheelTrajectories] = getEventTimes(expInfo, {'stimulusOnTimes' 'interactiveOnTimes' 'stimulusOffTimes'});

%% load traces

[allFcell, expInfo] = loadCellData(expInfo);
[cellResps, respTimes] = getCellResps(expInfo, allFcell);
cellResps = zscore(cellResps);

%% align calcium traces to the event you want

% cut the trace into trial-by-trial traces, aligned to a particular event
events = {'stimulusOnTimes' 'prestimulusQuiescenceEndTimes' 'feedbackTimes'};
alignedResps = cell(1,length(events));
for e = 1:length(events)
    [alignedResps{e}, eventWindow] = alignResps(expInfo, cellResps, respTimes, eventTimes, events{e});
end


%% select cells with the properties you want

plotCells = chooseCellType('all', expInfo, cellResps, respTimes, eventTimes, 0.1);


%% report the mean response of each cell to your event, per trial, at a given timepoint
%output is a matrix of size trials x cells
Fs = 0.1;
timepoint = 0;
eventIdx = find(eventWindow == 0);
timeIdx = eventIdx + round(timepoint/Fs);

onsetResps = squeeze(alignedResps{whichETA}(:,timeIdx,plotCells));

% check for cells with nans and remove
someNaNs = find(isnan(onsetResps));
badCells = unique(ceil(someNaNs/size(onsetResps,1)));
onsetResps(:,badCells) = [];

%% set up trial conditions to compare

contrasts = unique(expInfo.block.events.contrastValues);

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

testTrials = 1:2:size(onsetResps,1);
trainTrials = 2:2:size(onsetResps,1);

d = 1;
for c = 1:length(contrastConditions)
    for t = 1:length(trialConditions)
        [~, condIdx{d}.all] = selectCondition(expInfo.block, contrastConditions{c}, eventTimes, trialConditions{t});
        condIdx{d}.test = intersect(testTrials,condIdx{d}.all);
        condIdx{d}.train = intersect(trainTrials,condIdx{d}.all);
        l = sort({contrastLabels{c},trialLabels{t}});
        labels{d} = strcat(l{1},l{2});
        d = d+1;
    end
end

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

% assign axis limits for later
ax1Lim = 1.2 * [min([cell2mat(dotProd_highSide(:)); cell2mat(dotProd_stimSide(:))]) ...
                max([cell2mat(dotProd_highSide(:)); cell2mat(dotProd_stimSide(:))])]; 

% compute diagonal histograms
for d = 1:4
    [density_high(:,d),value_high(:,d)] = ksdensity(dotProd_highSide{d,1}-dotProd_highSide{d,2},'Bandwidth', max(ax1Lim)/20);
    [density_stim(:,d),value_stim(:,d)] = ksdensity(dotProd_stimSide{d,1}-dotProd_stimSide{d,2},'Bandwidth', max(ax1Lim)/20);
end

% assign axis limits for later
ax2Lim = 1.2 * [-max([cell2mat(dotProd_highSide(:)); cell2mat(dotProd_stimSide(:))]) ...
                max([cell2mat(dotProd_highSide(:)); cell2mat(dotProd_stimSide(:))]) ...
                0 max(max([density_stim(:) density_high(:)]))];

% add histograms to global list
for h = 1:4
    allHistos_high{h, length(allHistos_high)+1} = dotProd_highSide{h,1}-dotProd_highSide{h,2};
    allHistos_stim{h, length(allHistos_stim)+1} = dotProd_highSide{h,1}-dotProd_highSide{h,2};
end

%% single-timepoint scatterplot

%initialize some plot values
colors = [0.1 0.7 0.1; 1 .6 0; 0 .4 1; 1 0 0];
mSize = 12;

for f = 1:min([length(mouseList) 2])
    % first plot single session
    if f == 1
        figure;
        set(gcf, 'Position', [186  554  850 700], 'Name', strcat(expInfo.mouseName,'_', expInfo.expDate));
    % then plot all sessions
    else
        figure(994);
        set(gcf, 'Position', [186  554  850 700], 'Name', 'all sessions');
    end
    
    % configure the subplots before adding any data
    for s = 1:4
        subplot(2,2,s)
        ax1 = gca;
        hold on;
        axis square;
        plotLine = line([ax1Lim(1) ax1Lim(2)],[ax1Lim(1) ax1Lim(2)]);
        set(plotLine,'Color', [.5 .5 .5],'LineStyle','--')
        axis([ax1Lim(1) ax1Lim(2) ax1Lim(1) ax1Lim(2)]);
        hold on;
        ax1.TickDir = 'out';
        if s < 3
            title('trials split by reward block');
        else
            title('trials split by stimulus side');
        end
    end
    
    % plot
    for s = 1:2
        %plot all trials compared to high-L and high-R mean activity
        %(sorted by reward or stim)
        subplot(2,2,2*s-1)
        plotHighL = scatter(dotProd_highSide{2*s-1,1},dotProd_highSide{2*s-1,2},mSize,'go');
        set(plotHighL,'MarkerEdgeColor',colors(2*s-1,:), 'Marker','o','MarkerFaceColor',colors(2*s-1,:),'MarkerFaceAlpha',.4)
        plotHighR = scatter(dotProd_highSide{2*s,1},dotProd_highSide{2*s,2},mSize,'ro');
        set(plotHighR,'MarkerEdgeColor',colors(2*s,:), 'Marker','o','MarkerFaceColor',colors(2*s,:),'MarkerFaceAlpha',.4)
        ylabel('similarity to high-right mean activity')
        xlabel('similarity to high-left mean activity')
            
        %plot all trials compared to stim-L and stim-R mean activity
        %(sorted by reward or stim)
        subplot(2,2,2*s)
        plotStimL = scatter(dotProd_stimSide{2*s-1,1},dotProd_stimSide{2*s-1,2},mSize,'go');
        set(plotStimL,'MarkerEdgeColor',colors(2*s-1,:), 'Marker','o','MarkerFaceColor',colors(2*s-1,:),'MarkerFaceAlpha',.4)
        plotStimR = scatter(dotProd_stimSide{2*s,1},dotProd_stimSide{2*s,2},mSize,'ro');
        set(plotStimR,'MarkerEdgeColor',colors(2*s,:), 'Marker','o','MarkerFaceColor',colors(2*s,:),'MarkerFaceAlpha',.4)       
        ylabel('similarity to stim-right mean activity')
        xlabel('similarity to stim-left mean activity')
    end
    
    % plot the diagonal histos for 4 conditions
    if f == 1
        % first plot single session
        figure;
        set(gcf, 'Position', [186  554  850 700], 'Name', strcat(expInfo.mouseName,'_', expInfo.expDate));
    else
        % then plot all sessions
        figure(1004);
        set(gcf, 'Position', [186  554  850 700], 'Name', 'all sessions');
    end
    
    for s = 1:2
        if s == 1
            plotValue = value_high;
            plotDensity = density_high;
        else
            plotValue = value_stim;
            plotDensity = density_stim;
        end
        
        %plot density compared to reward
        subplot(2,2,s)
        hold on;
        pHighL = plot(plotValue(:,1),plotDensity(:,1));
        set(pHighL, 'LineWidth', 2, 'Color', colors(1,:));
        pHighL.Color(4) = 0.3;
        pHighR = plot(plotValue(:,2),plotDensity(:,2));
        set(pHighR, 'LineWidth', 2, 'Color', colors(2,:));
        pHighR.Color(4) = 0.3;
        title('trials split by reward block');
        axis([ax2Lim(1) ax2Lim(2) ax2Lim(3) ax2Lim(4)]);
        
        %plot density compared to stimulus
        subplot(2,2,s+2)
        hold on;
        pHighL = plot(plotValue(:,3),plotDensity(:,3));
        set(pHighL, 'LineWidth', 2, 'Color', colors(3,:));
        pHighL.Color(4) = 0.3;
        pHighR = plot(plotValue(:,4),plotDensity(:,4));
        set(pHighR, 'LineWidth', 2, 'Color', colors(4,:));
        pHighR.Color(4) = 0.3;
        title('trials split by stimulus side');
        axis([ax2Lim(1) ax2Lim(2) ax2Lim(3) ax2Lim(4)]);
    end
    
end

%% trajectories

%compute trial trajectories
for p = 1:length(eventWindow)
    
    % extract the response of each cell at timepoint p (relative to event time)
    projResps = squeeze(alignedResps{whichETA}(:,p,plotCells));
    
    %get rid of bad cells
    someNaNs = find(isnan(projResps));
    badCells = unique(ceil(someNaNs/size(projResps,1)));
    projResps(:,badCells) = [];
    
    %compute (hiL - hiR) and (stimL - stimR) means from 'train' trials
    meanResps_reward = mean(projResps(condIdx{strcmp(labels,'allhighLeft')}.train,:),1) - ...
                            mean(projResps(condIdx{strcmp(labels,'allhighRight')}.train,:),1);
    meanResps_stim = mean(projResps(condIdx{strcmp(labels,'allstimLeft')}.train,:),1) - ...
                            mean(projResps(condIdx{strcmp(labels,'allstimRight')}.train,:),1);
                        
    % compute dot products between each 'test' trial and mean 'train' trials
    % (split per condition type)
    projDot_reward = dot(projResps, repmat(meanResps_reward,size(projResps,1),1),2)./(vectorNorm(projResps,2).*vectorNorm(repmat(meanResps_reward,size(projResps,1),1),2));
    projDot_stim = dot(projResps, repmat(meanResps_stim,size(projResps,1),1),2)./(vectorNorm(projResps,2).*vectorNorm(repmat(meanResps_stim,size(projResps,1),1),2));
    
    % compute the mean dot product per trial type
    trajectoryRewards_mean(p,1,1) = mean(projDot_stim(condIdx{strcmp(labels,'allhighLeft')}.test));
    trajectoryRewards_mean(p,2,1) = mean(projDot_reward(condIdx{strcmp(labels,'allhighLeft')}.test));
    trajectoryRewards_mean(p,1,2) = mean(projDot_stim(condIdx{strcmp(labels,'allhighRight')}.test));
    trajectoryRewards_mean(p,2,2) = mean(projDot_reward(condIdx{strcmp(labels,'allhighRight')}.test));
    
    trajectoryRewards_sem(p,1,1) = std(projDot_stim(condIdx{strcmp(labels,'allhighLeft')}.test))/sqrt(length(projDot_reward));
    trajectoryRewards_sem(p,2,1) = std(projDot_reward(condIdx{strcmp(labels,'allhighLeft')}.test))/sqrt(length(projDot_reward));
    trajectoryRewards_sem(p,1,2) = std(projDot_stim(condIdx{strcmp(labels,'allhighRight')}.test))/sqrt(length(projDot_reward));
    trajectoryRewards_sem(p,2,2) = std(projDot_reward(condIdx{strcmp(labels,'allhighRight')}.test))/sqrt(length(projDot_reward));

    trajectoryStims_mean(p,1,1) = mean(projDot_stim(condIdx{strcmp(labels,'allhighLeft')}.test));
    trajectoryStims_mean(p,2,1) = mean(projDot_reward(condIdx{strcmp(labels,'allstimLeft')}.test));
    trajectoryStims_mean(p,1,2) = mean(projDot_stim(condIdx{strcmp(labels,'allstimRight')}.test));
    trajectoryStims_mean(p,2,2) = mean(projDot_reward(condIdx{strcmp(labels,'allstimRight')}.test));

    trajectoryStims_sem(p,1,1) = std(projDot_stim(condIdx{strcmp(labels,'allhighLeft')}.test))/sqrt(length(projDot_stim));
    trajectoryStims_sem(p,2,1) = std(projDot_reward(condIdx{strcmp(labels,'allstimLeft')}.test))/sqrt(length(projDot_stim));
    trajectoryStims_sem(p,1,2) = std(projDot_stim(condIdx{strcmp(labels,'allstimRight')}.test))/sqrt(length(projDot_stim));
    trajectoryStims_sem(p,2,2) = std(projDot_reward(condIdx{strcmp(labels,'allstimRight')}.test))/sqrt(length(projDot_stim));

end

% set axis limits
axMax = max([max(abs(trajectoryRewards_mean(:))) max(abs(trajectoryStims_mean(:)))]);
axLim = [-axMax axMax] *1.2;

% plot
for f = 1:min([length(mouseList) 2])
    % first plot single session
    if f == 1
        figure;
        set(gcf,'Position',[100 100 1025 420], 'Name', strcat(expInfo.mouseName,'_', expInfo.expDate));
    % then plot all sessions
    else
        figure(1005);
        set(gcf,'Position',[100 100 1025 420], 'Name', 'all sessions');
    end
    
    subplot(1,2,1);
    hold on;
    p3 = plot(trajectoryRewards_mean(:,1,1),trajectoryRewards_mean(:,2,1),'Color',colors(1,:),'LineWidth',1);
    p3.Color(4) = 0.5;
    plot(trajectoryRewards_mean(eventIdx,1,1),trajectoryRewards_mean(eventIdx,2,1),'ko','MarkerFaceColor',colors(1,:));
    p4 = plot(trajectoryRewards_mean(:,1,2),trajectoryRewards_mean(:,2,2),'Color',colors(2,:),'LineWidth',1);
    p4.Color(4) = 0.5;
    plot(trajectoryRewards_mean(eventIdx,1,2),trajectoryRewards_mean(eventIdx,2,2),'ko','MarkerFaceColor',colors(2,:));
    xlim(axLim);
    ylim(axLim);
    axis square
    ylabel('highL - highR')
    % legend({'left high','right high'})
    xlabel('stimL - stimR')
    
    subplot(1,2,2)
    hold on;
    p1 = plot(trajectoryStims_mean(:,1,1),trajectoryStims_mean(:,2,1),'Color',colors(3,:),'LineWidth',1);
    p1.Color(4) = 0.5;
    plot(trajectoryStims_mean(eventIdx,1,1),trajectoryStims_mean(eventIdx,2,1),'ko','MarkerFaceColor',colors(3,:));
    p2 = plot(trajectoryStims_mean(:,1,2),trajectoryStims_mean(:,2,2),'Color',colors(4,:),'LineWidth',1);
    p2.Color(4) = 0.5;
    plot(trajectoryStims_mean(eventIdx,1,2),trajectoryStims_mean(eventIdx,2,2),'ko','MarkerFaceColor',colors(4,:));
    xlim(axLim);
    ylim(axLim);
    axis square
    ylabel('highL - highR')
    % legend({'left stim','right stim'})
    xlabel('stimL - stimR')

end
% 
% figure;
% hold on
% set(gcf,'Position',[100 100 1025 420]);
% set(gcf,'renderer','Painters');
% subplot(1,2,2)
% hold on;
% p1 = plot(trajectoryStims_mean(:,1,1),trajectoryStims_mean(:,2,1),'Color',colors(3,:),'LineWidth',1);
% p2 = plot(trajectoryStims_mean(:,1,2),trajectoryStims_mean(:,2,2),'Color',colors(4,:),'LineWidth',1);
% xlim(axLim);
% ylim(axLim);
% axis square
% ylabel('highL - highR')
% % legend({'left stim','right stim'})
% xlabel('stimL - stimR')
% 
% subplot(1,2,1);
% hold on;
% p3 = plot(trajectoryRewards_mean(:,1,1),trajectoryRewards_mean(:,2,1),'Color',colors(1,:),'LineWidth',1);
% p4 = plot(trajectoryRewards_mean(:,1,2),trajectoryRewards_mean(:,2,2),'Color',colors(2,:),'LineWidth',1);
% xlim(axLim);
% ylim(axLim);
% axis square
% ylabel('highL - highR')
% % legend({'left high','right high'})
% xlabel('stimL - stimR')

%% make a trajectory movie

% leftRewardsTrajectory_int(:,1) = interp1(eventWindow, leftRewardsTrajectory(:,1)', linspace(eventWindow(1),eventWindow(end), 50));
% leftRewardsTrajectory_int(:,2) = interp1(eventWindow, leftRewardsTrajectory(:,2)', linspace(eventWindow(1),eventWindow(end), 50));
% rightRewardsTrajectory_int(:,1) = interp1(eventWindow, rightRewardsTrajectory(:,1)', linspace(eventWindow(1),eventWindow(end), 50));
% rightRewardsTrajectory_int(:,2) = interp1(eventWindow, rightRewardsTrajectory(:,2)', linspace(eventWindow(1),eventWindow(end), 50));
% leftStimsTrajectory_int(:,1) = interp1(eventWindow, leftStimsTrajectory(:,1)', linspace(eventWindow(1),eventWindow(end), 50));
% leftStimsTrajectory_int(:,2) = interp1(eventWindow, leftStimsTrajectory(:,2)', linspace(eventWindow(1),eventWindow(end), 50));
% rightStimsTrajectory_int(:,1) = interp1(eventWindow, rightStimsTrajectory(:,1)', linspace(eventWindow(1),eventWindow(end), 50));
% rightStimsTrajectory_int(:,2) = interp1(eventWindow, rightStimsTrajectory(:,2)', linspace(eventWindow(1),eventWindow(end), 50));
% 
% figure(4);
% set(gcf,'color','w','Position',[100 100 1025 420]);
% 
% for s = 1:2
%     subplot(1,2,s)
%     plot(0,0)
%     hold on
%     xlim([-.25 .25]);
%     ylim([-.25 .25]);
%     axis square
%     box off
%     ax1 = gca;
%     ax1.TickDir = 'out';
%     ylabel('highL - highR')
%     xlabel('stimL - stimR')
% end
% 
% cd('C:\Users\Wool\Desktop\tempFigs')
% vidObj1 = VideoWriter('movie.avi');
% open(vidObj1);
% loops = 50;
% 
% subplot(1,2,1);
% for f = 1:loops-1
%     
%     plot(leftRewardsTrajectory_int(f:f+1,1),leftRewardsTrajectory_int(f:f+1,2),'Color',colors(1,:),'LineWidth',1);
%     plot(rightRewardsTrajectory_int(f:f+1,1),rightRewardsTrajectory_int(f:f+1,2),'Color',colors(2,:),'LineWidth',1);
%     if f == 17
%         plot(leftRewardsTrajectory_int(f,1),leftRewardsTrajectory_int(f,2),'ko','MarkerFaceColor',[0 .6 .12])
%         plot(rightRewardsTrajectory_int(f,1),rightRewardsTrajectory_int(f,2),'ko','MarkerFaceColor',[0 .6 .12])
%     end
% 
%     hold on
%     F1(f) = getframe(gcf);
%     writeVideo(vidObj1,F1(f));
% end
% 
% subplot(1,2,2)
% for f = 1:loops-1
%     plot(leftStimsTrajectory_int(f:f+1,1),leftStimsTrajectory_int(f:f+1,2),'Color',colors(3,:),'LineWidth',1);
%     plot(rightStimsTrajectory_int(f:f+1,1),rightStimsTrajectory_int(f:f+1,2),'Color',colors(4,:),'LineWidth',1);
%     if f == 17
%         plot(leftStimsTrajectory_int(f,1),leftStimsTrajectory_int(f,2),'ko','MarkerFaceColor',[0 .6 .12])
%         plot(rightStimsTrajectory_int(f,1),rightStimsTrajectory_int(f,2),'ko','MarkerFaceColor',[0 .6 .12])
%     end
% 
%     hold on
%     F1(f) = getframe(gcf);
%     writeVideo(vidObj1,F1(f));
% end




end


%% WORKBENCH
% 
% 
% % pick a 'mother time' to compare all trial-by-trial activity
% %output is a matrix of size trials x cells
% Fs = 0.1;
% motherTime = -.2; %100 ms before stim onset
% eventIdx = find(eventWindow == 0);
% motherIdx = eventIdx + round(motherTime/Fs);
% 
% motherResps = squeeze(alignedResps{2}(:,motherIdx,plotCells));
% 
% % check for cells with nans and remove
% someNaNs = find(isnan(motherResps));
% badCells = unique(ceil(someNaNs/size(motherResps,1)));
% motherResps(:,badCells) = [];
% 
% % compute the mean vectors upon which to project trial-by-trial vectors
% motherResps_reward = nanmean(motherResps(condIdx{strcmp(labels,'allhighLeft')}.train,:),1) - ...
%                         nanmean(motherResps(condIdx{strcmp(labels,'allhighRight')}.train,:),1);
% motherResps_stim = nanmean(motherResps(condIdx{strcmp(labels,'allstimLeft')}.train,:),1) - ...
%                         nanmean(motherResps(condIdx{strcmp(labels,'allstimRight')}.train,:),1);
% 

% plot(leftStimsTrajectory(1,1),leftStimsTrajectory(1,2),'ko','MarkerFaceColor','k')
% plot(leftStimsTrajectory(eventIdx,1),leftStimsTrajectory(eventIdx,2),'ko','MarkerFaceColor',[0 .6 .12])
% plot(leftStimsTrajectory(end,1),leftStimsTrajectory(end,2),'ko','MarkerFaceColor',colorDarkRed)
% plot(rightStimsTrajectory(1,1),rightStimsTrajectory(1,2),'ko','MarkerFaceColor','k')
% plot(rightStimsTrajectory(eventIdx,1),rightStimsTrajectory(eventIdx,2),'ko','MarkerFaceColor',[0 .6 .12])
% plot(rightStimsTrajectory(end,1),rightStimsTrajectory(end,2),'ko','MarkerFaceColor',colorDarkRed)

% plot(leftRewardsTrajectory(1,1),leftRewardsTrajectory(1,2),'ko','MarkerFaceColor','k')
% plot(leftRewardsTrajectory(eventIdx,1),leftRewardsTrajectory(eventIdx,2),'ko','MarkerFaceColor',[0 .6 .12])
% plot(leftRewardsTrajectory(end,1),leftRewardsTrajectory(end,2),'ko','MarkerFaceColor',colorDarkRed)
% plot(rightRewardsTrajectory(1,1),rightRewardsTrajectory(1,2),'ko','MarkerFaceColor','k')
% plot(rightRewardsTrajectory(eventIdx,1),rightRewardsTrajectory(eventIdx,2),'ko','MarkerFaceColor',[0 .6 .12])
% plot(rightRewardsTrajectory(end,1),rightRewardsTrajectory(end,2),'ko','MarkerFaceColor',colorDarkRed)

% plot(leftStimsTrajectory(1,1),leftStimsTrajectory(1,2),'ko','MarkerFaceColor','k')
% plot(leftStimsTrajectory(eventIdx,1),leftStimsTrajectory(eventIdx,2),'ko','MarkerFaceColor',[0 .6 .12])
% plot(leftStimsTrajectory(end,1),leftStimsTrajectory(end,2),'ko','MarkerFaceColor',colorDarkRed)
% plot(rightStimsTrajectory(1,1),rightStimsTrajectory(1,2),'ko','MarkerFaceColor','k')
% plot(rightStimsTrajectory(eventIdx,1),rightStimsTrajectory(eventIdx,2),'ko','MarkerFaceColor',[0 .6 .12])
% plot(rightStimsTrajectory(end,1),rightStimsTrajectory(end,2),'ko','MarkerFaceColor',colorDarkRed)

% plot(leftRewardsTrajectory(1,1),leftRewardsTrajectory(1,2),'ko','MarkerFaceColor','k')
% plot(leftRewardsTrajectory(eventIdx,1),leftRewardsTrajectory(eventIdx,2),'ko','MarkerFaceColor',[0 .6 .12])
% plot(leftRewardsTrajectory(end,1),leftRewardsTrajectory(end,2),'ko','MarkerFaceColor',colorDarkRed)
% plot(rightRewardsTrajectory(1,1),rightRewardsTrajectory(1,2),'ko','MarkerFaceColor','k')
% plot(rightRewardsTrajectory(eventIdx,1),rightRewardsTrajectory(eventIdx,2),'ko','MarkerFaceColor',[0 .6 .12])
% plot(rightRewardsTrajectory(end,1),rightRewardsTrajectory(end,2),'ko','MarkerFaceColor',colorDarkRed)







