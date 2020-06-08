clear all;
% close all;
expInfo = initExpInfo({{'LEW031'}},{{'2020-02-03',1,[1]}});
% expInfo = initExpInfo('LEW031');
matched = 0;
%% load data & extract some variables to make this code still work

if matched == 1
    [expInfo, neuralData, behavioralData] = processExperiment(expInfo,'matched');
    [neuralData] = alignResps(expInfo, neuralData, behavioralData);
    [neuralData] = getSignificantActivity(expInfo, behavioralData, neuralData);
    combinedNeuralData = combineNeuralData(expInfo, behavioralData, neuralData,'matched');

    alignedResps = combinedNeuralData.matched.eta.alignedResps;
    eventWindow = combinedNeuralData.matched.eta.eventWindow;
    bfcH = combinedNeuralData.matched.stats.bfcH;
    pLabels = combinedNeuralData.matched.stats.labels;
    events = combinedNeuralData.matched.eta.events;
elseif matched == 0 
    [expInfo, neuralData, behavioralData] = processExperiment(expInfo);
    [neuralData] = alignResps(expInfo, neuralData, behavioralData);
    [neuralData] = getSignificantActivity(expInfo, behavioralData, neuralData);
    combinedNeuralData = combineNeuralData(expInfo, behavioralData, neuralData);

    alignedResps = neuralData.eta.alignedResps;
    eventWindow = neuralData.eta.eventWindow;
    bfcH = neuralData.stats.bfcH;
    pLabels = neuralData.stats.labels;
    events = neuralData.eta.events;
end

%% get cell responses at a particular timepoint

Fs = 0.1;
%%%%%%%%%%%%%%%% compute baseline activity

% align traces to stim onset
event = 'stimulusOnTimes';
stim_alignedTraces = alignedResps{strcmp(events,event)};
stim_eventWindow = eventWindow;

%designate a baseline window
stim_eventIdx = find(stim_eventWindow == 0);
stim_preTime = [-0.5 0] / Fs;
baselineIdx = stim_eventIdx + stim_preTime(1) : stim_eventIdx;

%compute the mean baseline activity per cell, per trial (trials x neurons)
baselineResps = squeeze(mean(stim_alignedTraces(:,baselineIdx,:),2));

%%%% compute peristimulus activity

%designate a peristimulus window
stimTime = [0 0.4] / Fs;
stimIdx = stim_eventIdx + stimTime(1) :stim_eventIdx + stimTime(2);

%compute the mean peristimulus activity per cell, per trial (trials x neurons)
stimResps = squeeze(mean(stim_alignedTraces(:,stimIdx,:),2));

%%%%%%%%%%%%%%%% compute perimovement activity

% align traces to movement onset
event = 'prestimulusQuiescenceEndTimes';
mov_alignedTraces = alignedResps{strcmp(events,event)};
mov_eventWindow = eventWindow;

%designate a movement window
mov_eventIdx = find(mov_eventWindow == 0);
movTime = [-0.2 0.1] / Fs;
movIdx = mov_eventIdx + movTime(1) : mov_eventIdx + movTime(2);

%compute the mean perimovement activity per cell, per trial (trials x neurons)
movResps = squeeze(mean(mov_alignedTraces(:,movIdx,:),2));

%%%%%%%%%%%%%%%% compute premovement activity

%designate a movement window
pmov_eventIdx = find(mov_eventWindow == 0);
pmovTime = [-0.7 -0.3] / Fs;
pmovIdx = pmov_eventIdx + pmovTime(1) : pmov_eventIdx + pmovTime(2);

%compute the mean perimovement activity per cell, per trial (trials x neurons)
pmovResps = squeeze(mean(mov_alignedTraces(:,pmovIdx,:),2));

%%%%%%%%%%%%%%%% compute perireward activity

% align traces to movement onset
event = 'feedbackTimes';
rew_alignedTraces = alignedResps{strcmp(events,event)};
rew_eventWindow = eventWindow;

%designate a movement window
rew_eventIdx = find(rew_eventWindow == 0);
rewTime = [0 0.2] / Fs;
rewIdx = rew_eventIdx + rewTime(1) : rew_eventIdx + rewTime(2);

%compute the mean perireward activity per cell, per trial (trials x neurons)
rewResps = squeeze(mean(rew_alignedTraces(:,rewIdx,:),2));

%%
resps = stimResps;
whichCells = 'leftStim'; %choose from 'pLabels' array
if strcmp(whichCells, 'all')
    cellActivity = mean(resps, 2);
elseif strcmp(whichCells, 'PCs')
    cellActivity = mean(resps(:,1:5),2);
else
    cellActivity = mean(resps(:,bfcH(:,strcmp(pLabels,whichCells)) > 0), 2);
end
respColor = [0 .4 1];

trials = 1:length(behavioralData.eventTimes(1).daqTime);
windowSize = 10;
windowIdx = [true(1,windowSize), false(1,max(trials)-windowSize)];

figure;hold on
subplot(1,2,1);
hold on;

fakeActivity = randsample(cellActivity,length(cellActivity));
for i = 1:length(cellActivity)
    shiftedActivity = circshift(cellActivity,i);
    shiftedFake = circshift(fakeActivity,i);
    randActivity = randsample(cellActivity,windowSize);
    fMean(i) = mean(shiftedFake(windowIdx));
    fError(i) = std(shiftedFake(windowIdx))/sqrt(windowSize);
    wMean(i) = mean(shiftedActivity(windowIdx));
    wError(i) = std(shiftedActivity(windowIdx))/sqrt(windowSize);
    sMean(i) = mean(randActivity);
    sError(i) = std(randActivity)/sqrt(windowSize);
end


wMean_smoothed = flipud(smooth(wMean,50));
wError_smoothed = flipud(smooth(wError,50));

fMean_smoothed = flipud(smooth(fMean,50));
fError_smoothed = flipud(smooth(fError,50));

sMean_smoothed = smooth(sMean,50);
sError_smoothed = smooth(sError,50);

yMin = 0.999*min([min(wMean_smoothed-wError_smoothed)]);
yMax = 1.001*max([max(wMean_smoothed+wError_smoothed)]);

blockTrials = expInfo.block.events.highRewardSideValues;
blockBounds = [1 find(diff(expInfo.block.events.blockSwitchesValues) ~= 0) max(trials)];
for b = 2:length(blockBounds)
    if mean(blockTrials(blockBounds(b-1):blockBounds(b))) > 0
        color = 'r';
    elseif mean(blockTrials(blockBounds(b-1):blockBounds(b))) < 0
        color = [0 .4 1];
    end
    plotBlock = fill(...
                    [blockBounds(b-1) blockBounds(b) blockBounds(b) blockBounds(b-1)],...
                    [yMin yMin yMax yMax],color,'LineStyle','none');
    alpha(0.2)
end

plot_wMean = plot(wMean_smoothed,'color','k','LineWidth',2);
plot_wError  = fill([trials';fliplr(trials)'],[(wMean_smoothed-wError_smoothed);flipud(wMean_smoothed+wError_smoothed)],'k', 'LineStyle', 'none');
alpha(0.2)

% plot_fMean = plot(fMean_smoothed,'color','k','LineWidth',2);
% plot_fError  = fill([trials';fliplr(trials)'],[(fMean_smoothed-fError_smoothed);flipud(fMean_smoothed+fError_smoothed)],'k', 'LineStyle', 'none');
% alpha(0.2)
   
plot_sMean =  plot(sMean_smoothed,'color',[0.5 0.5 0.5],'LineWidth',2);
plot_sError = fill([trials';fliplr(trials)'],[(sMean_smoothed-sError_smoothed);flipud(sMean_smoothed+sError_smoothed)],[.5 .5 .5], 'LineStyle', 'none');
alpha(0.2)

ax = gca;
ax.TickDir = 'out';
axis([min(trials) max(trials) yMin yMax]);
xlabel('number of trial shifts')
ylabel('normalized mean activity')

subplot(1,2,2);
hold on;

windowIdx = expInfo.block.events.highRewardSideValues == 1 ...
    .* expInfo.block.events.contrastValues < 0;
windowIdx = windowIdx(trials);

for i = 1:length(cellActivity)
    shiftedActivity = circshift(cellActivity,i);
    shiftedFake = circshift(fakeActivity,i);
    randActivity = randsample(cellActivity,length(windowIdx));
    wMean(i) = mean(shiftedActivity(windowIdx));
    wError(i) = std(shiftedActivity(windowIdx))/sqrt(sum(windowIdx));
    sMean(i) = mean(randActivity(windowIdx));
    sError(i) = std(randActivity(windowIdx))/sqrt(sum(windowIdx));
    fMean(i) = mean(shiftedFake(windowIdx));
    fError(i) = std(shiftedFake(windowIdx))/sqrt(windowSize);
end


wMean_smoothed = flipud(smooth(wMean,50));
wError_smoothed = flipud(smooth(wError,50));
sMean_smoothed = smooth(sMean,50);
sError_smoothed = smooth(sError,50);
fMean_smoothed = flipud(smooth(fMean,50));
fError_smoothed = flipud(smooth(fError,50));


yMin = 0.999*min([min(wMean_smoothed-wError_smoothed)]);
yMax = 1.001*max([max(wMean_smoothed+wError_smoothed)]);

blockTrials = expInfo.block.events.highRewardSideValues;
blockBounds = [1 find(diff(expInfo.block.events.blockSwitchesValues) ~= 0) max(trials)];
for b = 2:length(blockBounds)
    if mean(blockTrials(blockBounds(b-1):blockBounds(b))) > 0
        color = 'r';
    elseif mean(blockTrials(blockBounds(b-1):blockBounds(b))) < 0
        color = [0 .4 1];
    end
    plotBlock = fill(...
                    [blockBounds(b-1) blockBounds(b) blockBounds(b) blockBounds(b-1)],...
                    [yMin yMin yMax yMax],color,'LineStyle','none');
    alpha(0.2)
end

plot_wMean = plot(wMean_smoothed,'color','k','LineWidth',2);
plot_wError  = fill([trials';fliplr(trials)'],[(wMean_smoothed-wError_smoothed);flipud(wMean_smoothed+wError_smoothed)],'k', 'LineStyle', 'none');
alpha(0.2)
   
plot_sMean =  plot(sMean_smoothed,'color',[0.5 0.5 0.5],'LineWidth',2);
plot_sError = fill([trials';fliplr(trials)'],[(sMean_smoothed-sError_smoothed);flipud(sMean_smoothed+sError_smoothed)],[.5 .5 .5], 'LineStyle', 'none');
alpha(0.2)

% plot_fMean = plot(fMean_smoothed,'color','k','LineWidth',2);
% plot_fError  = fill([trials';fliplr(trials)'],[(fMean_smoothed-fError_smoothed);flipud(fMean_smoothed+fError_smoothed)],'k', 'LineStyle', 'none');
% alpha(0.2)


ax = gca;
ax.TickDir = 'out';
axis([min(trials) max(trials) yMin yMax]);
xlabel('number of trial shifts')
ylabel('normalized mean activity')
 
%%
% subplot(1,3,3);
% hold on;
% 
% initial_PC_activity = resps(1,1:5);
% for i = 1:length(cellActivity)
%     shiftedActivity = resps(i,1:5);
%     dotCorrelation(i) = dot(initial_PC_activity, shiftedActivity)/(norm(initial_PC_activity)*norm(shiftedActivity));
%     
% end
% 
% blockTrials = expInfo.block.events.highRewardSideValues;
% blockBounds = [1 find(diff(expInfo.block.events.blockSwitchesValues) ~= 0) max(trials)];
% for b = 2:length(blockBounds)
%     if mean(blockTrials(blockBounds(b-1):blockBounds(b))) > 0
%         color = 'r';
%     elseif mean(blockTrials(blockBounds(b-1):blockBounds(b))) < 0
%         color = [0 .4 1];
%     end
%     plotBlock = fill(...
%                     [blockBounds(b-1) blockBounds(b) blockBounds(b) blockBounds(b-1)],...
%                     [0.85 .85 .95 .95],color,'LineStyle','none');
%     alpha(0.2)
% end
% 
% 
% plot_dot = plot((smooth(dotCorrelation)),'color','k','LineWidth',2);
