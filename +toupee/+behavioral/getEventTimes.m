function [expInfo, eventTimes] = getEventTimes(expInfo, names)
% for signalsNames, use the exact fieldnames from the block file or you'll get an error
% e.g. signalsNames = {'stimulusOnTimes' 'interactiveOnTimes' 'stimulusOffTimes'};

% 2018 Nov 21: Now collects feedback times (as well as a rough estimate of when,
% on incorrect trials, the valve would have fired)
% 2019 May 09: Collects feedback times for biased-likelihood exps when
% valve didn't fire
% 2019 Nov 18: expInfo struct input arg
% add multi expInfos

for ex = 1:length(expInfo)
%% PULL OUT EXPINFO VARIABLES

block = expInfo(ex).block;
timeline = expInfo(ex).timeline;

%% DETECT PHD FLIPS
% since the photodiode doesn't always register a WHITE or BLACK signal
% during wheel movement, but also GRAY, we can use departures from GRAY to 
% BLACK or GRAY to WHITE to determine when threshold was reached and the
% feedback signal generated.

phdRaw = timeline.rawDAQData(:,strcmp({timeline.hw.inputs.name},'photoDiode'));

% apply a boxcar filter to get rid of some trace jaggedness
windowSize = 3;
b = (1/windowSize)*ones(1,windowSize);
a = 1;
phdFilt = filter(b, a, phdRaw);

%find the top three clusters of photodiode signal values 
[~, flipvals] = kmeans(phdFilt, 3); %[BLACK WHITE GRAY]

%interpolate the photodiode signal to one of these three cluster values
phdInterp = interp1(flipvals, flipvals, phdFilt, 'nearest', 'extrap');

%note when the phd flips occurred by capturing the diffs
whenFlips = abs(diff(phdInterp)) > 0;

% get rid of the timestamps when there was no flip
flipTimes = timeline.rawDAQTimestamps(2:end); %good enough for govt work
flipTimes(whenFlips < 1) = [];

%% MATCH UP BLOCK TIMES TO STIM ONSET TIMES
% determine the closest value of stimWindowUpdateTimes to
% each relevant block event, and then map this to a nearby flip of the
% photodiode

numCompleteTrials = numel(block.events.endTrialValues);

% find the closest stimWindowUpdate timepoint to every block event
for iName = 1:length(names)
    signalsUpdates{iName} = interp1(block.stimWindowUpdateTimes, block.stimWindowUpdateTimes, block.events.(names{iName}), 'nearest', 'extrap');
end

% set up the corrected-timestamp struct
eventTimes = struct('event',names);

for i = 1:length(signalsUpdates)
    updateTimes = signalsUpdates{i}(1:numCompleteTrials);
    
    for j = 1:numel(updateTimes)

        % choose a timepoint in signals
        timepoint = updateTimes(j);

        % find the earliest phd flip after the relevant stimWindowUpdate
        % timepoint
        firstFlip = first(flipTimes(flipTimes > timepoint));
        
        % collect these values for the minute
        realTimes(j) = firstFlip;

    end
    
    % save signals event times wrt to stimWindowUpdate times...
    eventTimes(i).signalsTime = getfield(block.events,names{i});
    eventTimes(i).signalsTime = eventTimes(i).signalsTime(1:numCompleteTrials);
    eventTimes(i).updateTime = updateTimes;
    
    % ... and also phd flip times
    eventTimes(i).daqTime = realTimes;
    
    % tidy up
    updateTimes = [];
    realTimes = [];
end


%% REWARD TIMINGS
% reward timings are also recorded on the Daq, we can take these at face
% value and add them to 'timelineTimes'

%gather the reward pulse from timeline
reward = timeline.rawDAQData(:,strcmp({timeline.hw.inputs.name},'rewardCommand'));

rewardThresh = max(reward)/2;
rewardTrace = reward > rewardThresh;
realRewardTimes = timeline.rawDAQTimestamps(find(rewardTrace(2:end) & ~rewardTrace(1:end-1))+1);

% isolate keyboard rewards by first interpolating them to the Daq times
allRewardTimes = realRewardTimes;
keyboardRewards = interp1(realRewardTimes,realRewardTimes,block.inputs.keyboardTimes,'nearest');

% then identify them within the 'all rewards' vector, and remove them
[~,kIdx] = intersect(realRewardTimes,keyboardRewards);
realRewardTimes(kIdx) = [];

feedbackRewards = interp1(block.outputs.rewardTimes, block.outputs.rewardTimes, block.events.feedbackTimes, 'nearest', 'extrap');
feedbackRewards = block.events.feedbackValues.*feedbackRewards;
feedbackRewards(feedbackRewards == 0) = NaN;

[~,rewardTrialIdx]=find(block.events.feedbackValues == 1);

% --------- if a likelihood experiment
allBlockRewardValues = block.outputs.rewardValues;
blockKeyboardRewards = interp1(block.outputs.rewardTimes,block.outputs.rewardTimes,block.inputs.keyboardTimes,'nearest');
[~,bkIdx] = intersect(block.outputs.rewardTimes, blockKeyboardRewards);
allBlockRewardValues(bkIdx) = [];
try 
    isProbExp = block.paramsValues(1).rewardProbability(1) > 0;
    if exist('isProbExp')
        rewardTrialIdx(allBlockRewardValues == 0) = []; %if biased likelihood experiment, take out 'rewards' of 0
        rewardTrialIdx(find(rewardTrialIdx > numCompleteTrials)) = [];
    end
catch
end
% ---------

rewardTrialTimes = nan(1,numCompleteTrials);
rewardTrialTimes(rewardTrialIdx) = realRewardTimes;

% record these values
l = length(names);
eventTimes(l+1).event = 'rewardOnTimes';
eventTimes(l+1).signalsTime = feedbackRewards(1:numCompleteTrials);
eventTimes(l+1).daqTime = rewardTrialTimes(1:numCompleteTrials);

% compute average valve delay
meanRewDiff = nanmean(eventTimes(l+1).daqTime - eventTimes(l+1).signalsTime);

% record the feedback timings so you can compare timepoints between
% rewarded and unrewarded trials
eventTimes(l+2).event = 'feedbackTimes';
eventTimes(l+2).signalsTime = block.events.feedbackTimes(1:numCompleteTrials);
%a rough guess at when the mouse would expect the valve to fire, based on
%the mean valve delay
eventTimes(l+2).daqTime = interp1(timeline.rawDAQTimestamps,timeline.rawDAQTimestamps,block.events.feedbackTimes(1:numCompleteTrials) + meanRewDiff,'nearest');

eventTimes{ex} = eventTimes;
end


