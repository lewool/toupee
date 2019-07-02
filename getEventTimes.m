function [eventTimes, wheelTrajectories] = getEventTimes(block, Timeline, signalsNames)
% for signalsNames, use the exact fieldnames from the block file or you'll get an error
% e.g. signalsNames = {'stimulusOnTimes' 'interactiveOnTimes' 'stimulusOffTimes'};

%% DETECT PHD FLIPS
% since the photodiode doesn't always register a WHITE or BLACK signal
% during wheel movement, but also GRAY, we can use departures from GRAY to 
% BLACK or GRAY to WHITE to determine when threshold was reached and the
% feedback signal generated.

phdRaw = Timeline.rawDAQData(:,strcmp({Timeline.hw.inputs.name},'photoDiode'));

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
flipTimes = Timeline.rawDAQTimestamps(2:end); %good enough for govt work
flipTimes(whenFlips < 1) = [];

%% MATCH UP BLOCK TIMES TO STIM ONSET TIMES
% determine the closest value of stimWindowUpdateTimes to
% each relevant block event, and then map this to a nearby flip of the
% photodiode

numCompleteTrials = numel(block.events.endTrialValues);

% find the closest stimWindowUpdate timepoint to every block event
for iName = 1:length(signalsNames)
    signalsUpdates{iName} = interp1(block.stimWindowUpdateTimes, block.stimWindowUpdateTimes, block.events.(signalsNames{iName}), 'nearest', 'extrap');
end

% set up the corrected-timestamp struct
eventTimes = struct('event',signalsNames);

for i = 1:length(signalsUpdates)
    updateTimes = signalsUpdates{i};
    
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
    eventTimes(i).signalsTime = getfield(block.events,signalsNames{i});
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

%gather the reward pulse from Timeline
reward = Timeline.rawDAQData(:,strcmp({Timeline.hw.inputs.name},'rewardCommand'));

rewardThresh = max(reward)/2;
rewardTrace = reward > rewardThresh;
realRewardTimes = Timeline.rawDAQTimestamps(find(rewardTrace(2:end) & ~rewardTrace(1:end-1))+1);

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
% rewardTrialIdx(find(rewardTrialIdx > numCompleteTrials)) = [];
rewardTrialTimes = nan(1,numCompleteTrials);
rewardTrialTimes(rewardTrialIdx) = realRewardTimes;

% record these values
l = length(signalsNames);
eventTimes(l+1).event = 'rewardOnTimes';
eventTimes(l+1).signalsTime = feedbackRewards(1:numCompleteTrials);
eventTimes(l+1).daqTime = rewardTrialTimes(1:numCompleteTrials);

%% RAW WHEEL TRACES

rawPos = block.inputs.wheelValues;
rawPos = wheel.correctCounterDiscont(rawPos); % correction because sometimes negative wheel positions wrap around
rawTimes = block.inputs.wheelTimes;

% sample/interpolate the trace
Fs = 1000;
t = rawTimes(1):1/Fs:rawTimes(end);
t = interp1(Timeline.rawDAQTimestamps, Timeline.rawDAQTimestamps, t, 'nearest', 'extrap');
pos = interp1(rawTimes, rawPos, t, 'linear');

wheelRadius = 31; % mm (burgess wheel) (measured by Chris)
rotaryEncoderResolution = 100*4; % number of ticks for one revolution (factor of 4 is according to CB)
wheelGain = block.paramsValues(1).wheelGain;

pos = pos./(rotaryEncoderResolution)*2*pi*wheelRadius; % convert to mm
rawPos = rawPos./(rotaryEncoderResolution)*2*pi*wheelRadius; % convert to mm

deg = pos.*wheelGain;
rawDeg = rawPos.*wheelGain;

%compute wheel velocity
[vel, ~] = wheel.computeVelocity(pos, 100, Fs);

%% WHEEL TRAJECTORIES

% block events
newTrialTimes = interp1(t, t, block.events.newTrialTimes, 'nearest', 'extrap');
interactiveOnTimes = interp1(t, t, block.events.interactiveOnTimes, 'nearest', 'extrap');

for i = 1:numCompleteTrials
    
    %consider the start of the trace when block says the trial starts. this is a bit
    %arbitrary but there's not much going on at this part of any trial so a
    %good place to make an arbitrary cut
    idxStart = find(t == newTrialTimes(i)); 
    
    %and so the end of the trace should be the last timepoint before the next trial
    try idxEnd = find(t == newTrialTimes(i+1))-1; 
    catch
        %or, it's the end of the whole trace (re: the last trial)
        idxEnd = numel(t); 
    end
    
    %find the timepoint at which the interactive period started
    interIdx = find(t == interactiveOnTimes(i));
    
    %collect the timepoints for a single trial's wheel trace
    trialIdx = idxStart:idxEnd;
    trialWheelTimes = t(trialIdx);
    
    % record the raw times in their own struct
    wheelTrajectories(i).signalsTimes = trialWheelTimes;
    wheelTrajectories(i).pos = deg(trialIdx);
    wheelTrajectories(i).vel = vel(trialIdx);
    
    % this offset value lets us zero the wheel position at the
    % start of the interactive period
    wheelTrajectories(i).interactivePosOffset = deg(interIdx);
end

%% DETERMINE MOVEMENT PERIODS

%threshold the velocity to determine when movements started and stopped
vel(isnan(vel)) = 0; %weird
env = envelope(vel,100);
isMoving = env > block.paramsValues(1).quiescenceThreshold * 5; %?
timesStartedMovement = t(2:end) .* (diff(isMoving) > 0);
timesStartedMovement(timesStartedMovement == 0) = [];
timesEndedMovement = t(2:end) .* (diff(isMoving) < 0);
timesEndedMovement(timesEndedMovement == 0) = [];

for iTime = 1:numCompleteTrials
    % find when quiescence ends (usually after stimulus onset)
    timepoint = block.events.stimulusOnTimes(iTime);
    try 
        firstMove = timesStartedMovement(find(timesStartedMovement > (timepoint-.5)));
        firstMove = firstMove(1);
        prestimulusQuiescenceEnded(iTime) = firstMove;
    catch
        prestimulusQuiescenceEnded(iTime) = NaN; 
    end
    
    
    
    % but flag any movements that occur 100ms around stimulus onset (these
    % can't be real responses)
    if abs(firstMove - timepoint) < 0.100
        goodQuiescenceFlag(iTime) = false;
        wheelTrajectories(iTime).goodQuiescenceFlag = false;
    else
        goodQuiescenceFlag(iTime) = true;
        wheelTrajectories(iTime).goodQuiescenceFlag = true;
    end
    
    % find when quiescence starts (before stimulus onset)
    try
        lastMove = timesEndedMovement(find(timesEndedMovement < (timepoint-.5)));
        lastMove = lastMove(end);
        prestimulusQuiescenceStarted(iTime) = lastMove;
    catch
        prestimulusQuiescenceStarted(iTime) = NaN;
    end
    
end

%% MOVEMENT TIMINGS

l = length(eventTimes);
eventTimes(l+1).event = 'prestimulusQuiescenceStartTimes';
eventTimes(l+1).daqTime = prestimulusQuiescenceStarted;
eventTimes(l+2).event = 'prestimulusQuiescenceEndTimes';
eventTimes(l+2).daqTime = prestimulusQuiescenceEnded;

%% MAX VELOCITY PER TRIAL

for iTrial = 1:length(eventTimes(strcmp({eventTimes.event},'feedbackTimes')).daqTime)
    startMove = interp1(wheelTrajectories(iTrial).signalsTimes, wheelTrajectories(iTrial).signalsTimes, eventTimes(strcmp({eventTimes.event},'prestimulusQuiescenceEndTimes')).daqTime(iTrial), 'nearest', 'extrap');
    endMove = interp1(wheelTrajectories(iTrial).signalsTimes, wheelTrajectories(iTrial).signalsTimes, eventTimes(strcmp({eventTimes.event},'feedbackTimes')).daqTime(iTrial), 'nearest', 'extrap');
    mvIdx = find(endMove == wheelTrajectories(iTrial).signalsTimes);
    try
        [~,maxIdx] = find(abs(wheelTrajectories(iTrial).vel) == max(abs(wheelTrajectories(iTrial).vel(mvIdx-50:mvIdx))));
        wheelTrajectories(iTrial).maxVel = wheelTrajectories(iTrial).vel(maxIdx);
    catch
        wheelTrajectories(iTrial).maxVel = NaN;
    end
end



end


