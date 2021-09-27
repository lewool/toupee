function allEventTimes = getEventTimes(expInfo, signalsNames)
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
Timeline = expInfo(ex).Timeline;

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

% external TL now needs an offset to bring it close to signals times.
% do this by aligning the first 'real' flip to the first stim presentation
% plus a little buffer so TL is not brought ahead
tlOffset = flipTimes(3) - signalsUpdates{1}(1) - 0.1;

% set up the corrected-timestamp struct
eventTimes = struct('event',signalsNames);

for i = 1:length(signalsUpdates)
    updateTimes = signalsUpdates{i}(1:numCompleteTrials);
    updateTimes_offset = signalsUpdates{i}(1:numCompleteTrials) + tlOffset;
    
    for j = 1:numel(updateTimes)

        % choose a timepoint in signals
        timepoint = updateTimes_offset(j);

        % find the earliest phd flip after the relevant stimWindowUpdate
        % timepoint
        firstFlip = first(flipTimes(flipTimes > timepoint));
        
         % collect these values for the minute
        try
            realTimes(j) = firstFlip;
        catch
            disp(j)
            error('oops')
        end

    end
    
    % save signals event times wrt to stimWindowUpdate times...
    eventTimes(i).signalsTime = getfield(block.events,signalsNames{i});
    eventTimes(i).signalsTime = eventTimes(i).signalsTime(1:numCompleteTrials);
    eventTimes(i).updateTime = updateTimes;
    
    % ... and also phd flip times
    eventTimes(i).daqTime = realTimes;
    
    % tidy up
    updateTimes_offset = [];
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
% allRewardTimes = realRewardTimes;
% keyboardRewards = interp1(realRewardTimes,realRewardTimes,block.inputs.keyboardTimes+tlOffset,'nearest');
% 
% % then identify them within the 'all rewards' vector, and remove them
% [~,kIdx] = intersect(realRewardTimes,keyboardRewards);
% realRewardTimes(kIdx) = [];

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
    isProbExp = block.paramsValues(1).highProbability(1) > 0;
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
l = length(signalsNames);
eventTimes(l+1).event = 'rewardOnTimes';
eventTimes(l+1).signalsTime = feedbackRewards(1:numCompleteTrials);
eventTimes(l+1).daqTime = rewardTrialTimes(1:numCompleteTrials);

% compute average valve delay
meanRewDiff = nanmean(eventTimes(l+1).daqTime - eventTimes(l+1).signalsTime - tlOffset);

% record the feedback timings so you can compare timepoints between
% rewarded and unrewarded trials
eventTimes(l+2).event = 'feedbackTimes';
eventTimes(l+2).signalsTime = block.events.feedbackTimes(1:numCompleteTrials);
%a rough guess at when the mouse would expect the valve to fire, based on
%the mean valve delay
eventTimes(l+2).daqTime = interp1(Timeline.rawDAQTimestamps,Timeline.rawDAQTimestamps,block.events.feedbackTimes(1:numCompleteTrials)+ tlOffset + meanRewDiff,'nearest');
 

% %% RAW WHEEL TRACES
% 
% rawPos = block.inputs.wheelValues;
% rawPos = wheel.correctCounterDiscont(rawPos); % correction because sometimes negative wheel positions wrap around
% rawPos = rawPos(2:end); %the first timestamp seems to be messed up, so ignore it
% rawTimes = block.inputs.wheelTimes(2:end);
% 
% % sample/interpolate the trace
% Fs = 1000;
% t = rawTimes(1):1/Fs:rawTimes(end);
% t = interp1(Timeline.rawDAQTimestamps, Timeline.rawDAQTimestamps, t, 'nearest', 'extrap');
% pos = interp1(rawTimes, rawPos, t, 'linear');
% 
% wheelRadius = 31; % mm (burgess wheel) (measured by Chris)
% rotaryEncoderResolution = 100*4; % number of ticks for one revolution (factor of 4 is according to CB)
% wheelGain = block.paramsValues(1).wheelGain;
% 
% pos = pos./(rotaryEncoderResolution)*2*pi*wheelRadius; % convert to mm
% rawPos = rawPos./(rotaryEncoderResolution)*2*pi*wheelRadius; % convert to mm
% 
% deg = pos.*wheelGain;
% rawDeg = rawPos.*wheelGain;
% 
% %% WHEEL TRAJECTORIES
% 
% % block events
% newTrialTimes = interp1(t, t, block.events.newTrialTimes, 'nearest', 'extrap');
% interactiveOnTimes = interp1(t, t, block.events.interactiveOnTimes, 'nearest', 'extrap');
% 
% for i = 1:numCompleteTrials
%     
%     %consider the start of the trace when block says the trial starts. this is a bit
%     %arbitrary but there's not much going on at this part of any trial so a
%     %good place to make an arbitrary cut
%     idxStart = find(t == newTrialTimes(i)); 
%     
%     %and so the end of the trace should be the last timepoint before the next trial
%     try idxEnd = find(t == newTrialTimes(i+1))-1; 
%     catch
%         %or, it's the end of the whole trace (re: the last trial)
%         idxEnd = numel(t); 
%     end
%     
%     %find the timepoint at which the interactive period started
%     interIdx = find(t == interactiveOnTimes(i));
%     
%     %collect the timepoints for a single trial's wheel trace
%     trialIdx = idxStart:idxEnd;
%     trialWheelTimes = t(trialIdx);
%     
%     % record the raw times in their own struct
%     wheelTrajectories{ex}(i).signalsTimes = trialWheelTimes;
%     wheelTrajectories{ex}(i).values = deg(trialIdx);
%     
%     % this offset value lets us zero the wheel position at the
%     % start of the interactive period
%     wheelTrajectories{ex}(i).interactivePosOffset = deg(interIdx);
% end
% 
% %% DETERMINE MOVEMENT PERIODS
% 
% %compute wheel velocity
% [vel, ~] = wheel.computeVelocity(pos, 100, Fs);
% 
% %threshold the velocity to determine when movements started and stopped
% vel(isnan(vel)) = 0; %weird
% env = envelope(vel,100);
% isMoving = env > 50; %?
% 
% %look back in time a bit to where velocity actually started changing
% timeAround = 0.075;
% intAround = timeAround * Fs;
% 
% idxIsGoing = find(diff(isMoving)>0);
% for tu = 1:length(idxIsGoing)
%     if idxIsGoing(tu) < intAround
%         [~, gIdx] = min(env(1:idxIsGoing(tu)));
%         timesStartedMovement(tu) = t(gIdx);
%     else
%         [~, gIdx] = min(env(idxIsGoing(tu)-intAround:idxIsGoing(tu)));
%         timesStartedMovement(tu) = t(idxIsGoing(tu)-intAround+gIdx);
%     end
% end
% 
% idxIsStopping = find(diff(isMoving)<0);
% for td = 1:length(idxIsStopping)
%     if idxIsStopping(td) + intAround > length(env)
%         [~, sIdx] = min(env(idxIsStopping(td):end));
%         timesEndedMovement(td) = t(idxIsStopping(td)+sIdx-1);
%     else
%     	[~, sIdx] = min(env(idxIsStopping(td):idxIsStopping(td)+intAround));
%         timesEndedMovement(td) = t(idxIsStopping(td)+sIdx-1);
%     end
%     
% end
% 
% for iTime = 1:numCompleteTrials
%     % find when quiescence ends (usually after stimulus onset)
%     timepoint = block.events.stimulusOnTimes(iTime);
%     outcomePoint = block.events.feedbackTimes(iTime);
%     try 
%         firstMove = timesStartedMovement(find(timesStartedMovement > (timepoint-.5)));
%         firstMove = firstMove(firstMove -eventTimes(1).daqTime(iTime) > 0);
%         firstMove = firstMove(1);
%         prestimulusQuiescenceEnded(iTime) = firstMove;
%     catch
%         prestimulusQuiescenceEnded(iTime) = NaN; 
%     end
%     
%     
%     
%     % but flag any movements that occur 100ms around stimulus onset (these
%     % can't be real responses)
%     if abs(firstMove - timepoint) < 0.100
%         goodQuiescenceFlag(iTime) = false;
%         wheelTrajectories{ex}(iTime).goodQuiescenceFlag = false;
%     else
%         goodQuiescenceFlag(iTime) = true;
%         wheelTrajectories{ex}(iTime).goodQuiescenceFlag = true;
%     end
%     
%     % find when quiescence starts (before stimulus onset)
%     try
%         lastMove = timesEndedMovement(find(timesEndedMovement < (timepoint-.5)));
%         lastMove = lastMove(end);
%         prestimulusQuiescenceStarted(iTime) = lastMove;
%     catch
%         prestimulusQuiescenceStarted(iTime) = NaN;
%     end
%     
% end

% %% MOVEMENT TIMINGS
% 
% l = length(eventTimes);
% eventTimes(l+1).event = 'prestimulusQuiescenceStartTimes';
% eventTimes(l+1).daqTime = prestimulusQuiescenceStarted;
% eventTimes(l+2).event = 'prestimulusQuiescenceEndTimes';
% eventTimes(l+2).daqTime = prestimulusQuiescenceEnded;
% 
allEventTimes{ex} = eventTimes;
end


