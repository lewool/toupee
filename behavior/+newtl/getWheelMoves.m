function allWheelMoves = getWheelMoves(expInfo, allEventTimes)

for ex = 1:length(expInfo)

    %% load data

    block = expInfo(ex).block;
    Timeline = expInfo(ex).Timeline;
    eventTimes = allEventTimes{ex};

    numCompleteTrials = numel(block.events.endTrialValues);

    %% feed in the wheel trace from the experiment
    rawPos = Timeline.rawDAQData(:,strcmp({Timeline.hw.inputs.name},'rotaryEncoder'))';
    rawPos = wheel.correctCounterDiscont(rawPos); % correction because sometimes negative wheel positions wrap around
    rawTimes = Timeline.rawDAQTimestamps;

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

    % %compute wheel velocity
    [vel, ~] = wheel.computeVelocity2(pos, 1, Fs);

    %% trial-by-trial traces

    % event times
    tlOffset = eventTimes(1).daqTime(1) - eventTimes(1).signalsTime(1);
    newTrialTimes = interp1(t, t, block.events.newTrialTimes + tlOffset, 'nearest', 'extrap');
    prestimTimes = interp1(t, t, eventTimes(1).daqTime - 0.5, 'nearest','extrap');
    stimTimes = interp1(t, t, eventTimes(1).daqTime, 'nearest','extrap');

    for i = 1:numCompleteTrials

        %consider the start of the trace when block says the trial starts. this
        %doesn't need to be precise
        idxStart = find(t == newTrialTimes(i));

        %identify a prestim index
        idxPre = find(t > prestimTimes(i) & t <= stimTimes(i));

        %and so the end of the trace should be the last timepoint before the next trial
        try idxEnd = find(t == newTrialTimes(i+1))-1; 
        catch
            %or, it's the end of the whole trace (re: the last trial)
            idxEnd = numel(t); 
        end

        %collect the timepoints for a single trial's wheel trace
        trialIdx = idxStart:idxEnd;
        trialWheelTimes = t(trialIdx);

        % record the traces/times
        wheelTraceTimes{i} = trialWheelTimes;
        wheelTraceValues{i} = pos(trialIdx) - mean(pos(idxPre));

    end

    %% find wheel moves

    [inMove, moveOnsets, moveOffsets, moveAmps, peakVelTimes] = wheel.findWheelMoves3(pos, t, Fs);
    peakVelTimes = interp1(t,t,peakVelTimes,'nearest','extrap');

    %% record movement in each epoch of each trial 
    % three epochs: prestimulus period, immediately post-stimulus, or between
    % go cue and feedback
    
    isMoving_prestim = zeros(1,numCompleteTrials);
    isMoving_early = zeros(1,numCompleteTrials);
    isMoving_feedback = zeros(1,numCompleteTrials);
    
    moveOnset_prestim = zeros(1,numCompleteTrials);
    moveOffset_prestim = zeros(1,numCompleteTrials);
    moveDir_prestim = zeros(1,numCompleteTrials);
    peakVel_prestim = zeros(1,numCompleteTrials);
    
    moveOnset_early = zeros(1,numCompleteTrials);
    moveOffset_early = zeros(1,numCompleteTrials);
    moveDir_early = zeros(1,numCompleteTrials);
    peakVel_early = zeros(1,numCompleteTrials);
    
    moveOnset_goCue = zeros(1,numCompleteTrials);
    moveOffset_goCue = zeros(1,numCompleteTrials);
    moveDir_goCue = zeros(1,numCompleteTrials);
    peakVel_goCue = zeros(1,numCompleteTrials);
    
    moveOnset_feedback = zeros(1,numCompleteTrials);
    moveOffset_feedback = zeros(1,numCompleteTrials);
    moveDir_feedback = zeros(1,numCompleteTrials);
    peakVel_feedback = zeros(1,numCompleteTrials);

    for iTrial = 1:numCompleteTrials

        %find bounding times
        prestimTime = eventTimes(1).daqTime(iTrial) - 0.5;
        stimOnTime = eventTimes(1).daqTime(iTrial);
        goCueTime = eventTimes(2).daqTime(iTrial);
        feedbackTime = eventTimes(5).daqTime(iTrial);

        %is the mouse moving during any these epochs?
        isMoving_prestim(iTrial) = any(inMove(t > prestimTime & t <= stimOnTime));
        isMoving_early(iTrial) = any(inMove(t > stimOnTime & t <= goCueTime));
        isMoving_feedback(iTrial) = any(inMove(t > goCueTime & t <= feedbackTime));

        %if there is movement in the 'prestim' epoch...
        if isMoving_prestim(iTrial)

            %characterize the last prestim movement, if an onset exists
            moveIdx_prestim = last(find(moveOnsets > prestimTime & moveOnsets <= stimOnTime));
            if ~isempty(moveIdx_prestim)
                moveOnset_prestim(iTrial) = moveOnsets(moveIdx_prestim);
                moveOffset_prestim(iTrial) = moveOffsets(moveIdx_prestim);
                moveDir_prestim(iTrial) = -sign(moveAmps(moveIdx_prestim));
                peakVel_prestim(iTrial) = vel(t == peakVelTimes(moveIdx_prestim));
            else
                moveOnset_prestim(iTrial) = NaN;
                moveOffset_prestim(iTrial) = NaN;
                moveDir_prestim(iTrial) = NaN;
                peakVel_prestim(iTrial) = NaN;
            end
        else
            moveOnset_prestim(iTrial) = NaN;
            moveOffset_prestim(iTrial) = NaN;
            moveDir_prestim(iTrial) = NaN;
            peakVel_prestim(iTrial) = NaN;
        end

        %if there is movement in the 'early' epoch...
        if isMoving_early(iTrial)

            %characterize the first early movement, if an onset exists
            moveIdx_early = first(find(moveOnsets > stimOnTime & moveOnsets <= goCueTime));
            if ~isempty(moveIdx_early)
                moveOnset_early(iTrial) = moveOnsets(moveIdx_early);
                moveOffset_early(iTrial) = moveOffsets(moveIdx_early);
                moveDir_early(iTrial) = -sign(moveAmps(moveIdx_early));
                peakVel_early(iTrial) = vel(t == peakVelTimes(moveIdx_early));
            else
                moveOnset_early(iTrial) = NaN;
                moveOffset_early(iTrial) = NaN;
                moveDir_early(iTrial) = NaN;
                peakVel_early(iTrial) = NaN;
            end
        else
            moveOnset_early(iTrial) = NaN;
            moveOffset_early(iTrial) = NaN;
            moveDir_early(iTrial) = NaN;
            peakVel_early(iTrial) = NaN;
        end

        %if there is movement in the 'feedback' epoch...
        if isMoving_feedback(iTrial)

            %characterize the first movement after the go cue, if an onset
            %exists
            moveIdx_goCue = first(find(moveOnsets > goCueTime & moveOnsets <= feedbackTime));
            if ~isempty(moveIdx_goCue)
                moveOnset_goCue(iTrial) = moveOnsets(moveIdx_goCue);
                moveOffset_goCue(iTrial) = moveOffsets(moveIdx_goCue);
                moveDir_goCue(iTrial) = -sign(moveAmps(moveIdx_goCue));
                peakVel_goCue(iTrial) = vel(t == peakVelTimes(moveIdx_goCue));
            else
                moveOnset_goCue(iTrial) = NaN;
                moveOffset_goCue(iTrial)= NaN;
                moveDir_goCue(iTrial) = NaN;
                peakVel_goCue(iTrial) = NaN;
            end

            %characterize the movement that triggered feedback, if an onset
            %exists
            moveIdx_feedback = last(find(moveOnsets > goCueTime & moveOnsets <= feedbackTime));
            if ~isempty(moveIdx_feedback)
                moveOnset_feedback(iTrial) = moveOnsets(moveIdx_feedback);
                moveOffset_feedback(iTrial) = moveOffsets(moveIdx_feedback);
                moveDir_feedback(iTrial) = -sign(moveAmps(moveIdx_feedback));
                peakVel_feedback(iTrial) = vel(t == peakVelTimes(moveIdx_feedback));
            else
                moveOnset_feedback(iTrial) = NaN;
                moveOffset_feedback(iTrial) = NaN;
                moveDir_feedback(iTrial) = NaN;
                peakVel_feedback(iTrial) = NaN;
            end

        else
            moveOnset_goCue(iTrial) = NaN;
            moveOffset_goCue(iTrial)= NaN;
            moveDir_goCue(iTrial) = NaN;
            peakVel_goCue(iTrial) = NaN;
            moveOnset_feedback(iTrial) = NaN;
            moveOffset_feedback(iTrial) = NaN;
            moveDir_feedback(iTrial) = NaN;
            peakVel_feedback(iTrial) = NaN;
        end

    end
    
    %identify the first movement of the trial, redundant info but easy to
    %do
    firstMoveOnset = nan(1,numCompleteTrials);
    firstMoveOnset((~isnan(moveOnset_goCue))) = moveOnset_goCue(~isnan(moveOnset_goCue));
    firstMoveOnset((~isnan(moveOnset_early))) = moveOnset_early(~isnan(moveOnset_early));
    
    firstMoveOffset = nan(1,numCompleteTrials);
    firstMoveOffset((~isnan(moveOffset_goCue))) = moveOffset_goCue(~isnan(moveOffset_goCue));
    firstMoveOffset((~isnan(moveOffset_early))) = moveOffset_early(~isnan(moveOffset_early));
    
    firstMoveDir = nan(1,numCompleteTrials);
    firstMoveDir((~isnan(moveDir_goCue))) = moveDir_goCue(~isnan(moveDir_goCue));
    firstMoveDir((~isnan(moveDir_early))) = moveDir_early(~isnan(moveDir_early));
    
    firstMoveVel = nan(1,numCompleteTrials);
    firstMoveVel((~isnan(peakVel_goCue))) = peakVel_goCue(~isnan(peakVel_goCue));
    firstMoveVel((~isnan(peakVel_early))) = peakVel_early(~isnan(peakVel_early));

    
    %% load into struct

    wheelMoves.epochs(1) = struct(...
        'epoch', 'prestimulus',...
        'isMoving', isMoving_prestim,... 
        'onsetTimes', moveOnset_prestim,...
        'offsetTimes', moveOffset_prestim,...
        'moveDir', moveDir_prestim,...
        'peakVel', peakVel_prestim ...
        );

    wheelMoves.epochs(2) = struct(...
        'epoch', 'early',...
        'isMoving', isMoving_early,... 
        'onsetTimes', moveOnset_early,...
        'offsetTimes', moveOffset_early,...
        'moveDir', moveDir_early,...
        'peakVel', peakVel_early ...
        );

    wheelMoves.epochs(3) = struct(...
        'epoch', 'goCue',...
        'isMoving', isMoving_feedback,... 
        'onsetTimes', moveOnset_goCue,...
        'offsetTimes', moveOffset_goCue,...
        'moveDir', moveDir_goCue,...
        'peakVel', peakVel_goCue ...
        );

    wheelMoves.epochs(4) = struct(...
        'epoch', 'feedback',...
        'isMoving', isMoving_feedback,... 
        'onsetTimes', moveOnset_feedback,...
        'offsetTimes', moveOffset_feedback,...
        'moveDir', moveDir_feedback,...
        'peakVel', peakVel_feedback ...
        );
    
    wheelMoves.epochs(5) = struct(...
        'epoch', 'firstMove',...
        'isMoving', [],... 
        'onsetTimes', firstMoveOnset,...
        'offsetTimes', firstMoveOffset,...
        'moveDir', firstMoveDir,...
        'peakVel', firstMoveVel ...
        );

     wheelMoves.traces.pos = wheelTraceValues;
     wheelMoves.traces.time = wheelTraceTimes;
     
     allWheelMoves{ex} = wheelMoves;
     clearvars -except ex expInfo allWheelMoves allEventTimes
     
     
end