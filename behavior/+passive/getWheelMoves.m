function allWheelMoves = getWheelMoves(expInfo, allEventTimes)

for ex = 1:length(expInfo)
    %% load data

    block = expInfo(ex).block;
    Timeline = expInfo(ex).Timeline;
    eventTimes = allEventTimes{ex};
    
    numCompleteTrials = numel(block.events.endTrialValues);
    
    if isstruct(eventTimes) && numCompleteTrials > 0
        
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
        wheelGain = 1;

        pos = pos./(rotaryEncoderResolution)*2*pi*wheelRadius; % convert to mm
        rawPos = rawPos./(rotaryEncoderResolution)*2*pi*wheelRadius; % convert to mm

        % %compute wheel velocity
        [vel, ~] = wheel.computeVelocity2(pos, 1, Fs);

        %% trial-by-trial traces

        % event times
        newTrialTimes = interp1(t, t, block.events.newTrialTimes, 'nearest', 'extrap');
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
        isMoving_stim = zeros(1,numCompleteTrials);

        moveOnset_prestim = zeros(1,numCompleteTrials);
        moveOffset_prestim = zeros(1,numCompleteTrials);
        moveDir_prestim = zeros(1,numCompleteTrials);
        peakVel_prestim = zeros(1,numCompleteTrials);

        moveOnset_stim = zeros(1,numCompleteTrials);
        moveOffset_stim = zeros(1,numCompleteTrials);
        moveDir_stim = zeros(1,numCompleteTrials);
        peakVel_stim = zeros(1,numCompleteTrials);

        for iTrial = 1:numCompleteTrials

            %find bounding times
            prestimTime = eventTimes(1).daqTime(iTrial) - 0.5;
            stimOnTime = eventTimes(1).daqTime(iTrial);
            stimOffTime = eventTimes(2).daqTime(iTrial);

            %is the mouse moving during any these epochs?
            isMoving_prestim(iTrial) = any(inMove(t > prestimTime & t <= stimOnTime));
            isMoving_stim(iTrial) = any(inMove(t > stimOnTime & t <= stimOffTime));

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

            %if there is movement in the 'stim' epoch...
            if isMoving_stim(iTrial)

                %characterize the first early movement, if an onset exists
                moveIdx_stim = first(find(moveOnsets > stimOnTime & moveOnsets <= stimOffTime));
                if ~isempty(moveIdx_stim)
                    moveOnset_stim(iTrial) = moveOnsets(moveIdx_stim);
                    moveOffset_stim(iTrial) = moveOffsets(moveIdx_stim);
                    moveDir_stim(iTrial) = -sign(moveAmps(moveIdx_stim));
                    peakVel_stim(iTrial) = vel(t == peakVelTimes(moveIdx_stim));
                else
                    moveOnset_stim(iTrial) = NaN;
                    moveOffset_stim(iTrial) = NaN;
                    moveDir_stim(iTrial) = NaN;
                    peakVel_stim(iTrial) = NaN;
                end
            else
                moveOnset_stim(iTrial) = NaN;
                moveOffset_stim(iTrial) = NaN;
                moveDir_stim(iTrial) = NaN;
                peakVel_stim(iTrial) = NaN;
            end

        end

        %identify the first movement of the trial, redundant info but easy to
        %do
        firstMoveOnset = nan(1,numCompleteTrials);
        firstMoveOnset((~isnan(moveOnset_stim))) = moveOnset_stim(~isnan(moveOnset_stim));
        firstMoveOnset((~isnan(moveOnset_prestim))) = moveOnset_prestim(~isnan(moveOnset_prestim));

        firstMoveOffset = nan(1,numCompleteTrials);
        firstMoveOffset((~isnan(moveOffset_stim))) = moveOffset_stim(~isnan(moveOffset_stim));
        firstMoveOffset((~isnan(moveOffset_prestim))) = moveOffset_prestim(~isnan(moveOffset_prestim));

        firstMoveDir = nan(1,numCompleteTrials);
        firstMoveDir((~isnan(moveDir_stim))) = moveDir_stim(~isnan(moveDir_stim));
        firstMoveDir((~isnan(moveDir_prestim))) = moveDir_prestim(~isnan(moveDir_prestim));

        firstMoveVel = nan(1,numCompleteTrials);
        firstMoveVel((~isnan(peakVel_stim))) = peakVel_stim(~isnan(peakVel_stim));
        firstMoveVel((~isnan(peakVel_prestim))) = peakVel_prestim(~isnan(peakVel_prestim));


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
            'epoch', 'stimulus',...
            'isMoving', isMoving_stim,... 
            'onsetTimes', moveOnset_stim,...
            'offsetTimes', moveOffset_stim,...
            'moveDir', moveDir_stim,...
            'peakVel', peakVel_stim ...
            );

        wheelMoves.epochs(3) = struct(...
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

         disp(char(strcat({'session '},num2str(ex),{' completed'})))
    else
        % % TODO: put workbench (below) here to replace line 206 (nan) % %
        allWheelMoves{ex} = NaN;
        clearvars -except ex expInfo allWheelMoves allEventTimes
        disp(char(strcat({'session '},num2str(ex),{' completed'})))
    end
end



%% WORKBENCH

% %% chop up traces around move on/offsets (for sessions without trial structure)
% 
% timeBefore = 1; %s
% timeAfter = 1;
% 
% for m = 1:length(moveOnsets)
%     onIdx = find(t == moveOnsets(m));
%     offIdx = find(t == moveOffsets(m));
%     try
%         ticIdx = onIdx - timeBefore*Fs;
%     catch
%         ticIdx = 1;
%     end
%     try
%         tocIdx = offIdx + timeAfter*Fs;
%     catch
%         tocIdx = numel(t);
%     end
%     wheelTraceTimes{m} = t(ticIdx:tocIdx);
%     wheelTraceTimes{m} = pos(ticIdx:tocIdx);
% end