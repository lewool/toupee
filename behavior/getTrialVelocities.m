function behavioralData = getTrialVelocities(expInfo, behavioralData)
nt = numel(expInfo.block.events.endTrialValues);
stimOnTimes = behavioralData.eventTimes(1).daqTime;
moveOnTimes = behavioralData.wheelMoves.epochs(5).onsetTimes;
feedbackTimes = behavioralData.eventTimes(5).daqTime;

choiceDir = expInfo.block.events.responseValues(1:nt);
contrasts = expInfo.block.events.contrastValues(1:nt);
uniqueContrasts = getUniqueContrasts(expInfo);

feedbackLatencies = behavioralData.eventTimes(5).daqTime - behavioralData.eventTimes(1).daqTime;
firstMoveLatencies = behavioralData.eventTimes(5).daqTime - behavioralData.wheelMoves.epochs(5).onsetTimes;

plotWindow = -2:0.001:2;
eventIdx = find(plotWindow == 0);
fromStart = plotWindow(1)*1000;
toEnd = plotWindow(end)*1000;
%
for ETA = 1:3
    for t = 1:length(behavioralData.wheelMoves.traces.time)
        if ETA == 1
            traceTime = behavioralData.wheelMoves.traces.time{1,t}  - stimOnTimes(t);
        elseif ETA == 2
            traceTime = behavioralData.wheelMoves.traces.time{1,t}  - moveOnTimes(t);
        elseif ETA == 3
            traceTime = behavioralData.wheelMoves.traces.time{1,t}  - feedbackTimes(t);
        end

        tracePos = behavioralData.wheelMoves.traces.pos{1,t};
        velTime = traceTime;
        vel = behavioralData.wheelMoves.traces.vel{1,t};
        [~,zeroIdx] = min(abs(velTime));
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
    alignedVels{ETA} = bigVel;
    clear bigVel
end

behavioralData.eta.alignedVels = alignedVels;
behavioralData.eta.eventWindow = plotWindow;