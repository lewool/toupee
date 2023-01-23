function [alignedWheel, plotWindow] = getAlignedWheelTraces(expInfo, behavioralData, eta)

stimOnTimes = behavioralData.eventTimes(1).daqTime;
firstMoveTimes = behavioralData.wheelMoves.epochs(5).onsetTimes;
feedbackTimes = behavioralData.eventTimes(5).daqTime;
cueOnTimes = behavioralData.eventTimes(2).daqTime;

if eta == 1
    plotWindow = -.5:0.001:3.5;
    timeCorrect = stimOnTimes;
elseif eta == 2
    plotWindow = -2:0.001:2;
    timeCorrect = firstMoveTimes;
elseif eta == 3
    plotWindow = -2:0.001:2;
    timeCorrect = feedbackTimes;
elseif eta == 4
    plotWindow = -2:0.001:2;
    timeCorrect = cueOnTimes;
end
    
eventIdx = find(plotWindow == 0);
fromStart = plotWindow(1)*1000;
toEnd = plotWindow(end)*1000;
%
for t = 1:length(behavioralData.wheelMoves.traces.time)
    traceTime = behavioralData.wheelMoves.traces.time{1,t}  - timeCorrect(t);
    velTime = traceTime;
    vel = behavioralData.wheelMoves.traces.vel{1,t};
%     zeroIdx = find(velTime == 0);
    [~, zeroIdx] = min(abs(velTime - 0));
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
    
    alignedWheel(t,:) = truncVel;
end