%%
%Calculate an approximate time to the first wheel movement. This is different from the "timeToFeedback" in that it is based on wheel movement, rather
%than the time when the threshold was reached. WheelMove is an interpolation of the wheel movement (to get it's nearest position at every ms).

%Define a sample rate (sR--used the same as timeline) and resample wheelValues at that rate using 'pchip' and 'extrap' to get wheelPos. Get the
%indices for stimOnset and feedback based on event times and sR.  
sR = 1000;
stimOnsetIdx = round(stimPeriodStart(~timeOuts)*sR);
feedbackIdx = round(feedbackTimes(~timeOuts)*sR);
wheelTime = 1/sR:1/sR:max(x.standardizedBlock.inputs.wheelTimes);
wheelPos = interp1(x.standardizedBlock.inputs.wheelTimes, x.standardizedBlock.inputs.wheelValues', wheelTime', 'pchip', 'extrap')';

%Caluculate "wheelThresh" based on 30% of the median difference in wheel position between closed loop onset and feedback time (30% is a bit arbitrary
%but it seemed to work well). We eliminate trials when the mouse timed out (timeOuts) and didn't move for this. We find the first time that the wheel
%crossed this threshold after stimulus onset ("threshCrsIdx") and the "sign" of that crossing (based on the difference in position from onset)
wheelThresh = median(abs(wheelPos(round(feedbackTimes(~timeOuts)*sR))-wheelPos(round(closedLoopStart(~timeOuts)*sR))))*0.3;
threshCrsIdx = arrayfun(@(x,y) max([nan find(abs(wheelPos(x:y)-wheelPos(x))>wheelThresh,1)+x]), stimOnsetIdx, feedbackIdx);

%Define a summation window (sumWin--51ms) and velocity threshhold. sR*3/sumWin means the wheel will need to move at least 3 "units" in 50ms for this
%to count as a movement initiation. Obviously, the physical distance of a "unit" depends on the rotary encoder. 3 seems to work well for 360 and 1024
%encoders (the only ones I have used). I don't know about 100 encoders.
sumWin = 51;
velThresh  = sR*3/sumWin;

%Get wheel velocity ("wheelVel") from the interpolated wheel, and then use "posVelScan" and "negVelScan" to detect continuous velocity for sumWin
%movements in either direction. Note, here we use forward smoothing, so we are looking for times when the velocity initiates, and then continues for
%the duration of "sumWin". Also note, we introduce huge values whenever the wheel is moving in the opposite direction so that any "sumWin" including a
%move in the opposite direction is eliminated from the scan. Times when the velocity is zero cannot be movement inititations either, so we multiply by
%(wheelVel~=0)
wheelVel = diff([wheelPos(1); wheelPos'])*sR;
posVelScan = conv(wheelVel.*double(wheelVel>0) - double(wheelVel<0)*1e6, [ones(1,sumWin) zeros(1,sumWin-1)]./sumWin, 'same').*(wheelVel~=0);
negVelScan = conv(wheelVel.*double(wheelVel<0) + double(wheelVel>0)*1e6, [ones(1,sumWin) zeros(1,sumWin-1)]./sumWin, 'same').*(wheelVel~=0);

%Identify onsets in borth directions that exceed "velThresh", sort them, and record their sign. Also, extract all times the mouse is "moving"
moveOnsetIdx = [strfind((posVelScan'>=velThresh), [0,1]) -1*strfind((-1*negVelScan'>=velThresh), [0,1])];
movingIdx = sort([find(posVelScan'>=velThresh) find((-1*negVelScan'>=velThresh))]);
[~, srtIdx] = sort(abs(moveOnsetIdx));
moveOnsetSign = sign(moveOnsetIdx(srtIdx));
moveOnsetIdx = abs(moveOnsetIdx(srtIdx));

%"firstMoveTimes" are the first onsets occuring after stimOnsetIdx. Eliminate any that are longer than 1.5s, as these would be timeouts. Also, remove
%onsets when the mouse was aready moving at the time of the stimulus onset (impossible to get an accurate movement onset time in this case)
firstMoveIdx =  arrayfun(@(x) max([nan moveOnsetIdx(find(moveOnsetIdx>x, 1))]), stimOnsetIdx);
errorDetect = (firstMoveIdx-stimOnsetIdx)>1.5*sR | ismember(stimOnsetIdx, movingIdx) | firstMoveIdx>threshCrsIdx;
firstMoveIdx(errorDetect) = nan;

firstMoveDir = arrayfun(@(x) max([nan moveOnsetSign(find(moveOnsetIdx>x, 1))]), stimOnsetIdx);
firstMoveDir(isnan(firstMoveIdx)) = nan;

if p.wheelGain<0
    firstMoveDir = firstMoveDir*-1;
    r.wheelTimeValue = cellfun(@(x) [x(:,1) x(:,2)*-1], r.wheelTimeValue, 'uni', 0);
    warning('Wheel connected in reverse... adjusting.');
end
firstMoveDir = ((firstMoveDir==-1)+1).*(abs(firstMoveDir));

preThreshOnsets = arrayfun(@(x,y) moveOnsetIdx(moveOnsetIdx>=x & moveOnsetIdx<=y), stimOnsetIdx, threshCrsIdx, 'uni', 0);
preThreshDirs = arrayfun(@(x,y) moveOnsetSign(moveOnsetIdx>=x & moveOnsetIdx<=y), stimOnsetIdx, threshCrsIdx, 'uni', 0);
noNAN = ~isnan(firstMoveDir);

reliableMoves = firstMoveDir;
reliableMoves(noNAN) = cellfun(@(x) max(diff([x(1);x(:)]))<50, preThreshOnsets(noNAN)) & cellfun(@(x) length(unique(x))==1, preThreshDirs(noNAN));

theshDirection  = arrayfun(@(x) max([nan moveOnsetSign(find(moveOnsetIdx<x, 1, 'last'))]), threshCrsIdx);
if p.wheelGain<0; theshDirection = theshDirection*-1; warning('Wheel connected in reverse... adjusting.'); end
theshDirection = ((theshDirection==-1)+1).*(abs(theshDirection));
threshReliable = firstMoveDir;
threshReliable(noNAN) = (threshCrsIdx(noNAN)-firstMoveIdx(noNAN))<250;

[firstMoveTime, firstMoveDirection, firstMoveReliable, threshMoveDirection, threshMoveReliable, threshMoveTime] = deal(feedbackValues*nan);
firstMoveTime(~timeOuts) = (firstMoveIdx - stimOnsetIdx)/sR;
firstMoveDirection(~timeOuts) = firstMoveDir;
firstMoveReliable(~timeOuts) = reliableMoves;
threshMoveTime(~timeOuts) = (threshCrsIdx - stimOnsetIdx)/sR;
threshMoveDirection(~timeOuts) = theshDirection;
threshMoveReliable(~timeOuts) = threshReliable;