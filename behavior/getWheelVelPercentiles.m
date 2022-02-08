function [velTrialIdx, vels] = getWheelVelPercentiles(expInfo, behavioralData,time)
pt = 4;
if nargin < 3
    time = 'late';
end

[~, timeTrials] = selectCondition(expInfo(1), getUniqueContrasts(expInfo(1)), behavioralData(1), ...
    initTrialConditions('movementTime',time,'responseType','correct'));

allVels = behavioralData(1).wheelMoves.epochs(5).peakVel;

[sortVels, sortIdx] = sort(allVels(timeTrials));

prctLength = floor(length(timeTrials)/pt);
for p = 1:pt
    velTrialIdx{p} = sortIdx(p*prctLength-prctLength+1:p*prctLength-prctLength+prctLength);
    vels{p} = sortVels(p*prctLength-prctLength+1:p*prctLength-prctLength+prctLength);
end
    