function [cellResps, respTimes] = interpCellTimes(expInfo, allFcell, Fs)

%% LOAD DATA FROM EXPINFO

mouseName = expInfo.mouseName;
expDate = expInfo.expDate;
expNum = expInfo.expNum;
expSeries = expInfo.expSeries;
block = expInfo.block;
Timeline = expInfo.Timeline;
numPlanes = expInfo.numPlanes;

%% SET UPSAMPLING RATE (if not an input)

if nargin < 3
    Fs = 0.1; 
end

%% GET FRAME TIMES
planeInfo = getPlaneFrameTimes(Timeline, numPlanes);

for p = 1:length(planeInfo)
    minTs(p) = min(planeInfo(p).frameTimes);
    maxTs(p) = max(planeInfo(p).frameTimes);
end

startT = min(minTs);
endT = max(maxTs);

globalTime = startT:Fs:endT;

%%
cellResps = [];
for iPlane = 2:numPlanes
    planeTime = planeInfo(iPlane).frameTimes;
    C = double(allFcell(iPlane).spikes{1}');
    if size(planeTime,2) ~= size(C,1)
        planeTime = planeTime(1:size(C,1));
    end
    samplePoints = {planeTime, 1:size(C,2)};
    F = griddedInterpolant(samplePoints,C);
    queryPoints = {globalTime, 1:size(C,2)};
    Cq = F(queryPoints);
    cellResps(:, size(cellResps,2)+1:size(cellResps,2)+size(Cq,2)) = Cq;
end

respTimes = globalTime;