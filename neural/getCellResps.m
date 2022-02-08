function [cellResps, respTimes] = getCellResps(expInfo, allFcell, Fs)

% 26 Mar 2020 flyback now omitted during loadCellData, not here

for ex = 1:length(expInfo)
%% LOAD DATA FROM EXPINFO

mouseName = expInfo(ex).mouseName;
expDate = expInfo(ex).expDate;
expNum = expInfo(ex).expNum;
expSeries = expInfo(ex).expSeries;
block = expInfo(ex).block;
Timeline = expInfo(ex).Timeline;
numPlanes = expInfo(ex).numPlanes;

%% SET UPSAMPLING RATE (if not an input)

if nargin < 3
    Fs = 0.1; 
end

%% GET FRAME TIMES
planeInfo = getPlaneFrameTimes(expInfo(ex));

for p = 1:length(planeInfo)
    minTs(p) = min(planeInfo(p).frameTimes);
    maxTs(p) = max(planeInfo(p).frameTimes);
end

startT = min(minTs);
endT = max(maxTs);

globalTime = startT:Fs:endT;

%%
resps = [];
for iPlane = 1:numPlanes-1
    planeTime = planeInfo(iPlane).frameTimes;
    try
        %matlab suite2p output
        C = double(allFcell{ex}(iPlane).spikes{1}');
    catch
        %python suite2p output
        C = double(allFcell{ex}(iPlane).spikes');
    end
    if size(planeTime,2) ~= size(C,1)
        planeTime = planeTime(1:size(C,1));
    end
    samplePoints = {planeTime, 1:size(C,2)};
    F = griddedInterpolant(samplePoints,C);
    queryPoints = {globalTime, 1:size(C,2)};
    Cq = F(queryPoints);
    resps(:, size(resps,2)+1:size(resps,2)+size(Cq,2)) = Cq;
end

    
% cellResps{ex} = zscore(resps);
% cellResps{ex} = resps./max(resps);
cellResps{ex} = resps./prctile(resps,99);

%get rid of bad cells
someNaNs = find(isnan(cellResps{ex}));
badCells = unique(ceil(someNaNs/size(cellResps{ex},1)));
cellResps{ex}(:,badCells) = [];
respTimes{ex} = globalTime;
end