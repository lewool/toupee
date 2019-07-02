function planeInfo = getPlaneFrameTimes(Timeline, numPlanes)

% A simple function for retrieving the plane-by-plane frame times from
% Timeline's 'neuralFrames'
% LEW 16 June 2018

%sampling rate of Timeline
timelineRate = 1000;

%sampling rate of the scope
imagingRate = 30;

%roughly, how many Timeline samples per neuralFrame
neuralFrameRate = timelineRate/imagingRate;

allNeuralFrames = Timeline.rawDAQData(:,strcmp({Timeline.hw.inputs.name},'neuralFrames'));
totalTimelineSamples = length(allNeuralFrames);

[~, neuralFrameIdx, ~] = unique(allNeuralFrames);

% discard the neuralFrames before and after acquisition
discardFrames = find(diff([neuralFrameIdx;totalTimelineSamples]) > ceil(neuralFrameRate));
neuralFrameIdx(discardFrames) = [];

allNeuralFrameTimes = Timeline.rawDAQTimestamps(neuralFrameIdx);

for iPlane = 1:numPlanes
    planeIdx = iPlane:numPlanes:numel(allNeuralFrameTimes);
    planeInfo(iPlane).planeFrames = planeIdx';
    planeInfo(iPlane).frameTimes = allNeuralFrameTimes(planeIdx);
end
    