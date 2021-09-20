function planeInfo = getPlaneFrameTimes(expInfo)

% A simple function for retrieving the plane-by-plane frame times from
% Timeline's 'neuralFrames'
% LEW 16 June 2018
% 20 July 2021: Massive update to accommodate multiExp data in suite2p

mouseName = expInfo.mouseName;
expDate = expInfo.expDate;
expNum = expInfo.expNum;
numPlanes = double(expInfo.numPlanes);

[expRef, expLog] = data.constructExpRef(mouseName,expDate,expNum);

% load all possible data locations
paths = data.dataPaths();
dataLocations = paths.server;

% load the timeline file
for i = 1:length(dataLocations)
    timelineFilePath = data.makeExpFilePath(expRef, expLog, dataLocations{i}, 'timeline');
    try
        load(timelineFilePath);
        if exist('Timeline')
            server = char(dataLocations{i});
            break
        end
    catch
    end
end

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

%preliminary plane frame times
for iPlane = 1:numPlanes
    planeIdx = iPlane:numPlanes:numel(allNeuralFrameTimes);
    planeInfo(iPlane).planeFrames = planeIdx';
    planeInfo(iPlane).frameTimes = allNeuralFrameTimes(planeIdx);
end

% now we go to the tiffs themselves to adjust exact number...
if ~isa(numPlanes,'double')
    double(numPlanes)
end

% go to session folder
expDir = data.constructExpDir(expRef, expLog, server);
cd(expDir)

% get total number of actual tiff frames for the session
ntiffs = length(dir('*.tif'));
maxFramesPerTiff = 2000;

tfirst = dir(strcat(expRef,'_2P_00001_00001.tif'));
maxFileSize = tfirst.bytes;
try
    tlast = dir(strcat(expRef,'_2P_00001_000',num2str(ntiffs),'.tif'));
catch
    tlast = dir(strcat(expRef,'_2P_00001_0000',num2str(ntiffs),'.tif'));
end
lastFileSize = tlast.bytes;
nframes = maxFramesPerTiff * (ntiffs-1) + floor(maxFramesPerTiff*(lastFileSize/maxFileSize));


% make a vector of frame labels (assuming start at plane 1 and no drops)
framelabels = repmat([1:numPlanes],[1 ceil(ntiffs*maxFramesPerTiff/numPlanes)]);

% figure out the plane of the last frame of the session
lastFramePlane = framelabels(nframes);

baseLength = floor(nframes/numPlanes);
for iPlane = 1:numPlanes
    if iPlane <= lastFramePlane
        planeLength(iPlane) = baseLength + 1;
    else
        planeLength(iPlane) = baseLength;
    end
end

%fine adjust (if possible...old single sessions fail)
try
    for iPlane = 1:numPlanes
        planeInfo(iPlane).planeFrames = planeInfo(iPlane).planeFrames(1:planeLength(iPlane));
        planeInfo(iPlane).frameTimes = planeInfo(iPlane).frameTimes(1:planeLength(iPlane));
    end
catch
end