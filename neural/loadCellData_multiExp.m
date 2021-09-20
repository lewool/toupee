function [allFcell,expInfo] = loadCellData_multiExp(expInfo)
%Load the calcium for curated ROIs across
%a series of experiments processed together in Suite2P
%
%DOES NOT LOAD FLYBACK PLANE
%
% 24 Feb 2020 LEW added option to import processed data from python suite2p
% 26 Mar 2020 LEW omitted flyblack from load
% 27 Mar 2020 LEW adapted for multi expInfo load

for ex = 1:length(expInfo)
%% LOAD EXPINFO VARIABLES

mouseName = expInfo(ex).mouseName;
expDate = expInfo(ex).expDate;
expSeries = expInfo(ex).expSeries;

paths = data.dataPaths();
procDir = paths.local{1}{1};
dataLocations = paths.server;

fileName = strcat('Fall.mat');

dat = load(fullfile(procDir,mouseName,expDate,'suite2p','plane0',fileName));
expInfo(ex).numPlanes = dat.ops.nplanes;
expInfo(ex).numChannels = dat.ops.nchannels;
procType = 'python';
%% segment concatenated tiffs into individual sessions

seriesFolder = strrep(num2str(expSeries),'  ',' ');
seriesFolder = strrep(num2str(seriesFolder),'  ',' ');
seriesFolder = strrep(num2str(seriesFolder),' ','_');

    [expRef, expLog] = data.constructExpRef(mouseName,expDate,expSeries);
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

    if ~exist('Timeline')
        warning('No timeline file was found')
        Timeline = [];
        server = [];
    end
    
    % get the plane frame times per session
    planeInfo = getPlaneFrameTimes(expInfo(ex));
    for p = 1:length(planeInfo)
        sessionFrames(p,ex) = length(planeInfo(p).frameTimes);
    end
end

%%
for ex = 1:length(expInfo)

sessionCounter = [ones(expInfo(ex).numPlanes,1) sessionFrames];
sessionStart = sum(sessionCounter(:,1:ex),2);
sessionEnd = sessionStart + sessionCounter(:,ex+1)-1;

for iPlane = 2:expInfo(ex).numPlanes
    folderName = strcat('plane',num2str(iPlane-1));
    dat = load(fullfile(procDir,mouseName,expDate,'suite2p',folderName,fileName));

    FcellAll = dat.F(:,sessionStart(iPlane):sessionEnd(iPlane));
    FcellNeuAll = dat.Fneu(:,sessionStart(iPlane):sessionEnd(iPlane));
    neuropilCoeff = dat.ops.neucoeff;
    spikesAll = dat.spks(:,sessionStart(iPlane):sessionEnd(iPlane));
    meanImg = dat.ops.meanImgE;

    %extract the dat.stat.iscell values as logicals to create an indexing vector
    iscellIdx = logical(dat.iscell(:,1));

    %reshape Fcell and FcellNeu to include only ROIs you selected manually during Suite2P
    Fcell = FcellAll(iscellIdx,:);
    FcellNeu = FcellNeuAll(iscellIdx,:);
    FcellCorrected = Fcell - (neuropilCoeff * FcellNeu);

    % compute dF/F for each ROI trace/neuropil, plus all neuropil
    for j = 1:size(Fcell,1)
        FcellCorrectedDFF(j,:) = dff(FcellCorrected(j,:),0.1,20,30,expInfo(ex).numPlanes);
        FcellNeuDFF(j,:) = dff(FcellNeu(j,:),0.1,20,30,expInfo(ex).numPlanes);
    end

    for k = 1:size(FcellNeuAll,1)
        FcellNeuAllDFF(k,:) = dff(FcellNeuAll(k,:),0.1,20,30,expInfo(ex).numPlanes);
    end

    % extract the deconvolved spikes as well
    spikes = spikesAll(iscellIdx,:);

    planeTraces.Fcell = Fcell;
    planeTraces.FcellNeu = FcellNeu;
    planeTraces.neuropilCoefficient = neuropilCoeff;
    planeTraces.FcellCorrected = FcellCorrected;
    planeTraces.FcellCorrectedDFF = FcellCorrectedDFF;
    planeTraces.FcellNeuDFF = FcellNeuDFF;
    planeTraces.FcellNeuAvg = mean(FcellNeuAll);
    planeTraces.FcellNeuAvgDFF = mean(FcellNeuAllDFF);
    planeTraces.spikes = spikes;
    planeTraces.meanImage = meanImg;

    clearvars FcellAll FcellNeuAll Fcell FcellNeu neuropilCoeff FcellCorrected FcellCorrectedDFF FcellNeuDFF FcellNeuAllDFF spikes

    %a single (corrected) cell trace is found in allPlaneTraces(iPlane).FcellCorrected(iCell,:)
    allFcellInd(iPlane-1) = planeTraces;
    clear planeTraces;
end
allFcell{ex} = allFcellInd;

end
