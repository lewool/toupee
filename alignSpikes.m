function [alignedTraces, alignedSpikes, eventWindow] = alignSpikes(mouseName, expDate, expNum, expSeries, allFcell, eventTimes, ops, event, cellIdx, day)
%Extract the trial-by-trial activity for ROIs in each imaging plane
%and align to a particular trial event (stim on, movement on, etc.)
%
%'event' is a string that is taken from the 'events' field of the eventTimes 
%structure generated in getEventTimes.m
%
%'alignedTraces' is a 1 x n struct of calcium, neuropil, and spike activity 
%for each ROI, where n is the number of planes. Each cell in the struct has
%dimension trial x time x ROI no. (ROI calcium/spikes) or trial x time (neuropil)
%
%A small amount of upsampling is done so all ROIs across planes have the
%same timepoints (prestimTimes, periEventTimes)
%
% OPTIONAL: 'cellIdx' and 'day'
%'cellIdx' is an indexing struct (generated via registers2p.m and 
%daisyChainDays.m) that can be used to pull out only cells imaged across
%all days of a multiday, registered experiment. 'day' is specified with a 
%scalar that refers to a column of cellIdx{plane}. CAUTION: Use indexing 
%that assumes only 'iscell' ROIs from Suite2P as we already excluded 
%non-iscell ROI traces in 'loadExpTraces.m'
%
% 9 July 2018 Added cell indexing
% 7 Dec 2018 Edited interpolation computation
% 10 June 2019 Stripped out everything but session spike alignment



%% LOAD DATA
try
    cd('C:\Users\Wool\Documents\MATLAB\expPipeline\expPipeline');
    [block, Timeline] = loadData(mouseName, expDate, expNum);
catch
    cd('C:\Users\Wool\Documents\GitHub\toupee\behavior')
    [block, Timeline] = loadData(mouseName, expDate, expNum);
end


%% GET FRAME TIMES
 
planeInfo = getPlaneFrameTimes(Timeline, ops.numPlanes);

%% SPLIT THE CA TRACES BY PLANE / TRIAL / CELL

numCompleteTrials = numel(block.events.endTrialTimes);

%Upsampling rate
Fs = 15;

% prestimulus quiescence window
framesBefore = 15; %seconds
framesAfter = -1;
prestimWindow = -framesBefore/Fs:1/Fs:framesAfter/Fs;
stimOn = 'stimulusOnTimes';
prestimTimes = bsxfun(@plus,eventTimes(strcmp({eventTimes.event},stimOn)).daqTime',prestimWindow);
prestimTimes = prestimTimes(1:numCompleteTrials,:);

% event window
framesBefore = 30; %seconds
framesAfter = 30;
eventWindow = -framesBefore/Fs:1/Fs:framesAfter/Fs;

% event = 'stimulusOnTimes'; %stimulusOnTimes interactiveOnTimes stimulusOffTimes rewardOnTimes prestimulusQuiescenceEndTimes
periEventTimes = bsxfun(@plus,eventTimes(strcmp({eventTimes.event},event)).daqTime',eventWindow);
periEventTimes = periEventTimes(1:numCompleteTrials,:);

for iPlane = 1:ops.numPlanes
    % retrieve the relevant plane's cell traces
    planeSpikes = zscore(double(allFcell(iPlane).spikes{1,find(expSeries == expNum)})')';

    % retrieve the frame times for this plane's cells
    planeFrameTimes = planeInfo(iPlane).frameTimes;
    if size(planeFrameTimes,2) ~= size(planeSpikes,2)
        planeFrameTimes = planeFrameTimes(1:size(planeSpikes,2));
    end
    
    for iCell = 1:size(planeSpikes,1)        
        alignedTraces{iPlane}.eventSpikes(:,:,iCell) = interp1(planeFrameTimes,planeSpikes(iCell,:),periEventTimes,'previous');
    end
    
end

alignedSpikes = [];
for iPlane = 2:ops.numPlanes
    alignedSpikes = cat(3,alignedSpikes,alignedTraces{iPlane}.eventSpikes);
end

%% NORMALIZE THE TRACES TRIAL BY TRIAL


% for iPlane = 1:ops.numPlanes
%     normNeuropil = zeros(size(alignedTraces{iPlane}.eventNeuropil));
%     normCellCalcium = zeros(size(alignedTraces{iPlane}.eventCalcium));
%     
%     for iTrial = 1:size(alignedTraces{iPlane}.preStimulusCalcium,1)
%         %normalize neuropil
%         singleTrialNeuropil = alignedTraces{iPlane}.eventNeuropil(iTrial,:);
%         baselineNeuropil = mean(mean(alignedTraces{iPlane}.preStimulusNeuropil(iTrial,:),1));
%         normNeuropil(iTrial,:) = (singleTrialNeuropil - baselineNeuropil)/baselineNeuropil;
%         
%         %normalize cell calcium
%         for iCell = 1:size(alignedTraces{iPlane}.preStimulusCalcium,3)
%             singleTrialCalcium = alignedTraces{iPlane}.eventCalcium(iTrial,:,iCell);
%             baselineCalcium = mean(mean(alignedTraces{iPlane}.preStimulusCalcium(iTrial,:,iCell),1));
%             normCellCalcium(iTrial,:,iCell) = (singleTrialCalcium - baselineCalcium)/baselineCalcium;
%         end
%     end
% 
%     alignedTraces{iPlane}.normEventNeuropil = normNeuropil;
%     alignedTraces{iPlane}.normEventCalcium = normCellCalcium;
% end

%% WORKBENCH

% this is the old way to retrieve plane frame time information, but it's
% very slow and generates huge files. see line 8 for the replacement
% function 'getPlaneFrameTimes'
% 
% infoLoc = fullfile('D:\Data\2P\',mouseName,expDate,num2str(expNum));
% infoFile = fullfile(infoLoc,'allPlaneInfo.mat');
% 
% try 
%     cd(infoLoc);
%     load(infoFile);
% catch
%     try 
%         %catch legacy filename 'allinfo'
%         cd(infoLoc);
%         allPlaneInfo = load(fullfile(infoLoc,'allinfo.mat'),'allinfo');
%         allPlaneInfo = allPlaneInfo.allinfo;
%     catch
%         %generate new file if none exists
%         fprintf('File does not yet exist for this experiment. Generating it now: \n\n');
%         allPlaneInfo = getPlaneFrames(mouseName, expDate, expNum, numPlanes, numChannels);
%         cd(infoLoc);
%         load(infoFile);
%     end
% end
