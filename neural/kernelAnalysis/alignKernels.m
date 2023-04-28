function [kernelAnalysis] = alignKernels(expInfo, kernelAnalysis, neuralData, behavioralData, timeBeforeAndAfter, events)
%Extract the trial-by-trial activity for ROIs in each imaging plane
%and align to a particular trial event (stim on, movement on, etc.)
%
%'event' is a string that is taken from the 'events' field of the eventTimes 
%structure generated in getEventTimes.m


%set defaults
Fs = 0.1;

if nargin < 4
    timeBefore = 2;
    timeAfter = 2;
    events = {'stimulusOnTimes' 'firstMoveTimes' 'feedbackTimes' 'interactiveOnTimes'};
elseif nargin < 5
    events = {'stimulusOnTimes' 'firstMoveTimes' 'feedbackTimes' 'interactiveOnTimes'};
    
    % call window range; if fails, set default
    try
        timeBefore = timeBeforeAndAfter(1);
        timeAfter = timeBeforeAndAfter(2);
    catch
        try
            timeBefore = timeBeforeAndAfter;
            timeAfter = timeBeforeAndAfter;
        catch
            timeBefore = 2;
            timeAfter = 2;
        end
    end
else
    % call window range; if fails, set default
    try
        timeBefore = timeBeforeAndAfter(1);
        timeAfter = timeBeforeAndAfter(2);
    catch
        try
            timeBefore = timeBeforeAndAfter;
            timeAfter = timeBeforeAndAfter;
        catch
            timeBefore = 2;
            timeAfter = 2;
        end
    end
end

if ~isstruct(behavioralData.eventTimes) || ~isstruct(behavioralData.wheelMoves)
    alignedResps = NaN;
else

    nt = numel(expInfo.block.events.endTrialTimes);

    for ev = 1:length(events)

        % get event time window for each trial
        eventWindow = -timeBefore:Fs:timeAfter;
        if strcmp(events(ev), 'firstMoveTimes')
            firstMoveTimes = min([behavioralData.wheelMoves.epochs(2).onsetTimes; behavioralData.wheelMoves.epochs(3).onsetTimes]);
            periEventTimes = bsxfun(@plus,firstMoveTimes',eventWindow);
        else
            periEventTimes = bsxfun(@plus,behavioralData.eventTimes(strcmp({behavioralData.eventTimes.event},events{ev})).daqTime',eventWindow);
        end
        periEventTimes = periEventTimes(1:nt,:);

        %initialize alignedResp cell
        alignedResps{ev} = zeros(nt, length(eventWindow), size(kernelAnalysis.predictedActivity,2));

        %grab cell responses associated with the event time windows 
        %(size is nTrials x windowLength x nCells)
        for iCell = 1:size(kernelAnalysis.predictedActivity,2)        
            alignedResps{ev}(:,:,iCell) = interp1(neuralData.respTimes,kernelAnalysis.predictedActivity(:,iCell),periEventTimes,'previous');
        end

%         %reshape the matrix to cells x time x trials
%         M = permute(alignedResps{ev},[3 2 1]);
%         [i,~] = find(isnan(M));
%         nanCells = unique(i);
% 
%         for c = fliplr(nanCells')
%             M(c,:,:) = [];
%         end
% 
%         % center the data
%         Mcent = zeros(size(M,1), size(M,2),size(M,3));
%         for iCell = 1:size(M,1)
%             Mslice = squeeze(M(iCell,:,:));
%             MsliceMean = mean(Mslice,1);
%             Mcent(iCell,:,:) = Mslice - MsliceMean;
%         end
%         % Mcent = M;
% 
%         Mresh = reshape(Mcent,size(Mcent,1), size(Mcent,2)*size(Mcent,3));
%         [U,S,V] = svd(Mresh,'econ');
%         newM = S*abs(V');
%         newMresh = reshape(newM,size(M,1),size(M,2), size(M,3));
% 
%         % reshape back into the original dimensions
%         %output is a matrix of size trials x time x PCs
% 
%         alignedPCs{ev} = permute(newMresh,[3 2 1]);

    end
end

% save data into the neuralData struct
kernelAnalysis.eta.alignedResps = alignedResps;
%     neuralData.eta.alignedPCs = alignedPCs;
kernelAnalysis.eta.events = events;
kernelAnalysis.eta.eventWindow = eventWindow;
        
    
    
    
    
    
    
    
    
    
    
    