function [neuralData] = alignResps(expInfo, neuralData, behavioralData, timeBeforeAndAfter, events)
%Extract the trial-by-trial activity for ROIs in each imaging plane
%and align to a particular trial event (stim on, movement on, etc.)
%
%'event' is a string that is taken from the 'events' field of the eventTimes 
%structure generated in getEventTimes.m

%28 March 2020 changed inputs and outputs to take multiple experiments and
%deliver tidier structs

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

for ex = 1:length(expInfo)
    if ~isstruct(behavioralData(ex).eventTimes) || ~isstruct(behavioralData(ex).wheelMoves)
        alignedResps = NaN;
    else
        
        nt = numel(expInfo(ex).block.events.endTrialTimes);

        for ev = 1:length(events)

            % get event time window for each trial
            eventWindow = -timeBefore:Fs:timeAfter;
            if strcmp(events(ev), 'firstMoveTimes')
                firstMoveTimes = min([behavioralData(ex).wheelMoves.epochs(2).onsetTimes; behavioralData(ex).wheelMoves.epochs(3).onsetTimes]);
                periEventTimes = bsxfun(@plus,firstMoveTimes',eventWindow);
            else
                periEventTimes = bsxfun(@plus,behavioralData(ex).eventTimes(strcmp({behavioralData(ex).eventTimes.event},events{ev})).daqTime',eventWindow);
            end
            periEventTimes = periEventTimes(1:nt,:);

            %initialize alignedResp cell
            alignedResps{ev} = zeros(nt, length(eventWindow), size(neuralData(ex).cellResps,2));

            %grab cell responses associated with the event time windows 
            %(size is nTrials x windowLength x nCells)
            for iCell = 1:size(neuralData(ex).cellResps,2)        
                alignedResps{ev}(:,:,iCell) = interp1(neuralData(ex).respTimes,neuralData(ex).cellResps(:,iCell),periEventTimes,'previous');
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
    neuralData(ex).eta.alignedResps = alignedResps;
%     neuralData(ex).eta.alignedPCs = alignedPCs;
    neuralData(ex).eta.events = events;
    neuralData(ex).eta.eventWindow = eventWindow;
    
end    
    
    
    
    
    
    
    
    
    
    
    