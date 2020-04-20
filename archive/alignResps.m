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
    events = {'stimulusOnTimes' 'prestimulusQuiescenceEndTimes' 'feedbackTimes'};
elseif nargin < 5
    events = {'stimulusOnTimes' 'prestimulusQuiescenceEndTimes' 'feedbackTimes'};
    
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
    
    nt = numel(expInfo(ex).block.events.endTrialTimes);
    
    for ev = 1:length(events)
        
        % get event time window for each trial
        eventWindow = -timeBefore:Fs:timeAfter;
        periEventTimes = bsxfun(@plus,behavioralData(ex).eventTimes(strcmp({behavioralData(ex).eventTimes.event},events{ev})).daqTime',eventWindow);
        periEventTimes = periEventTimes(1:nt,:);

        %initialize alignedResp cell
        alignedResps{ev} = zeros(nt, length(eventWindow), size(neuralData(ex).cellResps,2));
        
        %grab cell responses associated with the event time windows 
        %(size is nTrials x windowLength x nCells)
        for iCell = 1:size(neuralData(ex).cellResps,2)        
            alignedResps{ev}(:,:,iCell) = interp1(neuralData(ex).respTimes,neuralData(ex).cellResps(:,iCell),periEventTimes,'previous');
        end
    end
    
    % save data into the neuralData struct
    neuralData(ex).eta.alignedResps = alignedResps;
    neuralData(ex).eta.events = events;
    neuralData(ex).eta.eventWindow = eventWindow;
    
end    
    
    
    
    
    
    
    
    
    
    
    