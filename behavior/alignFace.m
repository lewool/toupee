function [eyeData] = alignFace(expInfo, eyeData, behavioralData, timeBeforeAndAfter, events)
%Extract the trial-by-trial activity for ROIs in each imaging plane
%and align to a particular trial event (stim on, movement on, etc.)
%
%'event' is a string that is taken from the 'events' field of the eventTimes 
%structure generated in getEventTimes.m

Fs = 0.1;

if nargin < 4
    timeBefore = 2;
    timeAfter = 2;
    events = {'stimulusOnTimes' 'firstMoveTimes' 'feedbackTimes'};
elseif nargin < 5
    events = {'stimulusOnTimes' 'firstMoveTimes' 'feedbackTimes'};
    
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
    
    et = behavioralData(ex).eventTimes;
    wm = behavioralData(ex).wheelMoves;
    nt = numel(expInfo(ex).block.events.endTrialTimes);
    nr = 0;
    if isfield(eyeData.proc,'pupil')
        nr = nr + 1;
        faceResps(nr,:) = eyeData.proc.pupil.area;
    end
    if isfield(eyeData.proc,'motSVD')
        for t = 1:length(eyeData.proc.motSVD)
            if ~isempty(eyeData.proc.motSVD{t})
                nr = nr + 1;
                faceResps(nr,:) = eyeData.proc.motSVD{t}(:,1);
            end
        end
    end
    
    for f = 1:size(faceResps, 1)
        faceResps_int(f,:) = interp1(eyeData.time,faceResps(f,:), expInfo(ex).Timeline.rawDAQTimestamps,'nearest','extrap');
    end
    
    for ev = 1:length(events)
        if strcmp(events(ev), 'firstMoveTimes')
            startTimes = wm.epochs(5).onsetTimes - timeBefore;
            endTimes = wm.epochs(5).onsetTimes + timeAfter;
        else
            startTimes = et(strcmp({et.event},events{ev})).daqTime - timeBefore;
            endTimes = et(strcmp({et.event},events{ev})).daqTime - timeAfter;
        end
        
        % get event time window for each trial
        eventWindow = -timeBefore:Fs:timeAfter;
        if strcmp(events(ev), 'firstMoveTimes')
            firstMoveTimes = wm.epochs(5).onsetTimes;
            periEventTimes = bsxfun(@plus,firstMoveTimes',eventWindow);
        else
            periEventTimes = bsxfun(@plus,et(strcmp({et.event},events{ev})).daqTime',eventWindow);
        end
        periEventTimes = periEventTimes(1:nt,:);

        %initialize alignedResp cell
        alignedResps{ev} = zeros(nt, length(eventWindow), nr);
        
        %grab cell responses associated with the event time windows 
        %(size is nTrials x windowLength x nCells)
        for iCell = 1:nr     
            alignedResps{ev}(:,:,iCell) = interp1(expInfo(ex).Timeline.rawDAQTimestamps,faceResps_int(iCell,:),periEventTimes,'nearest','extrap');
        end   

    end
    
    % save data into the neuralData struct
    eyeData(ex).eta.alignedResps = alignedResps;
%     neuralData(ex).eta.alignedPCs = alignedPCs;
    eyeData(ex).eta.events = events;
    eyeData(ex).eta.eventWindow = eventWindow;
    
end    
    
    
    
    
    
    
    
    
    
    
    