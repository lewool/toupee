function [eyeData] = alignFaceSVD(expInfo, eyeData, behavioralData, timeBeforeAndAfter, events)
%Extract the trial-by-trial activity for ROIs in each imaging plane
%and align to a particular trial event (stim on, movement on, etc.)
%
%'event' is a string that is taken from the 'events' field of the eventTimes 
%structure generated in getEventTimes.m

Fs = 0.02;

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
    %preallocate epmty matrices of correct size 
    faceROIs = zeros(length(expInfo),length(eyeData(ex).proc.face{1}.motionSVD(1,:)),length(eyeData(ex).proc.face{1}.motionSVD(1,:)));
    faceROIs_int = [];
    
    %extract traces from the processed movie file
    if isfield(eyeData(ex).proc,'pupil')
        nr = nr + 1;
        faceROIs(nr,:,:) = eyeData(ex).proc.pupil.area;
    end
    if isfield(eyeData(ex).proc,'face')
        for t = 1:length(eyeData(ex).proc.face)
            if ~isempty(eyeData(ex).proc.face{t})
                nr = nr + 1;
                faceROIs(nr,:,:) = eyeData(ex).proc.face{t}.motionSVD;
            end
        end
    end
    
    for f = 1:size(faceROIs, 1)
        faceROIs_int(f,:) = interp1(eyeData(ex).timeAligned,faceROIs(f,:), expInfo(ex).Timeline.rawDAQTimestamps,'nearest','extrap');
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
        alignedFace{ev} = zeros(nt, length(eventWindow), nr);
        
        %grab cell responses associated with the event time windows 
        %(size is nTrials x windowLength x nCells)
        for iCell = 1:nr     
            alignedFace{ev}(:,:,iCell) = interp1(expInfo(ex).Timeline.rawDAQTimestamps,faceROIs_int(iCell,:),periEventTimes,'nearest','extrap');
        end   

    end
    
    % save data into the neuralData struct
    eyeData(ex).eta.alignedFace = alignedFace;
%     neuralData(ex).eta.alignedPCs = alignedPCs;
    eyeData(ex).eta.events = events;
    eyeData(ex).eta.eventWindow = eventWindow;
    
end    
    
    
    
    
    
    
    
    
    
    
    