function [eyeData] = alignMotionEnergy(expInfo, eyeData, behavioralData, timeBeforeAndAfter, events)
%Extract the trial-by-trial activity for ROIs in each imaging plane
%and align to a particular trial event (stim on, movement on, etc.)
%
%'event' is a string that is taken from the 'events' field of the eventTimes 
%structure generated in getEventTimes.m

Fs = 0.02;

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
end

for ex = 1:length(expInfo)
    
    expRef = data.constructExpRef(...
        expInfo(ex).mouseName,...
        expInfo(ex).expDate,...
        num2str(expInfo(ex).expNum));
    eyefile = strcat(expRef,'_eye_proc.mat');
    procdir = fullfile(...
        'G:\Data\Video\',...
        expInfo(ex).mouseName,...
        expInfo(ex).expDate,...
        num2str(expInfo(ex).expSeries),...
        eyefile);
    proc = load(procdir);
    
    vidIdx = 1:length(eyeData.timeInterp);
    
    et = behavioralData(ex).eventTimes;
    wm = behavioralData(ex).wheelMoves;
    nt = numel(expInfo(ex).block.events.endTrialTimes);
    nr = 0;
    
    %extract traces from the processed movie file
    
    vidIdx_int = interp1(eyeData(ex).timeAligned,vidIdx, expInfo(ex).Timeline.rawDAQTimestamps,'nearest','extrap');
    
    for f = 1:size(proc.motSVD_1, 2)
        motSVD_int(:,f) = interp1(eyeData(ex).timeAligned,proc.motSVD_1(:,f), expInfo(ex).Timeline.rawDAQTimestamps,'nearest','extrap');
        movSVD_int(:,f) = interp1(eyeData(ex).timeAligned,proc.movSVD_1(:,f), expInfo(ex).Timeline.rawDAQTimestamps,'nearest','extrap');
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
        alignedMot{ev} = zeros(nt, length(eventWindow), size(motSVD_int,2));
        alignedMov{ev} = zeros(nt, length(eventWindow), size(motSVD_int,2));
        
        %grab cell responses associated with the event time windows 
        %(size is nTrials x windowLength x nCells)
        for iCell = 1:size(motSVD_int,2)     
            alignedMot{ev}(:,:,iCell) = interp1(expInfo(ex).Timeline.rawDAQTimestamps,motSVD_int(:,iCell),periEventTimes,'nearest','extrap');
            alignedMov{ev}(:,:,iCell) = interp1(expInfo(ex).Timeline.rawDAQTimestamps,motSVD_int(:,iCell),periEventTimes,'nearest','extrap');
        end   

    end
    
    % save data into the neuralData struct
    eyeData(ex).eta.alignedMot = alignedMot;
    eyeData(ex).eta.alignedMov = alignedMov;    
end    

%% SVD
% for ex = 1:length(expInfo)
%     
%     et = behavioralData(ex).eventTimes;
%     wm = behavioralData(ex).wheelMoves;
%     nt = numel(expInfo(ex).block.events.endTrialTimes);
%     nr = 0;
%     npc = 50;
%     %extract traces from the processed movie file
%     faceROIs = [];
%     faceROIs_int = zeros(nr,length(expInfo(ex).Timeline.rawDAQTimestamps),npc); 
%     
%     if isfield(eyeData(ex).proc,'pupil')
%         nr = nr + 1;
%         faceROIs(nr,:,:) = nan(size(eyeData.proc.face{1, 1}.motionSVD,1),npc);
%     end
%     if isfield(eyeData(ex).proc,'face')
%         for t = 1:length(eyeData(ex).proc.face)
%             if ~isempty(eyeData(ex).proc.face{t})
%                 nr = nr + 1;
%                 faceROIs(nr,:,:) = eyeData(ex).proc.face{t}.motionSVD(:,1:npc);
%             end
%         end
%     end
%     
%     for f = 1:size(faceROIs, 1)
%         samplePoints = {eyeData(ex).timeAligned, 1:npc};
%         F = griddedInterpolant(samplePoints,squeeze(faceROIs(f,:,:)));
%         queryPoints = {expInfo(ex).Timeline.rawDAQTimestamps, 1:npc};
%         Cq = F(queryPoints);
%         faceROIs_int(f,:,:) = Cq;
%     end
%     
%     for ev = 1:length(events)
%         if strcmp(events(ev), 'firstMoveTimes')
%             startTimes = wm.epochs(5).onsetTimes - timeBefore;
%             endTimes = wm.epochs(5).onsetTimes + timeAfter;
%         else
%             startTimes = et(strcmp({et.event},events{ev})).daqTime - timeBefore;
%             endTimes = et(strcmp({et.event},events{ev})).daqTime - timeAfter;
%         end
%         
%         % get event time window for each trial
%         eventWindow = -timeBefore:Fs:timeAfter;
%         if strcmp(events(ev), 'firstMoveTimes')
%             firstMoveTimes = wm.epochs(5).onsetTimes;
%             periEventTimes = bsxfun(@plus,firstMoveTimes',eventWindow);
%         else
%             periEventTimes = bsxfun(@plus,et(strcmp({et.event},events{ev})).daqTime',eventWindow);
%         end
%         periEventTimes = periEventTimes(1:nt,:);
% 
%         %initialize alignedResp cell
%         alignedFaceSVD{ev} = zeros(nt, length(eventWindow), npc, nr);
%         
%         %grab cell responses associated with the event time windows 
%         %(size is nTrials x windowLength x nCells)
%         for iROI = 1:nr 
%             for iPC = 1:npc
%                 alignedFaceSVD{ev}(:,:,iPC,iROI) = interp1(expInfo(ex).Timeline.rawDAQTimestamps,faceROIs_int(iROI,:,iPC),periEventTimes,'nearest','extrap');
%             end
%         end   
% 
%     end
%     
%     % save data into the neuralData struct
%     eyeData(ex).eta.alignedFaceSVD = alignedFaceSVD;
%     
% end    
    
   
    
  
end
    
    
    