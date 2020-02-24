function [alignedResps, eventWindow] = alignResps(expInfo, cellResps, respTimes, eventTimes, event, Fs)
%Extract the trial-by-trial activity for ROIs in each imaging plane
%and align to a particular trial event (stim on, movement on, etc.)
%
%'event' is a string that is taken from the 'events' field of the eventTimes 
%structure generated in getEventTimes.m

numCompleteTrials = numel(expInfo.block.events.endTrialTimes);

%Upsampling rate
if nargin < 6
    Fs = 0.1;
end

% get event time window for each trial
timeBefore = 2; %seconds
timeAfter = 2;
eventWindow = -timeBefore:Fs:timeAfter;
periEventTimes = bsxfun(@plus,eventTimes(strcmp({eventTimes.event},event)).daqTime',eventWindow);
periEventTimes = periEventTimes(1:numCompleteTrials,:);

%grab cell responses associated with the event time windows 
%(size is nTrials x windowLength x nCells)
alignedResps = zeros(numCompleteTrials, length(eventWindow), size(cellResps,2));
for iCell = 1:size(cellResps,2)        
    alignedResps(:,:,iCell) = interp1(respTimes,cellResps(:,iCell),periEventTimes,'previous');
end