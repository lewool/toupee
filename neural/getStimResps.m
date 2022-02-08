function [baselineResps, stimResps] = getStimResps(eta)

Fs = 0.1;
alignedResps = eta.alignedResps;
events = eta.events;
eventWindow = eta.eventWindow;

%% compute baseline activity

% align traces to stim onset
event = 'stimulusOnTimes';
stim_alignedTraces = alignedResps{strcmp(events,event)};
stim_eventWindow = eventWindow;

%designate a baseline window
stim_eventIdx = find(stim_eventWindow == 0);
stim_preTime = [-0.5 0] / Fs;
baselineIdx = stim_eventIdx + stim_preTime(1) : stim_eventIdx - 1;

%compute the mean baseline activity per cell, per trial (trials x neurons)
baselineResps = squeeze(mean(stim_alignedTraces(:,baselineIdx,:),2));


%% compute peristimulus activity

%designate a peristimulus window
stimTime = [0 0.5] / Fs;
stimIdx = stim_eventIdx + stimTime(1) :stim_eventIdx + stimTime(2);

%compute the mean peristimulus activity per cell, per trial (trials x neurons)
stimResps = squeeze(mean(stim_alignedTraces(:,stimIdx,:),2));