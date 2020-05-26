function [baselineResps, stimResps, pmovResps, movResps, rewResps] = getEpochResps(eta)

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
baselineIdx = stim_eventIdx + stim_preTime(1) : stim_eventIdx;

%compute the mean baseline activity per cell, per trial (trials x neurons)
baselineResps = squeeze(mean(stim_alignedTraces(:,baselineIdx,:),2));


%% compute peristimulus activity

%designate a peristimulus window
stimTime = [0 0.3] / Fs;
stimIdx = stim_eventIdx + stimTime(1) :stim_eventIdx + stimTime(2);

%compute the mean peristimulus activity per cell, per trial (trials x neurons)
stimResps = squeeze(mean(stim_alignedTraces(:,stimIdx,:),2));


%% compute perimovement activity

% align traces to movement onset
event = 'firstMoveTimes';
mov_alignedTraces = alignedResps{strcmp(events,event)};
mov_eventWindow = eventWindow;

%designate a movement window
mov_eventIdx = find(mov_eventWindow == 0);
movTime = [-0.2 0.1] / Fs;
movIdx = mov_eventIdx + movTime(1) : mov_eventIdx + movTime(2);

%compute the mean perimovement activity per cell, per trial (trials x neurons)
movResps = squeeze(mean(mov_alignedTraces(:,movIdx,:),2));


%% compute premovement activity

%designate a movement window
pmov_eventIdx = find(mov_eventWindow == 0);
pmovTime = [-0.7 -0.1] / Fs;
pmovIdx = pmov_eventIdx + pmovTime(1) : pmov_eventIdx + pmovTime(2);

%compute the mean perimovement activity per cell, per trial (trials x neurons)
pmovResps = squeeze(mean(mov_alignedTraces(:,pmovIdx,:),2));


%% compute perireward activity

% align traces to movement onset
event = 'feedbackTimes';
rew_alignedTraces = alignedResps{strcmp(events,event)};
rew_eventWindow = eventWindow;

%designate a movement window
rew_eventIdx = find(rew_eventWindow == 0);
rewTime = [0 0.2] / Fs;
rewIdx = rew_eventIdx + rewTime(1) : rew_eventIdx + rewTime(2);

%compute the mean perireward activity per cell, per trial (trials x neurons)
rewResps = squeeze(mean(rew_alignedTraces(:,rewIdx,:),2));
