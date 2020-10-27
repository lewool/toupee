function neuralData = getSignificantActivity(expInfo, behavioralData, neuralData, matched)
% Determine which neurons are active for different aspects of the task:
% stimulus onset, movement onset, movement direction, reward onset, block
% identity (ITI)

Fs = 0.1;
    
if matched == 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%% UNMATCHED STATS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    for ex = 1:length(expInfo)
        
    %%%%%%%%%%%%%%%% compute baseline activity

    % align traces to stim onset
    event = 'stimulusOnTimes';
    stim_alignedTraces = neuralData(ex).eta.alignedResps{strcmp(neuralData(ex).eta.events,event)};
    stim_eventWindow = neuralData(ex).eta.eventWindow;

    %designate a baseline window
    stim_eventIdx = find(stim_eventWindow == 0);
    stim_preTime = [-0.5 0] / Fs;
    baselineIdx = stim_eventIdx + stim_preTime(1) : stim_eventIdx;

    %compute the mean baseline activity per cell, per trial (trials x neurons)
    baselineResps = squeeze(mean(stim_alignedTraces(:,baselineIdx,:),2));

    %%%% compute peristimulus activity

    %designate a peristimulus window
    stimTime = [0 0.3] / Fs;
    stimIdx = stim_eventIdx + stimTime(1) :stim_eventIdx + stimTime(2);

    %compute the mean peristimulus activity per cell, per trial (trials x neurons)
    stimResps = squeeze(mean(stim_alignedTraces(:,stimIdx,:),2));

    %%%%%%%%%%%%%%%% compute perimovement activity

    % align traces to movement onset
    event = 'firstMoveTimes';
    mov_alignedTraces = neuralData(ex).eta.alignedResps{strcmp(neuralData(ex).eta.events,event)};
    mov_eventWindow = neuralData(ex).eta.eventWindow;

    %designate a movement window
    mov_eventIdx = find(mov_eventWindow == 0);
    movTime = [-0.2 0.1] / Fs;
    movIdx = mov_eventIdx + movTime(1) : mov_eventIdx + movTime(2);

    %compute the mean perimovement activity per cell, per trial (trials x neurons)
    movResps = squeeze(mean(mov_alignedTraces(:,movIdx,:),2));

    %%%%%%%%%%%%%%%% compute premovement activity

    %designate a movement window
    pmov_eventIdx = find(mov_eventWindow == 0);
    pmovTime = [-0.7 -0.2] / Fs;
    pmovIdx = pmov_eventIdx + pmovTime(1) : pmov_eventIdx + pmovTime(2);

    %compute the mean perimovement activity per cell, per trial (trials x neurons)
    pmovResps = squeeze(mean(mov_alignedTraces(:,pmovIdx,:),2));

    %%%%%%%%%%%%%%%% compute perireward activity

    % align traces to feedback onset
    event = 'feedbackTimes';
    rew_alignedTraces = neuralData(ex).eta.alignedResps{strcmp(neuralData(ex).eta.events,event)};
    rew_eventWindow = neuralData(ex).eta.eventWindow;

    %designate a movement window
    rew_eventIdx = find(rew_eventWindow == 0);
    rewTime = [0 0.2] / Fs;
    rewIdx = rew_eventIdx + rewTime(1) : rew_eventIdx + rewTime(2);

    %compute the mean perireward activity per cell, per trial (trials x neurons)
    rewResps = squeeze(mean(rew_alignedTraces(:,rewIdx,:),2));

    %%%%%%%%%%%%%%%% organize trial types for subsequent comparison

    contrasts = unique(expInfo(ex).block.events.contrastValues);
    et = behavioralData(ex);
    %stimulus
    [~, stimTrials] = selectCondition(expInfo(ex), contrasts, et, initTrialConditions('movementTime','late'));
    [~, leftstimTrials] = selectCondition(expInfo(ex), contrasts(contrasts<0), et, initTrialConditions('movementTime','late'));
    [~, rightstimTrials] = selectCondition(expInfo(ex), contrasts(contrasts>0), et, initTrialConditions('movementTime','late'));
    [~, zerostimTrials] = selectCondition(expInfo(ex), 0, et, initTrialConditions('movementTime','late'));

    %movement
    [~, movleftTrials] = selectCondition(expInfo(ex), contrasts, et, initTrialConditions('movementDir','cw'));
    [~, movrightTrials] = selectCondition(expInfo(ex), contrasts, et, initTrialConditions('movementDir','ccw'));

    %reward
    [~, correctTrials] = selectCondition(expInfo(ex), contrasts, et, initTrialConditions('responseType','correct'));
    [~, incorrectTrials] = selectCondition(expInfo(ex), contrasts, et, initTrialConditions('responseType','incorrect'));

    %value
    [~, highleftTrials] = selectCondition(expInfo(ex), contrasts, et, initTrialConditions('movementTime','late','highRewardSide','left'));
    [~, highrightTrials] = selectCondition(expInfo(ex), contrasts, et, initTrialConditions('movementTime','late','highRewardSide','right'));


    %%%%%%%%%%%%%%%% Wilcoxon tests

    labels = {'stim', 'leftStim', 'rightStim', 'mov', 'leftMov', 'rightMov', 'reward', 'value','advanceMov'};
    pValues = zeros(size(baselineResps,2),length(labels));
    for iCell = 1:size(baselineResps,2)
        pValues(iCell,1) = signrank(stimResps(stimTrials,iCell),baselineResps(stimTrials,iCell),'tail','right');
        pValues(iCell,2) = ranksum(stimResps(leftstimTrials,iCell),stimResps(zerostimTrials,iCell),'tail','right');
        pValues(iCell,3) = ranksum(stimResps(rightstimTrials,iCell),stimResps(zerostimTrials,iCell),'tail','right');
        pValues(iCell,4) = signrank(baselineResps(:,iCell),movResps(:,iCell));
        pValues(iCell,5) = ranksum(movResps(movleftTrials,iCell),movResps(movrightTrials,iCell),'tail','right');
        pValues(iCell,6) = ranksum(movResps(movrightTrials,iCell),movResps(movleftTrials,iCell),'tail','right');
        pValues(iCell,7) = ranksum(rewResps(correctTrials,iCell),rewResps(incorrectTrials,iCell));
        pValues(iCell,8) = ranksum(stimResps(highleftTrials,iCell),stimResps(highrightTrials,iCell)); 
        pValues(iCell,9) = signrank(pmovResps(:,iCell),baselineResps(:,iCell)); 
    end
    
    %%%%%%%%%%%%%%%% Bonferroni correction

    alpha = 0.05;
    bfc = size(pValues,2);
    bfcAlpha = alpha/bfc;
    bfcH = pValues < bfcAlpha;
    
    %%%%%%%%%%%%%%%% collect into struct
    neuralData(ex).stats = struct(...
        'pValues',pValues,...
        'labels',{labels},...
        'alpha',alpha,...
        'bfcAlpha',bfcAlpha,...
        'bfcH',bfcH...
        );
    
    end

elseif matched == 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%% MATCHED STATS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%% compute baseline activity
    
    % align traces to stim onset
    event = 'stimulusOnTimes';
    stim_alignedTraces = neuralData.eta.alignedResps{strcmp(neuralData.eta.events,event)};
    stim_eventWindow = neuralData.eta.eventWindow;

    %designate a baseline window
    stim_eventIdx = find(stim_eventWindow == 0);
    stim_preTime = [-0.5 0] / Fs;
    baselineIdx = stim_eventIdx + stim_preTime(1) : stim_eventIdx;

    %compute the mean baseline activity per cell, per trial (trials x neurons)
    baselineResps = squeeze(mean(stim_alignedTraces(:,baselineIdx,:),2));


    %%%%%%%%%%%%%%%% compute peristimulus activity

    %designate a peristimulus window
    stimTime = [0 0.3] / Fs;
    stimIdx = stim_eventIdx + stimTime(1) :stim_eventIdx + stimTime(2);

    %compute the mean peristimulus activity per cell, per trial (trials x neurons)
    stimResps = squeeze(mean(stim_alignedTraces(:,stimIdx,:),2));

    %%%% compute perimovement activity

    % align traces to movement onset
    event = 'firstMoveTimes';
    mov_alignedTraces = neuralData.eta.alignedResps{strcmp(neuralData.eta.events,event)};
    mov_eventWindow = neuralData.eta.eventWindow;

    %designate a movement window
    mov_eventIdx = find(mov_eventWindow == 0);
    movTime = [-0.2 0.1] / Fs;
    movIdx = mov_eventIdx + movTime(1) : mov_eventIdx + movTime(2);

    %compute the mean perimovement activity per cell, per trial (trials x neurons)
    movResps = squeeze(mean(mov_alignedTraces(:,movIdx,:),2));

    %%%%%%%%%%%%%%%% compute premovement activity

    %designate a movement window
    pmov_eventIdx = find(mov_eventWindow == 0);
    pmovTime = [-0.7 -0.2] / Fs;
    pmovIdx = pmov_eventIdx + pmovTime(1) : pmov_eventIdx + pmovTime(2);

    %compute the mean perimovement activity per cell, per trial (trials x neurons)
    pmovResps = squeeze(mean(mov_alignedTraces(:,pmovIdx,:),2));

    %%%%%%%%%%%%%%%% compute perireward activity

    % align traces to movement onset
    event = 'feedbackTimes';
    rew_alignedTraces = neuralData.eta.alignedResps{strcmp(neuralData.eta.events,event)};
    rew_eventWindow = neuralData.eta.eventWindow;

    %designate a movement window
    rew_eventIdx = find(rew_eventWindow == 0);
    rewTime = [0 0.2] / Fs;
    rewIdx = rew_eventIdx + rewTime(1) : rew_eventIdx + rewTime(2);

    %compute the mean perireward activity per cell, per trial (trials x neurons)
    rewResps = squeeze(mean(rew_alignedTraces(:,rewIdx,:),2));

    %%%%%%%%%%%%%%%% organize trial types for subsequent comparison

    contrasts = getUniqueContrasts(expInfo);
    et = behavioralData;
    
    %stimulus
    [~, stimTrials] = selectCondition(expInfo, contrasts, et, initTrialConditions('movementTime','late'));
    [~, leftstimTrials] = selectCondition(expInfo, contrasts(contrasts<0), et, initTrialConditions('movementTime','late'));
    [~, rightstimTrials] = selectCondition(expInfo, contrasts(contrasts>0), et, initTrialConditions('movementTime','late'));
    [~, zerostimTrials] = selectCondition(expInfo, 0, et, initTrialConditions('movementTime','late'));

    %movement
    [~, movleftTrials] = selectCondition(expInfo, contrasts, et, initTrialConditions('movementDir','cw'));
    [~, movrightTrials] = selectCondition(expInfo, contrasts, et, initTrialConditions('movementDir','ccw'));

    %reward
    [~, correctTrials] = selectCondition(expInfo, contrasts, et, initTrialConditions('responseType','correct'));
    [~, incorrectTrials] = selectCondition(expInfo, contrasts, et, initTrialConditions('responseType','incorrect'));

    %value
    [~, highleftTrials] = selectCondition(expInfo, contrasts, et, initTrialConditions('movementTime','late','highRewardSide','left'));
    [~, highrightTrials] = selectCondition(expInfo, contrasts, et, initTrialConditions('movementTime','late','highRewardSide','right'));


    %%%%%%%%%%%%%%%% Wilcoxon tests

    labels = {'stim', 'leftStim', 'rightStim', 'mov', 'leftMov', 'rightMov', 'reward', 'advanceMov'};
    pValues = zeros(size(baselineResps,2),length(labels));
    for iCell = 1:size(baselineResps,2)
        pValues(iCell,1) = signrank(stimResps(stimTrials,iCell),baselineResps(stimTrials,iCell),'tail','right');
        pValues(iCell,2) = ranksum(stimResps(leftstimTrials,iCell),stimResps(zerostimTrials,iCell),'tail','right');
        pValues(iCell,3) = ranksum(stimResps(rightstimTrials,iCell),stimResps(zerostimTrials,iCell),'tail','right');
        pValues(iCell,4) = signrank(baselineResps(:,iCell),movResps(:,iCell));
        pValues(iCell,5) = ranksum(movResps(movleftTrials,iCell),movResps(movrightTrials,iCell),'tail','right');
        pValues(iCell,6) = ranksum(movResps(movrightTrials,iCell),movResps(movleftTrials,iCell),'tail','right');
        pValues(iCell,7) = ranksum(rewResps(correctTrials,iCell),rewResps(incorrectTrials,iCell));
        pValues(iCell,8) = ranksum(pmovResps(movleftTrials,iCell),pmovResps(movrightTrials,iCell)); 
    end
    
    %%%%%%%%%%%%%%%% Bonferroni correction

    alpha = 0.05;
    bfc = size(pValues,2);
    bfcAlpha = alpha/bfc;
    bfcH = pValues < bfcAlpha;
    
    %%%%%%%%%%%%%%%% collect into struct
    neuralData.stats = struct(...
        'pValues',pValues,...
        'labels',{labels},...
        'alpha',alpha,...
        'bfcAlpha',bfcAlpha,...
        'bfcH',bfcH...
        );
    
end
end