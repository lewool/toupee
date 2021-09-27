function neuralData = getSignificantActivity(expInfo, behavioralData, neuralData, matched)
% Determine which neurons are active for different aspects of the task:
% stimulus onset, movement onset, movement direction, reward onset, block
% identity (ITI)

Fs = 0.1;
    
if matched == 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%% UNMATCHED STATS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    for ex = 1:length(expInfo)
        
    %%%%%%%%%%%%%%%% compute baseline activity
    et = behavioralData(ex);
    
    if isstruct(et.eventTimes) && isstruct(et.wheelMoves)
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

        %%%%%%%%%%%%%%%% organize trial types for subsequent comparison
        try
            contrasts = unique(expInfo(ex).block.events.contrastValues);
            leftScreens = find(sign(expInfo(ex).block.events.azimuthValues) == -1);
            centerScreens = find(sign(expInfo(ex).block.events.azimuthValues) == 0);
            rightScreens = find(sign(expInfo(ex).block.events.azimuthValues) == 1);

            %stimulus
            [~, stimTrials] = selectCondition(expInfo(ex), contrasts, et, initTrialConditions());
    %         [~, leftstimTrials] = selectCondition(expInfo(ex), contrasts(contrasts<0), et, initTrialConditions());
    %         [~, rightstimTrials] = selectCondition(expInfo(ex), contrasts(contrasts>0), et, initTrialConditions());
    %         [~, zerostimTrials] = selectCondition(expInfo(ex), 0, et, initTrialConditions());

            %%%%%%%%%%%%%%%% Wilcoxon tests

            labels = {'leftStim', 'centerStim', 'rightStim'};
            pValues = zeros(size(baselineResps,2),length(labels));
            for iCell = 1:size(baselineResps,2)
                pValues(iCell,1) = signrank(stimResps(leftScreens,iCell),baselineResps(leftScreens,iCell),'tail','right');
                pValues(iCell,2) = signrank(stimResps(centerScreens,iCell),stimResps(centerScreens,iCell),'tail','right');
                pValues(iCell,3) = signrank(stimResps(rightScreens,iCell),stimResps(rightScreens,iCell),'tail','right');
            end
        
        catch
            nt = numel(expInfo(ex).block.events.endTrialValues);
            
            ori0 = find(expInfo(ex).block.events.oriValues(1:nt) == 0);
            ori45 = find(expInfo(ex).block.events.oriValues(1:nt) == 45);
            ori90 = find(expInfo(ex).block.events.oriValues(1:nt) == 90);
            ori135 = find(expInfo(ex).block.events.oriValues(1:nt) == 135);
            ori180 = find(expInfo(ex).block.events.oriValues(1:nt) == 180);
            ori225 = find(expInfo(ex).block.events.oriValues(1:nt) == 225);
            ori270 = find(expInfo(ex).block.events.oriValues(1:nt) == 270);
            ori315 = find(expInfo(ex).block.events.oriValues(1:nt) == 315);
            ori360 = find(expInfo(ex).block.events.oriValues(1:nt) == 360);
            
            labels = {'ori0', 'ori45', 'ori90', 'ori135', 'ori180', 'ori225', 'ori270', 'ori315', 'ori360'};
            pValues = zeros(size(baselineResps,2),length(labels));
            for iCell = 1:size(baselineResps,2)
                pValues(iCell,1) = signrank(stimResps(ori0,iCell),baselineResps(ori0,iCell),'tail','right');
                pValues(iCell,2) = signrank(stimResps(ori45,iCell),baselineResps(ori45,iCell),'tail','right');
                pValues(iCell,3) = signrank(stimResps(ori90,iCell),baselineResps(ori90,iCell),'tail','right');
                pValues(iCell,4) = signrank(stimResps(ori135,iCell),baselineResps(ori135,iCell),'tail','right');
                pValues(iCell,5) = signrank(stimResps(ori180,iCell),baselineResps(ori180,iCell),'tail','right');
                pValues(iCell,6) = signrank(stimResps(ori225,iCell),baselineResps(ori225,iCell),'tail','right');
                pValues(iCell,7) = signrank(stimResps(ori270,iCell),baselineResps(ori270,iCell),'tail','right');
                pValues(iCell,8) = signrank(stimResps(ori315,iCell),baselineResps(ori315,iCell),'tail','right');
                pValues(iCell,9) = signrank(stimResps(ori360,iCell),baselineResps(ori360,iCell),'tail','right');
            end
            
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
    else
        neuralData(ex).stats = NaN;
    end
    
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

    %%%%%%%%%%%%%%%% organize trial types for subsequent comparison

    contrasts = getUniqueContrasts(expInfo);
    et = behavioralData;
    
    %stimulus
    [~, stimTrials] = selectCondition(expInfo, contrasts, et, initTrialConditions());
    [~, leftstimTrials] = selectCondition(expInfo, contrasts(contrasts<0), et, initTrialConditions());
    [~, rightstimTrials] = selectCondition(expInfo, contrasts(contrasts>0), et, initTrialConditions());
    [~, zerostimTrials] = selectCondition(expInfo, 0, et, initTrialConditions());

    %%%%%%%%%%%%%%%% Wilcoxon tests

    labels = {'stim', 'leftStim', 'rightStim'};
    pValues = zeros(size(baselineResps,2),length(labels));
    for iCell = 1:size(baselineResps,2)
        pValues(iCell,1) = signrank(stimResps(stimTrials,iCell),baselineResps(stimTrials,iCell),'tail','right');
        pValues(iCell,2) = ranksum(stimResps(leftstimTrials,iCell),stimResps(zerostimTrials,iCell),'tail','right');
        pValues(iCell,3) = ranksum(stimResps(rightstimTrials,iCell),stimResps(zerostimTrials,iCell),'tail','right');
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