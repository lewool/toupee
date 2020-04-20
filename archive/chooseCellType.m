function plotCells = chooseCellType(propType, expInfo, cellResps, respTimes, eventTimes, Fs)
% Chooses the types of cells you want to analyze based on their responses 
% propType = 'all', 'vis', 'mov', 'movleft', or 'movright'

% 20 Nov 2019 Takes new trialCondition struct for selectCondition.m
% 14 Jan 2020 Takes new cellResps array (all planes combined)
% 22 Jan 2020 Test vis responses from passive experiments
%% LOAD DATA FROM EXPINFO

block = expInfo.block;

%% select cells with the properties you want

clear allCells;
clear visCells;
clear movCells;

propContrasts = unique(block.events.contrastValues);
plotCells = [];

switch propType
    case 'all'
        
        plotCells = 1:size(cellResps,2);
        
    case 'vis_passive'
        
        % align traces to stim onset for testing visual responsiveness
        event = 'stimulusOnTimes';
        [stim_alignedTraces, stim_eventWindow] = alignResps(expInfo, cellResps, respTimes, eventTimes, event);

        %event window
        stim_eventIdx = find(stim_eventWindow == 0);

        %stimulus response index
        stim_periTime = [0 0.3] / Fs;
        stim_preTime = [-0.3 0] / Fs;
        stim_vPreIdx = stim_eventIdx + stim_preTime(1) : stim_eventIdx;
        stim_vPeriIdx = stim_eventIdx + stim_periTime(1) :stim_eventIdx + stim_periTime(2);
        sdTestIdx = stim_eventIdx : stim_eventIdx + 6;
        
        trialConditions = initTrialConditions;
        [~, condIdx] = selectCondition(block, propContrasts(propContrasts > .12), eventTimes, trialConditions);
        [~, condIdx_zero] = selectCondition(block, 0, eventTimes, trialConditions);
        leftScreenIdx = find(block.events.azimuthValues == -90);
        rightScreenIdx = find(block.events.azimuthValues == 90);
        condIdx_left = intersect(condIdx,leftScreenIdx);
        condIdx_right = intersect(condIdx,rightScreenIdx);
        
        for iCell = 1:size(stim_alignedTraces,3)

            %compare at a timepoint after stimulus onset
            leftStimResp = mean(stim_alignedTraces(condIdx_left,stim_vPeriIdx,iCell),2);
            rightStimResp = mean(stim_alignedTraces(condIdx_right,stim_vPeriIdx,iCell),2);
            zeroStimResp = mean(stim_alignedTraces(condIdx_zero,stim_vPeriIdx,iCell),2);
            if ~isnan(mean(leftStimResp)) && ~isnan(mean(zeroStimResp)) && ~isnan(mean(rightStimResp))
                [~,pL] = ttest2(leftStimResp,zeroStimResp);
                [~,pR] = ttest2(rightStimResp,zeroStimResp);
            else
                pL = 1;
                pR = 1;
            end

            %check if the response is >0.25 SD in the perievent window
            sdCheck = sum(mean(stim_alignedTraces(condIdx_left,sdTestIdx,iCell),1) > .25) > 1;

            % if it passes the statistical test and sdCheck = true
            if  (pL < 0.01) || (pR < 0.01)
                plotCells = [plotCells; iCell];
            end
        end
        
    case 'vis'
        % align traces to stim onset for testing visual responsiveness
        event = 'stimulusOnTimes';
        [stim_alignedTraces, stim_eventWindow] = alignResps(expInfo, cellResps, respTimes, eventTimes, event);

        %event window
        stim_eventIdx = find(stim_eventWindow == 0);

        %stimulus response index
        stim_periTime = [0 0.5] / Fs;
        stim_preTime = [-0.5 0] / Fs;
        stim_vPreIdx = stim_eventIdx + stim_preTime(1) : stim_eventIdx;
        stim_vPeriIdx = stim_eventIdx + stim_periTime(1) :stim_eventIdx + stim_periTime(2);
        sdTestIdx = stim_eventIdx : stim_eventIdx + 6;

        % cursory test for stimulus responsiveness
        % is response to left (or right) contrast trials statistically significantly different from
        % the response to 0% contrast trials?
        
        trialConditions_late = initTrialConditions('movementTime','late');
        [~, condIdx_left] = selectCondition(block, (propContrasts(propContrasts < 0)), eventTimes, trialConditions_late);
        [~, condIdx_right] = selectCondition(block, (propContrasts(propContrasts > 0)), eventTimes, trialConditions_late);
        [~, condIdx_zero] = selectCondition(block, 0, eventTimes, trialConditions_late);

        
        for iCell = 1:size(stim_alignedTraces,3)

            %compare at a timepoint after stimulus onset
            leftStimResp = mean(stim_alignedTraces(condIdx_left,stim_vPeriIdx,iCell),2);
            rightStimResp = mean(stim_alignedTraces(condIdx_right,stim_vPeriIdx,iCell),2);
            zeroStimResp = mean(stim_alignedTraces(condIdx_zero,stim_vPeriIdx,iCell),2);
            if ~isnan(mean(leftStimResp)) && ~isnan(mean(zeroStimResp)) && ~isnan(mean(rightStimResp))
                [~,pL] = ttest2(leftStimResp,zeroStimResp);
                [~,pR] = ttest2(rightStimResp,zeroStimResp);
            else
                pL = 1;
                pR = 1;
            end

            %check if the response is >0.25 SD in the perievent window
            sdCheck = sum(mean(stim_alignedTraces(condIdx_left,sdTestIdx,iCell),1) > .25) > 1;

            % if it passes the statistical test and sdCheck = true
            if  (pL < 0.01 ) || (pR < 0.01 )
                plotCells = [plotCells; iCell];
            end
        end
        
    case 'mov'

        % align traces to stim onset for testing movement responsiveness
        event = 'prestimulusQuiescenceEndTimes';
        [mov_alignedTraces, mov_eventWindow] = alignResps(expInfo, cellResps, respTimes, eventTimes, event);

        %event window
        mov_eventIdx = find(mov_eventWindow == 0);

        %movement response index
        mov_preTime= [-0.5 -0.2] / Fs;
        mov_periTime = [-0.2 0.3] / Fs;
        mov_preIdx = mov_eventIdx + mov_preTime(1) : mov_eventIdx + mov_preTime(2);
        mov_periIdx = mov_eventIdx + mov_periTime(1) : mov_eventIdx + mov_periTime(2);
        sdTestIdx = mov_eventIdx : mov_eventIdx + 6;

        % cursory test for movement responsiveness
        % is response during 0% contrast trials significantly different
        % between pre- and post-movement onset?
        
        trialConditions_default = initTrialConditions;
        [~, condIdx_mov] = selectCondition(block, 0, eventTimes, trialConditions_default);

        for iCell = 1:size(mov_alignedTraces,3)

            %compare at a timepoint after movement onset
            preMovResp = nanmean(mov_alignedTraces(condIdx_mov,mov_preIdx,iCell),2);
            periMovResp = nanmean(mov_alignedTraces(condIdx_mov,mov_periIdx,iCell),2);
            if ~isnan(nanmean(preMovResp)) && ~isnan(nanmean(periMovResp))
                [~,p] = ttest2(preMovResp,periMovResp);
            end

            %check if the response is >0.25 SD in the perievent window
            sdCheck = sum(mean(mov_alignedTraces(condIdx_mov,sdTestIdx,iCell),1) > 0.25) > 1;

            % if it a. passes the statistical test and b. post-move responses
            % are larger than pre-move responses (classical CRF shape)
            % and c. sdCheck = true
            if  p < 0.01 && (nanmean(preMovResp) < nanmean(periMovResp)) && sdCheck
                plotCells = [plotCells; iCell];
            end
        end
        
    case 'movleft'

        % align traces to stim onset for testing movement responsiveness
        event = 'prestimulusQuiescenceEndTimes';
        [mov_alignedTraces, mov_eventWindow] = alignResps(expInfo, cellResps, respTimes, eventTimes, event);

        %event window
        mov_eventIdx = find(mov_eventWindow == 0);

        %movement response index
        mov_preTime= [-0.5 -0.2] / Fs;
        mov_periTime = [-0.2 0.3] / Fs;
        mov_preIdx = mov_eventIdx + mov_preTime(1) : mov_eventIdx + mov_preTime(2);
        mov_periIdx = mov_eventIdx + mov_periTime(1) : mov_eventIdx + mov_periTime(2);
        sdTestIdx = mov_eventIdx : mov_eventIdx + 6;

        % cursory test for movement responsiveness
        % is response during 0% contrast trials significantly different
        % between pre- and post-movement onset?
        
        trialConditions_movLeft = initTrialConditions('movementTime','late','movementDir','cw');
        trialConditions_movRight = initTrialConditions('movementTime','late','movementDir','ccw');
        [~, condIdx_movleft] = selectCondition(block, propContrasts, eventTimes, trialConditions_movLeft);
        [~, condIdx_movright] = selectCondition(block, propContrasts, eventTimes, trialConditions_movRight);

        for iCell = 1:size(mov_alignedTraces,3)

            %compare at a timepoint after stimulus onset
            leftMovResp = mean(mov_alignedTraces(condIdx_movleft,mov_periIdx,iCell),2);
            leftPreResp =  mean(mov_alignedTraces(condIdx_movleft,mov_preIdx,iCell),2);
            rightMovResp = mean(mov_alignedTraces(condIdx_movright,mov_periIdx,iCell),2);
            rightPreResp =  mean(mov_alignedTraces(condIdx_movright,mov_preIdx,iCell),2);
            if ~isnan(mean(leftMovResp))
                [~,p] = ttest2([leftMovResp; rightMovResp],[leftPreResp; rightPreResp]);
                [~,pL] = ttest2(leftMovResp,rightMovResp);
            else
                p = 1;
                pL = 1;
            end

            %check if the response is >0.25 SD in the perievent window
            sdCheck = sum(mean(mov_alignedTraces([condIdx_movright condIdx_movleft],sdTestIdx,iCell),1) > .25) > 1;

            % if it a. passes the statistical test and b. left-stim responses
            % are larger than zero-stim responses (classical CRF
            % shape) and c. sdCheck = true
            if  (p < 0.01 && pL < 0.01 && mean(rightMovResp) < mean(leftMovResp))
                plotCells = [plotCells; iCell];
            end
        end
    
    case 'movright'

        % align traces to stim onset for testing movement responsiveness
        event = 'prestimulusQuiescenceEndTimes';
        [mov_alignedTraces, mov_eventWindow] = alignResps(expInfo, cellResps, respTimes, eventTimes, event);
        
        %event window
        mov_eventIdx = find(mov_eventWindow == 0);

        %movement response index
        mov_preTime= [-0.5 -0.2] / Fs;
        mov_periTime = [-0.2 0.3] / Fs;
        mov_preIdx = mov_eventIdx + mov_preTime(1) : mov_eventIdx + mov_preTime(2);
        mov_periIdx = mov_eventIdx + mov_periTime(1) : mov_eventIdx + mov_periTime(2);
        sdTestIdx = mov_eventIdx : mov_eventIdx + 6;

        % cursory test for movement responsiveness
        % is response during 0% contrast trials significantly different
        % between pre- and post-movement onset?
        trialConditions_movLeft = initTrialConditions('movementTime','late','movementDir','cw');
        trialConditions_movRight = initTrialConditions('movementTime','late','movementDir','ccw');
        [~, condIdx_movleft] = selectCondition(block, propContrasts, eventTimes, trialConditions_movLeft);
        [~, condIdx_movright] = selectCondition(block, propContrasts, eventTimes, trialConditions_movRight);
        
        for iCell = 1:size(mov_alignedTraces,3)

            %compare at a timepoint after stimulus onset
            leftMovResp = mean(mov_alignedTraces(condIdx_movleft,mov_periIdx,iCell),2);
            leftPreResp =  mean(mov_alignedTraces(condIdx_movleft,mov_preIdx,iCell),2);
            rightMovResp = mean(mov_alignedTraces(condIdx_movright,mov_periIdx,iCell),2);
            rightPreResp =  mean(mov_alignedTraces(condIdx_movright,mov_preIdx,iCell),2);
            if ~isnan(mean(rightMovResp))
                [~,p] = ttest2([leftMovResp; rightMovResp],[leftPreResp; rightPreResp]);
                [~,pR] = ttest2(leftMovResp,rightMovResp);
            else
                p = 1;
                pR = 1;
            end

            %check if the response is >0.25 SD in the perievent window
            sdCheck = sum(mean(mov_alignedTraces([condIdx_movright condIdx_movleft],sdTestIdx,iCell),1) > .25) > 1;

            % if it a. passes the statistical test and b. left-stim responses
            % are larger than zero-stim responses (classical CRF
            % shape) and c. sdCheck = true
            if  (p < 0.01 && pR < 0.01 && mean(rightMovResp) > mean(leftMovResp))
                plotCells = [plotCells; iCell];
            end
        end
end
