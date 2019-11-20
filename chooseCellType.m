function plotAll = chooseCellType(propType, expInfo, allFcell, eventTimes)
% Chooses the types of cells you want to analyze based on their responses 
% propType = 'all', 'vis', 'mov', 'movleft', or 'movright'

% 20 Nov 2019 Takes new trialCondition struct for selectCondition.m
%% LOAD DATA FROM EXPINFO

mouseName = expInfo.mouseName;
expDate = expInfo.expDate;
expNum = expInfo.expNum;
expSeries = expInfo.expSeries;
block = expInfo.block;
Timeline = expInfo.Timeline;
numPlanes = expInfo.numPlanes;

%% select cells with the properties you want

clear allCells;
clear visCells;
clear movCells;

propContrasts = unique(block.events.contrastValues);

switch propType
    case 'all'
        
        for iPlane = 1:numPlanes
            v = [];
            for iCell = 1:size(allFcell(iPlane).spikes{1,1},1)
                v = [v; iCell];
                allCells{iPlane} = v;
            end
        end
        
        plotAll = [];
        for iPlane = 1:numPlanes
            vv = [ones(length(allCells{iPlane}),1)*iPlane allCells{iPlane}];
            plotAll = [plotAll; vv];
        end
        
    case 'vis'
        % align traces to stim onset for testing visual responsiveness
        event = 'stimulusOnTimes';
        [stim_alignedTraces, stim_eventWindow] = getAlignedTraces(expInfo, allFcell, eventTimes, event);

        %imaging rate
        numPlanes = length(stim_alignedTraces);
        Fs = 15;%/ numPlanes;

        %event window
        stim_eventIdx = find(stim_eventWindow == 0);

        % stim_eventIdx = 3;
        %stimulus response index
        stim_crfTime = .5;
        stim_vPreIdx = stim_eventIdx - 3 : stim_eventIdx - 1;
        stim_vPeriIdx = stim_eventIdx : stim_eventIdx + ceil(stim_crfTime*Fs);
        sdTestIdx = stim_eventIdx : stim_eventIdx + 6;

        % cursory test for stimulus responsiveness
        % is response to left (or right) contrast trials statistically significantly different from
        % the response to 0% contrast trials?
        
        trialConditions_late = initTrialConditions('movementTime','late');
        [~, condIdx_left] = selectCondition(block, (propContrasts(propContrasts < 0)), eventTimes, trialConditions_late);
        [~, condIdx_right] = selectCondition(block, (propContrasts(propContrasts > 0)), eventTimes, trialConditions_late);
        [~, condIdx_zero] = selectCondition(block, 0, eventTimes, trialConditions_late);

        for iPlane = 1:numPlanes
            v = [];
            for iCell = 1:size(stim_alignedTraces{iPlane}.eventSpikes,3)

                %compare at a timepoint after stimulus onset
                leftStimResp = mean(stim_alignedTraces{iPlane}.eventSpikes(condIdx_left,stim_vPeriIdx,iCell),2);
                rightStimResp = mean(stim_alignedTraces{iPlane}.eventSpikes(condIdx_right,stim_vPeriIdx,iCell),2);
                zeroStimResp = mean(stim_alignedTraces{iPlane}.eventSpikes(condIdx_zero,stim_vPeriIdx,iCell),2);
                if ~isnan(mean(leftStimResp)) && ~isnan(mean(zeroStimResp)) && ~isnan(mean(rightStimResp))
                    [~,pL] = ttest2(leftStimResp,zeroStimResp);
                    [~,pR] = ttest2(rightStimResp,zeroStimResp);
                else
                    pL = 1;
                    pR = 1;
                end
                
                %check if the response is >0.25 SD in the perievent window
                sdCheck = sum(mean(stim_alignedTraces{iPlane}.eventSpikes(condIdx_left,sdTestIdx,iCell),1) > .25) > 1;

                % if it a. passes the statistical test and b. left-stim responses
                % are larger than zero-stim responses (classical CRF
                % shape) and c. sdCheck = true
                if  (pL < 0.01 && sdCheck) || (pR < 0.01 && sdCheck)
                    v = [v; iCell];
                end
                visCells{iPlane} = v;
            end
        end

        plotAll = [];
        for iPlane = 1:numPlanes
            vv = [ones(length(visCells{iPlane}),1)*iPlane visCells{iPlane}];
            plotAll = [plotAll; vv];
        end
        
    case 'mov'

        % align traces to stim onset for testing movement responsiveness
        event = 'prestimulusQuiescenceEndTimes';
        [mov_alignedTraces, mov_eventWindow] = getAlignedTraces(expInfo, allFcell, eventTimes, event);

        %imaging rate
        numPlanes = length(mov_alignedTraces);
        Fs = 15;%/ numPlanes;

        %event window
        mov_eventIdx = find(mov_eventWindow == 0);

        % stim_eventIdx = 3;
        %stimulus response index
        mov_window = 5;
        mov_vPreIdx = mov_eventIdx - mov_window : mov_eventIdx - 1;
        mov_vPeriIdx = mov_eventIdx + 1: mov_eventIdx + mov_window;
        sdTestIdx = mov_eventIdx : mov_eventIdx + 6;

        % cursory test for movement responsiveness
        % is response during 0% contrast trials significantly different
        % between pre- and post-movement onset?
        
        trialConditions_default = initTrialConditions;
        [~, condIdx_mov] = selectCondition(block, 0, eventTimes, trialConditions_default);

        for iPlane = 1:numPlanes
            m = [];
            for iCell = 1:size(mov_alignedTraces{iPlane}.eventSpikes,3)

                %compare at a timepoint after movement onset
                preMovResp = nanmean(mov_alignedTraces{iPlane}.eventSpikes(condIdx_mov,mov_vPreIdx,iCell),2);
                periMovResp = nanmean(mov_alignedTraces{iPlane}.eventSpikes(condIdx_mov,mov_vPeriIdx,iCell),2);
                if ~isnan(nanmean(preMovResp)) && ~isnan(nanmean(periMovResp))
                    [~,p] = ttest2(preMovResp,periMovResp);
                end
                
                %check if the response is >0.25 SD in the perievent window
                sdCheck = sum(mean(mov_alignedTraces{iPlane}.eventSpikes(condIdx_mov,sdTestIdx,iCell),1) > 0.25) > 1;

                % if it a. passes the statistical test and b. post-move responses
                % are larger than pre-move responses (classical CRF shape)
                % and c. sdCheck = true
                if  p < 0.01 && (nanmean(preMovResp) < nanmean(periMovResp)) && sdCheck
                    m = [m; iCell];
                end
                movCells{iPlane} = m;
            end
        end

        plotAll = [];
        for iPlane = 1:numPlanes
            vv = [ones(length(movCells{iPlane}),1)*iPlane movCells{iPlane}];
            plotAll = [plotAll; vv];
        end
        
    case 'movleft'

        % align traces to stim onset for testing movement responsiveness
        event = 'prestimulusQuiescenceEndTimes';
        [mov_alignedTraces, mov_eventWindow] = getAlignedTraces(expInfo, allFcell, eventTimes, event);

        %imaging rate
        numPlanes = length(mov_alignedTraces);
        Fs = 15;%/ numPlanes;

        %event window
        mov_eventIdx = find(mov_eventWindow == 0);

        % stim_eventIdx = 3;
        %stimulus response index
        mov_window = 5;
        mov_vPreIdx = mov_eventIdx - mov_window : mov_eventIdx - 1;
        mov_vPeriIdx = mov_eventIdx + 1: mov_eventIdx + mov_window;
        sdTestIdx = mov_eventIdx : mov_eventIdx + 6;

        % cursory test for movement responsiveness
        % is response during 0% contrast trials significantly different
        % between pre- and post-movement onset?
        
        trialConditions_movLeft = initTrialConditions('movementTime','late','movementDir','cw');
        trialConditions_movRight = initTrialConditions('movementTime','late','movementDir','ccw');
        [~, condIdx_movleft] = selectCondition(block, propContrasts, eventTimes, trialConditions_movLeft);
        [~, condIdx_movright] = selectCondition(block, propContrasts, eventTimes, trialConditions_movRight);

        for iPlane = 1:numPlanes
            v = [];
            for iCell = 1:size(mov_alignedTraces{iPlane}.eventSpikes,3)

                %compare at a timepoint after stimulus onset
                leftMovResp = mean(mov_alignedTraces{iPlane}.eventSpikes(condIdx_movleft,mov_vPeriIdx,iCell),2);
                leftPreResp =  mean(mov_alignedTraces{iPlane}.eventSpikes(condIdx_movleft,mov_vPreIdx,iCell),2);
                rightMovResp = mean(mov_alignedTraces{iPlane}.eventSpikes(condIdx_movright,mov_vPeriIdx,iCell),2);
                rightPreResp =  mean(mov_alignedTraces{iPlane}.eventSpikes(condIdx_movright,mov_vPreIdx,iCell),2);
                if ~isnan(mean(leftMovResp))
                    [~,p] = ttest2([leftMovResp; rightMovResp],[leftPreResp; rightPreResp]);
                    [~,pL] = ttest2(leftMovResp,rightMovResp);
                else
                    p = 1;
                    pL = 1;
                end
                
                %check if the response is >0.25 SD in the perievent window
                sdCheck = sum(mean(mov_alignedTraces{iPlane}.eventSpikes([condIdx_movright condIdx_movleft],sdTestIdx,iCell),1) > .25) > 1;

                % if it a. passes the statistical test and b. left-stim responses
                % are larger than zero-stim responses (classical CRF
                % shape) and c. sdCheck = true
                if  (p < 0.01 && pL < 0.01 && sdCheck && mean(rightMovResp) < mean(leftMovResp))
                    v = [v; iCell];
                end
                movCells{iPlane} = v;
            end
        end
        plotAll = [];
        for iPlane = 1:numPlanes
            vv = [ones(length(movCells{iPlane}),1)*iPlane movCells{iPlane}];
            plotAll = [plotAll; vv];
        end
    
        case 'movright'

        % align traces to stim onset for testing movement responsiveness
        event = 'prestimulusQuiescenceEndTimes';
        [mov_alignedTraces, mov_eventWindow] = getAlignedTraces(expInfo, allFcell, eventTimes, event);

        %imaging rate
        numPlanes = length(mov_alignedTraces);
        Fs = 15;%/ numPlanes;

        %event window
        mov_eventIdx = find(mov_eventWindow == 0);

        % stim_eventIdx = 3;
        %stimulus response index
        mov_window = 5;
        mov_vPreIdx = mov_eventIdx - mov_window : mov_eventIdx - 1;
        mov_vPeriIdx = mov_eventIdx + 1: mov_eventIdx + mov_window;
        sdTestIdx = mov_eventIdx : mov_eventIdx + 6;

        % cursory test for movement responsiveness
        % is response during 0% contrast trials significantly different
        % between pre- and post-movement onset?

        trialConditions_movLeft = initTrialConditions('movementTime','late','movementDir','cw');
        trialConditions_movRight = initTrialConditions('movementTime','late','movementDir','ccw');
        [~, condIdx_movleft] = selectCondition(block, propContrasts, eventTimes, trialConditions_movLeft);
        [~, condIdx_movright] = selectCondition(block, propContrasts, eventTimes, trialConditions_movRight);
        
        for iPlane = 1:numPlanes
            v = [];
            for iCell = 1:size(mov_alignedTraces{iPlane}.eventSpikes,3)

                %compare at a timepoint after stimulus onset
                leftMovResp = mean(mov_alignedTraces{iPlane}.eventSpikes(condIdx_movleft,mov_vPeriIdx,iCell),2);
                leftPreResp =  mean(mov_alignedTraces{iPlane}.eventSpikes(condIdx_movleft,mov_vPreIdx,iCell),2);
                rightMovResp = mean(mov_alignedTraces{iPlane}.eventSpikes(condIdx_movright,mov_vPeriIdx,iCell),2);
                rightPreResp =  mean(mov_alignedTraces{iPlane}.eventSpikes(condIdx_movright,mov_vPreIdx,iCell),2);
                if ~isnan(mean(rightMovResp))
                    [~,p] = ttest2([leftMovResp; rightMovResp],[leftPreResp; rightPreResp]);
                    [~,pR] = ttest2(leftMovResp,rightMovResp);
                else
                    p = 1;
                    pR = 1;
                end
                
                %check if the response is >0.25 SD in the perievent window
                sdCheck = sum(mean(mov_alignedTraces{iPlane}.eventSpikes([condIdx_movright condIdx_movleft],sdTestIdx,iCell),1) > .25) > 1;

                % if it a. passes the statistical test and b. left-stim responses
                % are larger than zero-stim responses (classical CRF
                % shape) and c. sdCheck = true
                if  (p < 0.01 && pR < 0.01 && sdCheck && mean(rightMovResp) > mean(leftMovResp))
                    v = [v; iCell];
                end
                movCells{iPlane} = v;
            end
        end

        plotAll = [];
        for iPlane = 1:numPlanes
            vv = [ones(length(movCells{iPlane}),1)*iPlane movCells{iPlane}];
            plotAll = [plotAll; vv];
        end
end

% remove the flyback plane cells from analyses
plotAll(plotAll(:,1)==1,:)=[];
