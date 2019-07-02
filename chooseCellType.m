function plotAll = chooseCellType(propType, mouseName, expDate, expNum, expSeries, block, allFcell, eventTimes, ops)
%% select cells with the properties you want
% propType = 'all'; % 'all' or 'vis' or 'mov'

clear allCells;
clear visCells;
clear movCells;

contrasts = unique(block.events.contrastValues);

switch propType
    case 'all'
        
        for iPlane = 1:ops.numPlanes
            v = [];
            for iCell = 1:size(allFcell(iPlane).spikes{1,1},1)
                v = [v; iCell];
                allCells{iPlane} = v;
            end
        end
        
        plotAll = [];
        for iPlane = 1:ops.numPlanes
            vv = [ones(length(allCells{iPlane}),1)*iPlane allCells{iPlane}];
            plotAll = [plotAll; vv];
        end
        
    case 'vis'
        % align traces to stim onset for testing visual responsiveness
        event = 'stimulusOnTimes';
        [stim_alignedTraces, stim_eventWindow] = getExpTraces(mouseName, expDate, expNum, expSeries, allFcell, eventTimes, ops, event);

        %imaging rate
        numPlanes = length(stim_alignedTraces);
        Fs = 15;%/ numPlanes;

        %event window
        stim_eventIdx = find(stim_eventWindow == 0);

        % stim_eventIdx = 3;
        %stimulus response index
        stim_crfTime = .2;
        stim_vPreIdx = stim_eventIdx - 3 : stim_eventIdx - 1;
        stim_vPeriIdx = stim_eventIdx : stim_eventIdx + ceil(stim_crfTime*Fs);

        % cursory test for stimulus responsiveness
        % is response to -100% contrast trials statistically significantly different from
        % the response to +100% contrast trials?

        [~, condIdx_left100] = selectCondition(block, contrasts(contrasts < 0), eventTimes, 'all', 'all', 'all', 'all', 'all', 'all', 'all', 'all', 'all');
        [~, condIdx_right100] = selectCondition(block, contrasts(contrasts > 0), eventTimes, 'all', 'all', 'all', 'all', 'all', 'all', 'all', 'all', 'all');

        for iPlane = 1:numPlanes
            v = [];
            for iCell = 1:size(stim_alignedTraces{iPlane}.eventSpikes,3)

                %compare at a timepoint after stimulus onset
                leftStimResp = mean(stim_alignedTraces{iPlane}.eventSpikes(condIdx_left100,stim_vPeriIdx,iCell),2);
                rightStimResp = mean(stim_alignedTraces{iPlane}.eventSpikes(condIdx_right100,stim_vPeriIdx,iCell),2);
                if ~isnan(mean(leftStimResp)) && ~isnan(mean(rightStimResp))
                    [~,p] = kstest2(leftStimResp,rightStimResp);
                end

                % if it a. passes the statistical test and b. left-stim responses
                % are larger than right-stim responses (classical CRF shape)
                if  p < 0.01 && (mean(leftStimResp) > mean(rightStimResp))
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
        [mov_alignedTraces, mov_eventWindow] = getExpTraces(mouseName, expDate, expNum, expSeries, allFcell, eventTimes, ops, event);

        %imaging rate
        numPlanes = length(mov_alignedTraces);
        Fs = 15;%/ numPlanes;

        %event window
        mov_eventIdx = find(mov_eventWindow == 0);

        % stim_eventIdx = 3;
        %stimulus response index
        mov_crfTime = .2;
        mov_vPreIdx = mov_eventIdx - 2 : mov_eventIdx ;
        mov_vPeriIdx = mov_eventIdx + 1: mov_eventIdx + 3;

        % cursory test for movement responsiveness
        % is response to -100% contrast trials statistically significantly different from
        % the response to +100% contrast trials?

        [~, condIdx_mov] = selectCondition(block, contrasts, eventTimes, 'all', 'all', 'all', 'all', 'all',  'all', 'all', 'all', 'all');

        for iPlane = 1:numPlanes
            m = [];
            for iCell = 1:size(mov_alignedTraces{iPlane}.eventSpikes,3)

                %compare at a timepoint after stimulus onset
                preMovResp = nanmean(mov_alignedTraces{iPlane}.eventSpikes(condIdx_mov,mov_vPreIdx,iCell),2);
                periMovResp = nanmean(mov_alignedTraces{iPlane}.eventSpikes(condIdx_mov,mov_vPeriIdx,iCell),2);
                if ~isnan(nanmean(preMovResp)) && ~isnan(nanmean(periMovResp))
                    [~,p] = kstest2(preMovResp,periMovResp);
                end

                % if it a. passes the statistical test and b. left-stim responses
                % are larger than right-stim responses (classical CRF shape)
                if  p < 0.05 && (nanmean(preMovResp) < nanmean(periMovResp))
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
end