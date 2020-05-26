function figName = rasterBrowser(expInfo, behavioralData, neuralData, whichCells, compareSides,k)

%% initialize experiment details

    alignedResps = neuralData.eta.alignedResps;
    eventWindow = neuralData.eta.eventWindow;
    bfcH = neuralData.stats.bfcH;
    pLabels = neuralData.stats.labels;
    et = behavioralData.eventTimes;
    wm = behavioralData.wheelMoves;


%% choose cells

% whichCells = 'leftStim'; %choose from 'pLabels' array
if strcmp(whichCells, 'all')
    plotCells = 1:size(alignedResps{1},3);
else
    plotCells = find(bfcH(:,strcmp(pLabels,whichCells)) > 0);
end

%% set up trial conditions to compare
contrasts = getUniqueContrasts(expInfo);
if strcmp(compareSides,'all')

    [~, stimConds{1}] = selectCondition(expInfo, contrasts(contrasts<0), behavioralData, initTrialConditions('movementTime','late'));
    [~, stimConds{2}] = selectCondition(expInfo, contrasts(contrasts==0), behavioralData, initTrialConditions('movementTime','late'));
    [~, stimConds{3}] = selectCondition(expInfo, contrasts(contrasts>0), behavioralData, initTrialConditions('movementTime','late'));
    allLabels{1} = {'Left' 'Zero' 'Right'};

    [~, moveConds{1}] = selectCondition(expInfo, contrasts, behavioralData, initTrialConditions('movementDir','cw','movementTime','late'));
    [~, moveConds{2}] = selectCondition(expInfo, contrasts, behavioralData, initTrialConditions('movementDir','ccw','movementTime','late'));
    allLabels{2} = {'Left' 'Right'};

    [~, rewConds{1}] = selectCondition(expInfo, contrasts, behavioralData, initTrialConditions('responseType','correct','movementTime','late'));
    [~, rewConds{2}] = selectCondition(expInfo, contrasts, behavioralData, initTrialConditions('responseType','incorrect','movementTime','late'));
    allLabels{3} = {'Correct' 'Incorrect'};

    allConds = {stimConds moveConds rewConds};
    
elseif strcmp(compareSides,'blocks')
    
    [~, stimConds{1}] = selectCondition(expInfo, contrasts(contrasts<0), behavioralData, initTrialConditions('highRewardSide','left','movementTime','late'));
    [~, stimConds{2}] = selectCondition(expInfo, contrasts(contrasts<0), behavioralData, initTrialConditions('highRewardSide','right','movementTime','late'));
    [~, stimConds{3}] = selectCondition(expInfo, contrasts(contrasts==0), behavioralData, initTrialConditions('highRewardSide','left','movementTime','late'));
    [~, stimConds{4}] = selectCondition(expInfo, contrasts(contrasts==0), behavioralData, initTrialConditions('highRewardSide','right','movementTime','late'));
    [~, stimConds{5}] = selectCondition(expInfo, contrasts(contrasts>0), behavioralData, initTrialConditions('highRewardSide','left','movementTime','late'));
    [~, stimConds{6}] = selectCondition(expInfo, contrasts(contrasts>0), behavioralData, initTrialConditions('highRewardSide','right','movementTime','late'));
    allLabels{1} = {'L' 'L' '0' '0' 'R' 'R'};

    [~, moveConds{1}] = selectCondition(expInfo, contrasts, behavioralData, initTrialConditions('highRewardSide','left','movementDir','cw','movementTime','late'));
    [~, moveConds{2}] = selectCondition(expInfo, contrasts, behavioralData, initTrialConditions('highRewardSide','right','movementDir','cw','movementTime','late'));
    [~, moveConds{3}] = selectCondition(expInfo, contrasts, behavioralData, initTrialConditions('highRewardSide','left','movementDir','ccw','movementTime','late'));
    [~, moveConds{4}] = selectCondition(expInfo, contrasts, behavioralData, initTrialConditions('highRewardSide','right','movementDir','ccw','movementTime','late'));
    allLabels{2} = {'L' 'L' 'R' 'R'};

    [~, rewConds{1}] = selectCondition(expInfo, contrasts, behavioralData, initTrialConditions('highRewardSide','left','responseType','correct','movementTime','late'));
    [~, rewConds{2}] = selectCondition(expInfo, contrasts, behavioralData, initTrialConditions('highRewardSide','right','responseType','correct','movementTime','late'));
    [~, rewConds{3}] = selectCondition(expInfo, contrasts, behavioralData, initTrialConditions('highRewardSide','left','responseType','incorrect','movementTime','late'));
    [~, rewConds{4}] = selectCondition(expInfo, contrasts, behavioralData, initTrialConditions('highRewardSide','right','responseType','incorrect','movementTime','late'));
    allLabels{3} = {'Cor' 'Cor' 'Inc' 'Inc'};

    allConds = {stimConds moveConds rewConds};
    
elseif strcmp(compareSides,'stimmove')
    
    [~, stimConds{1}] = selectCondition(expInfo, contrasts(contrasts<0), behavioralData, initTrialConditions('movementDir','cw','movementTime','late'));
    [~, stimConds{2}] = selectCondition(expInfo, contrasts(contrasts<0), behavioralData, initTrialConditions('movementDir','ccw','movementTime','late'));
    [~, stimConds{3}] = selectCondition(expInfo, contrasts(contrasts==0), behavioralData, initTrialConditions('movementDir','cw','movementTime','late'));
    [~, stimConds{4}] = selectCondition(expInfo, contrasts(contrasts==0), behavioralData, initTrialConditions('movementDir','ccw','movementTime','late'));
    [~, stimConds{5}] = selectCondition(expInfo, contrasts(contrasts>0), behavioralData, initTrialConditions('movementDir','cw','movementTime','late'));
    [~, stimConds{6}] = selectCondition(expInfo, contrasts(contrasts>0), behavioralData, initTrialConditions('movementDir','ccw','movementTime','late'));
    allLabels{1} = {'L' 'L' '0' '0' 'R' 'R'};

    [~, moveConds{1}] = selectCondition(expInfo, contrasts(contrasts<0), behavioralData, initTrialConditions('movementDir','cw','movementTime','late'));
    [~, moveConds{2}] = selectCondition(expInfo, contrasts(contrasts<0), behavioralData, initTrialConditions('movementDir','ccw','movementTime','late'));
    [~, moveConds{3}] = selectCondition(expInfo, contrasts(contrasts==0), behavioralData, initTrialConditions('movementDir','cw','movementTime','late'));
    [~, moveConds{4}] = selectCondition(expInfo, contrasts(contrasts==0), behavioralData, initTrialConditions('movementDir','ccw','movementTime','late'));
    [~, moveConds{5}] = selectCondition(expInfo, contrasts(contrasts>0), behavioralData, initTrialConditions('movementDir','cw','movementTime','late'));
    [~, moveConds{6}] = selectCondition(expInfo, contrasts(contrasts>0), behavioralData, initTrialConditions('movementDir','ccw','movementTime','late'));
    allLabels{2} = {'L' 'L' '0' '0' 'R' 'R'};

    [~, rewConds{1}] = selectCondition(expInfo, contrasts(contrasts<0), behavioralData, initTrialConditions('responseType','correct','movementTime','late'));
    [~, rewConds{2}] = selectCondition(expInfo, contrasts(contrasts>0), behavioralData, initTrialConditions('responseType','correct','movementTime','late'));
    [~, rewConds{3}] = selectCondition(expInfo, contrasts(contrasts<0), behavioralData, initTrialConditions('responseType','incorrect','movementTime','late'));
    [~, rewConds{4}] = selectCondition(expInfo, contrasts(contrasts>0), behavioralData, initTrialConditions('responseType','incorrect','movementTime','late'));
    allLabels{3} = {'Cor' 'Cor' 'Inc' 'Inc'};

    allConds = {stimConds moveConds rewConds};
    
end
%% set up some plot values

%compute the total size of the figure based on trials
% if numel(cell2mat(stimConds)) == numel(cell2mat(moveConds))
    for a = 1:length(allConds)
        for iCond = 1:length(allConds{a})
            tl{a}(iCond) = numel(allConds{a}{iCond});
        end
    end
% else
%     error('Trials got dropped...check trial conditions')
% end

buffer = 1;
total_raster = sum(tl{1}) + buffer*(length(stimConds)+1);
psth = round(total_raster*.5);
total_length = sum(tl{1})+ buffer*(length(stimConds)+1) + psth;

if strcmp(compareSides,'all')
    rasterColors = {[0 0 1; 0 0 0; 1 0 0],[0 0 1; 1 0 0],[.5 1 .5; 1 0 0]};
elseif strcmp(compareSides,'blocks')
    rasterColors = {[0 0 1; 0 0 1; 0 0 0; 0 0 0; 1 0 0; 1 0 0],[0 0 1; 0 0 1; 1 0 0; 1 0 0],[.5 1 .5; .5 1 .5; 1 0 0; 1 0 0]};
elseif strcmp(compareSides,'stimmove')
    rasterColors = {[0 0 1; 0 0 1; 0 0 0; 0 0 0; 1 0 0; 1 0 0],[0 0 1; 0 0 1; 0 0 0; 0 0 0; 1 0 0; 1 0 0],[.5 1 .5; .5 1 .5; 1 0 0; 1 0 0]};
end

%% plot (all trials)
fig = figure;
set(fig, 'Position', [680 540 1080 320]);
hold on;

if ~exist('k') == 1
    k = 1;
end

max_k = length(plotCells);

while k <= max_k

    %clear subplots
    for a = 1:length(allConds)
        for iCond = 1:length(allConds{a})
            [spidx1, spidx2] = goToRasterSubplot(length(allConds), total_length, tl{a}, a, iCond, 1, psth);
            subplot(total_length,length(allConds),[spidx1 spidx2])
            cla;
        end
        spidxA = sub2ind([length(allConds) total_length], a, 1);
        spidxB = sub2ind([length(allConds) total_length], a, psth);    
        subplot(total_length,length(allConds),[spidxA spidxB])
        cla;
    end
    
    iMin = [];
    iMax = [];
    yMin = [];
    yMax = [];
    
    for a = 1:length(allConds)
        
        %plot rasters
        for iCond = 1:length(allConds{a})

            % extract the relevant trials for that condition
            whichTrials = allConds{a}{iCond};
            numTrials = size(whichTrials,2);

            [relativeTimes, sortIdx] = sortTrialTimes(whichTrials, et, wm);

            [spidx1, spidx2] = goToRasterSubplot(length(allConds), total_length, tl{a}, a, iCond, 1, psth);
            ax = subplot(total_length,length(allConds),[spidx1 spidx2]);
            f = imagesc(eventWindow,1:numTrials,alignedResps{a}(whichTrials(sortIdx),:,plotCells(k)));
            colormap(ax,rasterColor(rasterColors{a}(iCond,:)));
            hold on;
%             line([0 0],[1 length(whichTrials)],'LineStyle','-','Color',[0 0 0]);
            mt = plot(relativeTimes(:,2,a),1:numTrials,'bo');
            st = plot(relativeTimes(:,1,a),1:numTrials,'ko');
            rt = plot(relativeTimes(:,3,a),1:numTrials,'ko');
            set([st mt rt], 'MarkerSize',1,'MarkerFaceColor','k','MarkerEdgeColor','none');
            box off
            set(gca,'ytick',[]);
            set(gca,'tickdir','out')
            ylabel(allLabels{a}(iCond))
            if iCond < length(allConds{a})
            set(gca,'xtick',[]);
            end

            iMax(end+1) = max(max(f.CData));
            iMin(end+1) = min(min(f.CData));
            
            if a > 1
                xlim([-2 1]);
            else
                xlim([-1 2]);
            end
        end
        
        %plot psth at the top
        spidxA = sub2ind([length(allConds) total_length], a, 1);
        spidxB = sub2ind([length(allConds) total_length], a, psth);    
        subplot(total_length,length(allConds),[spidxA spidxB])
        
        if strcmp(compareSides,'all')
            if a == 1
                [meanResp, semResp, plotColors] = computeContrastCurves(alignedResps{a}(:,:,plotCells(k)),expInfo,behavioralData,initTrialConditions('movementTime','late'),'side');
                plotContrastCurves(eventWindow, meanResp, semResp, plotColors,'-');
                xlim([-1 2]);
            elseif a == 2
                [meanResp, semResp, plotColors] = computeMoveDirCurves(alignedResps{a}(:,:,plotCells(k)),expInfo,behavioralData,initTrialConditions('movementTime','late'));
                plotContrastCurves(eventWindow, meanResp, semResp, plotColors,'-');
                xlim([-2 1]);
            elseif a == 3
                [meanResp, semResp, plotColors] = computeRewardCurves(alignedResps{a}(:,:,plotCells(k)),expInfo,behavioralData,initTrialConditions('movementTime','late'));
                plotContrastCurves(eventWindow, meanResp, semResp, plotColors,'-');
                xlim([-2 1]);
            end
        elseif strcmp(compareSides,'blocks')
            if a == 1
                [meanResp, semResp, plotColors] = computeContrastCurves(alignedResps{a}(:,:,plotCells(k)),expInfo,behavioralData,initTrialConditions('highRewardSide','left','movementTime','late'),'side');
                plotContrastCurves(eventWindow, meanResp, semResp, plotColors,'-');
                yMin(end+1) = min(min(min(meanResp-semResp)));
                yMax(end+1) = max([.1 max(max(max(meanResp+semResp)))]);
                [meanResp, semResp, plotColors] = computeContrastCurves(alignedResps{a}(:,:,plotCells(k)),expInfo,behavioralData,initTrialConditions('highRewardSide','right','movementTime','late'),'side');
                plotContrastCurves(eventWindow, meanResp, semResp, plotColors,':');
                xlim([-1 2]);
                yMin(end+1) = min(min(min(meanResp-semResp)));
                yMax(end+1) = max([.1 max(max(max(meanResp+semResp)))]);
            elseif a == 2
                [meanResp, semResp, plotColors] = computeMoveDirCurves(alignedResps{a}(:,:,plotCells(k)),expInfo,behavioralData,initTrialConditions('highRewardSide','left','movementTime','late'));
                plotContrastCurves(eventWindow, meanResp, semResp, plotColors,'-');
                yMin(end+1) = min(min(min(meanResp-semResp)));
                yMax(end+1) = max([.1 max(max(max(meanResp+semResp)))]);
                [meanResp, semResp, plotColors] = computeMoveDirCurves(alignedResps{a}(:,:,plotCells(k)),expInfo,behavioralData,initTrialConditions('highRewardSide','right','movementTime','late'));
                plotContrastCurves(eventWindow, meanResp, semResp, plotColors,':');
                xlim([-2 1]);
                yMin(end+1) = min(min(min(meanResp-semResp)));
                yMax(end+1) = max([.1 max(max(max(meanResp+semResp)))]);
            elseif a == 3
                [meanResp, semResp, plotColors] = computeRewardCurves(alignedResps{a}(:,:,plotCells(k)),expInfo,behavioralData,initTrialConditions('highRewardSide','left','movementTime','late'));
                plotContrastCurves(eventWindow, meanResp, semResp, plotColors,'-');
                yMin(end+1) = min(min(min(meanResp-semResp)));
                yMax(end+1) = max([.1 max(max(max(meanResp+semResp)))]);
                [meanResp, semResp, plotColors] = computeRewardCurves(alignedResps{a}(:,:,plotCells(k)),expInfo,behavioralData,initTrialConditions('highRewardSide','right','movementTime','late'));
                plotContrastCurves(eventWindow, meanResp, semResp, plotColors,':');
                xlim([-2 1]);
                yMin(end+1) = min(min(min(meanResp-semResp)));
                yMax(end+1) = max([.1 max(max(max(meanResp+semResp)))]);
            end
        elseif strcmp(compareSides,'stimmove')
            if a == 1
                [meanResp, semResp, plotColors] = computeContrastCurves(alignedResps{a}(:,:,plotCells(k)),expInfo,behavioralData,initTrialConditions('movementDir','cw','movementTime','late'),'balanced');
                plotContrastCurves(eventWindow, meanResp, semResp, plotColors,'-');
                yMin(end+1) = min(min(min(meanResp-semResp)));
                yMax(end+1) = max([.1 max(max(max(meanResp+semResp)))]);
                [meanResp, semResp, plotColors] = computeContrastCurves(alignedResps{a}(:,:,plotCells(k)),expInfo,behavioralData,initTrialConditions('movementDir','ccw','movementTime','late'),'balanced');
                plotContrastCurves(eventWindow, meanResp, semResp, plotColors,':');
                xlim([-1 2]);
                yMin(end+1) = min(min(min(meanResp-semResp)));
                yMax(end+1) = max([.1 max(max(max(meanResp+semResp)))]);
            elseif a == 2
                [meanResp, semResp, plotColors] = computeContrastCurves(alignedResps{a}(:,:,plotCells(k)),expInfo,behavioralData,initTrialConditions('movementDir','cw','movementTime','late'),'balanced');
                plotContrastCurves(eventWindow, meanResp, semResp, plotColors,'-');
                yMin(end+1) = min(min(min(meanResp-semResp)));
                yMax(end+1) = max([.1 max(max(max(meanResp+semResp)))]);
                [meanResp, semResp, plotColors] = computeContrastCurves(alignedResps{a}(:,:,plotCells(k)),expInfo,behavioralData,initTrialConditions('movementDir','ccw','movementTime','late'),'balanced');
                plotContrastCurves(eventWindow, meanResp, semResp, plotColors,':');
                xlim([-2 1]);
                yMin(end+1) = min(min(min(meanResp-semResp)));
                yMax(end+1) = max([.1 max(max(max(meanResp+semResp)))]);
            elseif a == 3
                [meanResp, semResp, plotColors] = computeRewardCurves(alignedResps{a}(:,:,plotCells(k)),expInfo,behavioralData,initTrialConditions('movementDir','cw','movementTime','late'));
                plotContrastCurves(eventWindow, meanResp, semResp, plotColors,'-');
                yMin(end+1) = min(min(min(meanResp-semResp)));
                yMax(end+1) = max([.1 max(max(max(meanResp+semResp)))]);
                [meanResp, semResp, plotColors] = computeRewardCurves(alignedResps{a}(:,:,plotCells(k)),expInfo,behavioralData,initTrialConditions('movementDir','ccw','movementTime','late'));
                plotContrastCurves(eventWindow, meanResp, semResp, plotColors,':');
                xlim([-2 1]);
                yMin(end+1) = min(min(min(meanResp-semResp)));
                yMax(end+1) = max([.1 max(max(max(meanResp+semResp)))]);
            end
        end
            
        ln = line([0 0],[0 1],'LineStyle','-','Color',[0 0 0],'linewidth',1);
        uistack(ln,'bottom');
        set(gca,'xtick',[]);
        set(gca,'TickDir','out');
        
        yMin(end+1) = min(min(min(meanResp-semResp)));
        yMax(end+1) = max([.1 max(max(max(meanResp+semResp)))]);
    end
    
    %unify axis limits across plots and add labels as appropriate
    for a = 1:length(allConds)
        for iCond = 1:length(allConds{a})
            [spidx1, spidx2] = goToRasterSubplot(length(allConds), total_length, tl{a}, a, iCond, 1, psth);
            subplot(total_length,length(allConds),[spidx1 spidx2])
            caxis([min(iMin) max(.4)]);
            if iCond == length(allConds{a})
                xlabel('Time (s)')
            end
        end
        spidxA = sub2ind([length(allConds) total_length], a, 1);
        spidxB = sub2ind([length(allConds) total_length], a, psth);    
        subplot(total_length,length(allConds),[spidxA spidxB])
        box off
        ylim([min(0) max(yMax)]);
        if a == 1 
            ylabel('Activity')
        end

    end    

    was_a_key = waitforbuttonpress;
    
    if was_a_key && strcmp(get(fig, 'CurrentKey'), 'leftarrow')
      k = max(1, k - 1);
    elseif was_a_key && strcmp(get(fig, 'CurrentKey'), 'rightarrow')
      k = min(max_k, k + 1);
    elseif was_a_key && strcmp(get(fig, 'CurrentKey'), 'return')
        disp(k)
        if strcmp(compareSides,'blocks')
            figName = strcat(expInfo.mouseName,'_',expInfo.expDate,'_cell_',num2str(plotCells(k)),'_blocks');
        elseif strcmp(compareSides,'stimmove')
            figName = strcat(expInfo.mouseName,'_',expInfo.expDate,'_cell_',num2str(plotCells(k)),'_stimmove');
        else
            figName = strcat(expInfo.mouseName,'_',expInfo.expDate,'_cell_',num2str(plotCells(k)));
        end
        printfig(gcf, figName)
        break
    end
    
end
    