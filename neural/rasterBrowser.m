function figName = rasterBrowser(expInfo, behavioralData, neuralData, whichCells,pickTrials, trialStruct, k)

% pickTrials = {'side_direction', 'side_direction', 'outcome_direction'};
% contrastOverride = 'contrast_direction';
% trialStruct = trialTypes.intVar.cb2D;
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

%%


for iA = 1:length(pickTrials)
    trialLists{iA} = reshape(getfield(trialStruct,pickTrials{iA})',[1 numel(getfield(trialStruct,pickTrials{iA}))])';
    trialArrays{iA} = getfield(trialStruct,pickTrials{iA});
end

%% set up some plot values

%compute the total size of the figure based on trials
totalTrials = sum(sum(cellfun(@length, getfield(trialStruct,pickTrials{1}))));
totalRasters = max(cellfun(@numel, trialLists));
border = 1;
total_raster = totalTrials + border*(totalRasters+1);
psth = round(total_raster*.5);
total_length = total_raster + psth;

for iETA = 1:3
    if contains(pickTrials{iETA},'side')
        if contains(pickTrials{iETA},'_direction') || contains(pickTrials{iETA},'_block')
            rasterColors{iETA} = [0 0 1; 0 0 1; 0 0 0; 0 0 0; 1 0 0; 1 0 0];
            rasterLabels{iETA} = {'L' 'L' '0' '0' 'R' 'R'};
            psthColors{iETA} = [0 .4 1; .75 .75 .75; 1 0 0];
        else
            rasterColors{iETA} = [0 0 1; 0 0 0; 1 0 0];
            rasterLabels{iETA} = {'Left' 'Zero' 'Right'};
            psthColors{iETA} = [0 .4 1; .75 .75 .75; 1 0 0];
        end
    elseif contains(pickTrials{iETA},'outcome')
        if contains(pickTrials{iETA},'_direction') || contains(pickTrials{iETA},'_block')
            rasterColors{iETA} = [.5 1 .5; .5 1 .5; 1 0 0; 1 0 0];
            rasterLabels{iETA} = {'Cor' 'Cor' 'Inc' 'Inc'};
            psthColors{iETA} = [.1 .7 .1; .75 0 0];
        else
            rasterColors{iETA} = [.5 1 .5; 1 0 0];
            rasterLabels{iETA} = {'Cor' 'Inc'};
            psthColors{iETA} = [.1 .7 .1; .75 0 0];
        end
    elseif contains(pickTrials{iETA},'direction')
        if contains(pickTrials{iETA},'_block')
            rasterColors{iETA} = [0 0 1; 0 0 1; 1 0 0; 1 0 0];
            rasterLabels{iETA} = {'L' 'L' 'R' 'R'};
            psthColors{iETA} = [0 .4 1; 1 0 0];
        else
            rasterColors{iETA} = [0 0 1; 1 0 0];
            rasterLabels{iETA} = {'L' 'R'};
            psthColors{iETA} = [0 .4 1; 1 0 0];
        end
    elseif contains(pickTrials{iETA},'block')
        rasterColors{iETA} = [0 0 1; 1 0 0];
        rasterLabels{iETA} = {'L' 'R'};
        psthColors{iETA} = [0 0 1; .5 0 0];
    end
end

%% plot (all trials)
fig = figure;
set(fig, 'Position', [80 250 870 660]);
hold on;

if ~exist('k') == 1
    k = 1;
end

max_k = length(plotCells);

while k <= max_k

    %clear subplots
    for a = 1:length(trialLists)
        for iCond = 1:length(trialLists{a})
            [spidx1, spidx2] = goToRasterSubplot(length(trialLists), total_length, cellfun(@length, trialLists{a})', a, iCond, 1, psth);
            subplot(total_length,length(trialLists),[spidx1 spidx2])
            cla;
        end
        spidxA = sub2ind([length(trialLists) total_length], a, 1);
        spidxB = sub2ind([length(trialLists) total_length], a, psth);    
        subplot(total_length,length(trialLists),[spidxA spidxB])
        cla;
    end
    
    iMin = [];
    iMax = [];
    yMin = [];
    yMax = [];

    for a = 1:3
        for iCond = 1:length(trialLists{a})
            
            % extract the relevant trials for that condition
            whichTrials = trialLists{a}{iCond};
            numTrials = size(whichTrials,2);
            
            % sort the trials by event time
            [relativeTimes, sortIdx] = sortTrialTimes(whichTrials, et, wm);
            
            [spidx1, spidx2] = goToRasterSubplot(length(trialLists), total_length, cellfun(@length, trialLists{a})', a, iCond, 1, psth);
            ax = subplot(total_length,length(trialLists),[spidx1 spidx2]);
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
            ylabel(rasterLabels{a}(iCond))
            if iCond < length(trialLists{a})
                set(gca,'xtick',[]);
            end
            iMax(end+1) = max(max(f.CData));
            iMin(end+1) = min(min(f.CData));
            if a > 1
                xlim([-1.5 1]);
            else
                xlim([-.5 1.5]);
            end
            set(gca,'xtick',[-1 0 1])
        end
    

        %plot psth at the top
        spidxA = sub2ind([length(trialLists) total_length], a, 1);
        spidxB = sub2ind([length(trialLists) total_length], a, psth);
        subplot(total_length,length(trialLists),[spidxA spidxB])
        [meanPSTH, semPSTH, rasters] = computePSTHs(alignedResps{a}(:,:,plotCells(k)),trialArrays{a});
        for d = 1:size(meanPSTH,2)
            if d == 1
                ls = '-';
            else
                ls = ':';
            end
            plotPSTHs(eventWindow, cell2mat(meanPSTH(:,d)), cell2mat(semPSTH(:,d)), psthColors{a},ls);
            if a > 1
                xlim([-1.5 1]);
            else
                xlim([-.5 1.5]);
            end
            yMin(end+1) = min(min(cell2mat(meanPSTH)-cell2mat(semPSTH)));
            yMax(end+1) = max([.01 max(max(cell2mat(meanPSTH)+cell2mat(semPSTH)))]);
        end
        
        ln = line([0 0],[-1 10],'LineStyle','-','Color',[0 0 0],'linewidth',1);
        uistack(ln,'bottom');
        set(gca,'xtick',[]);
        set(gca,'TickDir','out');
    end
    
    %unify axis limits across plots and add labels as appropriate
    for a = 1:3
        for iCond = 1:length(trialLists{a})
            [spidx1, spidx2] = goToRasterSubplot(length(trialLists), total_length, cellfun(@length, trialLists{a})', a, iCond, 1, psth);
            subplot(total_length,length(trialLists),[spidx1 spidx2])
            caxis([min(iMin) max(iMax)*.5]);
            if iCond == length(trialLists{a})
                xlabel('Time (s)')
            end
        end
        spidxA = sub2ind([length(trialLists) total_length], a, 1);
        spidxB = sub2ind([length(trialLists) total_length], a, psth);    
        subplot(total_length,length(trialLists),[spidxA spidxB])
        box off
        ylim([min(yMin) max(yMax)]);
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
        if contains(pickTrials,'block')
            figName = strcat(expInfo.mouseName,'_',expInfo.expDate,'_cell_',num2str(plotCells(k)),'_blocks');
        elseif contains(pickTrials,'direction')
            figName = strcat(expInfo.mouseName,'_',expInfo.expDate,'_cell_',num2str(plotCells(k)),'_stimmove');
        else
            figName = strcat(expInfo.mouseName,'_',expInfo.expDate,'_cell_',num2str(plotCells(k)));
        end
        printfig(gcf, figName)
        break
    end
end
