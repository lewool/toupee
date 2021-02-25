function figName = plotWhiskRasters(expInfo, behavioralData, neuralData, whichCells, whichTrials,whichSort k)
%whichSort can be either:
%'byEventTime' - this sorts trials by 1st move time
%or 'byWhisk' - this sorts trials by mean pre-stimulus whsiking (120-0ms before stim) 
%error needs fixing for sorting, change  on line 85

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
    plotCells = whichCells;
end

%%
for iA = 1:3
    for cond = 1:length(whichTrials)
        trialLists{iA}{cond,1} = whichTrials{cond};
    end
end

%% set up some plot values

%compute the total size of the figure based on trials
totalTrials = sum(cellfun(@length, trialLists{1}));
totalRasters = max(cellfun(@numel, trialLists));
border = 1;
total_raster = totalTrials + border*(totalRasters+1);
psth = round(total_raster*.5);
total_length = total_raster + psth;

rasterColors = [1 0 1; 0 1 1];
rasterLabels = {'Whisk' 'Quiesc'};
psthColors = [1 0 1; 0 1 1];

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
            
            if whichSort == 'byWhisk'
                %sort the trials by whisking
                [relativeTimes,sortIdx] = sortTrialByWhisk(whichTrials,eyeData,et,wm);
            elseif whichSort == 'byEventTimes'
                % sort the trials by event time, here 1st move
                [relativeTimes, sortIdx] = sortTrialTimes(whichTrials, et, wm);           
            else
                disp("Input sorting type: EITHER 'byEventTime' OR 'byWhisk'")
            end
            
            [spidx1, spidx2] = goToRasterSubplot(length(trialLists), total_length, cellfun(@length, trialLists{a})', a, iCond, 1, psth);
            ax = subplot(total_length,length(trialLists),[spidx1 spidx2]);
            f = imagesc(eventWindow,1:numTrials,alignedResps{a}(whichTrials(sortIdx),:,plotCells(k)));
            colormap(ax,rasterColor(rasterColors(iCond,:)));
            hold on;
            %             line([0 0],[1 length(whichTrials)],'LineStyle','-','Color',[0 0 0]);
            mt = plot(relativeTimes(:,2,a),1:numTrials,'bo');
            st = plot(relativeTimes(:,1,a),1:numTrials,'ko');
            rt = plot(relativeTimes(:,3,a),1:numTrials,'ko');
            set([st mt rt], 'MarkerSize',1,'MarkerFaceColor','k','MarkerEdgeColor','none');
            box off
            set(gca,'ytick',[]);
            set(gca,'tickdir','out')
            ylabel(rasterLabels(iCond))
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
        [meanPSTH, semPSTH, rasters] = computePSTHs(alignedResps{a}(:,:,plotCells(k)),trialLists{a});
        for d = 1:size(meanPSTH,2)
            if d == 1
                ls = '-';
            else
                ls = ':';
            end
            plotPSTHs(eventWindow, cell2mat(meanPSTH(:,d)), cell2mat(semPSTH(:,d)), psthColors,ls);
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
                if a == 1
                    xlabel('Time from stimulus (s)')
                elseif a ==2
                    xlabel('Time from movement (s)')
                elseif a == 3
                    xlabel('Time from outcome (s)')
                end
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
        saveName = strcat('C:\Users\Ella Svahn\Documents\eyedata\LEW031\Rasters\',expInfo.mouseName,'_',expInfo.expDate,'_cell_',num2str(plotCells(k)));
        print(gcf,'-dpng',saveName)
        break
    end
end