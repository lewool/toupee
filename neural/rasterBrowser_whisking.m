function figName = rasterBrowser_whisking(expInfo, behavioralData, neuralData, whichCells, whichTrials, eyeData, k)


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
if strcmp(whichTrials, 'all')
    whichTrials = 1:size(alignedResps{1},1);
end
for iA = 1:3
    trialLists{iA}{1,1} = whichTrials;
end

%% set up some plot values

%compute the total size of the figure based on trials
totalTrials = sum(cellfun(@length, trialLists{1}));
totalRasters = max(cellfun(@numel, trialLists));
border = 1;
total_raster = totalTrials + border*(totalRasters+1);
psth = round(total_raster*.5);
total_length = total_raster + psth;

rasterColors = [0 0 0];
<<<<<<< HEAD
rasterLabels = {'Trials ranked by pre-stim whisking'};
=======
rasterLabels = {'Trials ranked by pretrial whisking'};
>>>>>>> master
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
            
            % sort the trials by pre-trial whisking
            [relativeTimes, sortIdx] = sortTrialsByWhisking(whichTrials, eyeData, et, wm);
            
            [spidx1, spidx2] = goToRasterSubplot(length(trialLists), total_length, cellfun(@length, trialLists{a})', a, iCond, 1, psth);
            ax = subplot(total_length,length(trialLists),[spidx1 spidx2]);
            f = imagesc(eventWindow,1:numTrials,alignedResps{a}(whichTrials(sortIdx),:,plotCells(k)));
            colormap(ax,rasterColor(rasterColors(iCond,:)));
            hold on;
            %             line([0 0],[1 length(whichTrials)],'LineStyle','-','Color',[0 0 0]);
            mt = plot(relativeTimes(:,2,a),1:numTrials,'bo');
            st = plot(relativeTimes(:,1,a),1:numTrials,'ko');
            rt = plot(relativeTimes(:,3,a),1:numTrials,'ko');
<<<<<<< HEAD
            set([st], 'MarkerSize',1,'MarkerFaceColor','k','MarkerEdgeColor','none');
            set([mt], 'MarkerSize',1,'MarkerFaceColor','r','MarkerEdgeColor','none');
            set([rt], 'MarkerSize',1,'MarkerFaceColor','b','MarkerEdgeColor','none');
=======
            set([st mt rt], 'MarkerSize',1,'MarkerFaceColor','k','MarkerEdgeColor','none');
>>>>>>> master
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
        quartileLength = floor(length(whichTrials)/4);
        quartileColors = [...
            1 0 0;...
            1 .25 0;...
            1 .5 0;...
            1 .75 0];
        quartiles = [...
            whichTrials(sortIdx(1:quartileLength));...
            whichTrials(sortIdx(quartileLength+1:2*quartileLength));...
            whichTrials(sortIdx(2*quartileLength+1:3*quartileLength));...
            whichTrials(sortIdx(3*quartileLength+1:4*quartileLength))];
        
        for iQ = 1:4
            qIdx = quartiles(iQ,:);
        [meanPSTH, semPSTH, rasters] = computePSTHs(alignedResps{a}(:,:,plotCells(k)),qIdx);
        for d = 1:size(meanPSTH,2)
            if d == 1
                ls = '-';
            else
                ls = ':';
            end
            plotPSTHs(eventWindow, cell2mat(meanPSTH(:,d)), cell2mat(semPSTH(:,d)), quartileColors(iQ,:),ls);
            if a > 1
                xlim([-1.5 1]);
            else
                xlim([-.5 1.5]);
            end
            yMin(end+1) = min(min(cell2mat(meanPSTH)-cell2mat(semPSTH)));
            yMax(end+1) = max([.01 max(max(cell2mat(meanPSTH)+cell2mat(semPSTH)))]);
        end
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
<<<<<<< HEAD
            title(strcat(expInfo.mouseName,' / ',expInfo.expDate,' /c ',num2str(plotCells(k))),'FontSize',12);
=======
>>>>>>> master
        end

    end    

    
    was_a_key = waitforbuttonpress;
    
    if was_a_key && strcmp(get(fig, 'CurrentKey'), 'leftarrow')
      k = max(1, k - 1);
    elseif was_a_key && strcmp(get(fig, 'CurrentKey'), 'rightarrow')
      k = min(max_k, k + 1);
    elseif was_a_key && strcmp(get(fig, 'CurrentKey'), 'return')
        disp(strcat({'k = '},num2str(k)))
<<<<<<< HEAD
        saveName = strcat('C:\Users\Ella Svahn\Documents\eyedata\LEW031\Rasters\','whiskQuart',expInfo.mouseName,'_',expInfo.expDate,'_cell_',num2str(plotCells(k)));
        print(gcf,'-dpng',saveName)   
        break
    elseif was_a_key && strcmp(get(fig, 'CurrentKey'), 'escape')
        disp(strcat({'k = '},num2str(k)))
        saveName = strcat('C:\Users\Ella Svahn\Documents\eyedata\LEW031\Rasters\','whiskQuart',expInfo.mouseName,'_',expInfo.expDate,'_cell_',num2str(plotCells(k)));
        print(gcf,'-dpng',saveName)        
=======
        figName = strcat('whiskQuartiles_',expInfo.mouseName,'_',expInfo.expDate,'_cell_',num2str(plotCells(k)));
        printfig(gcf, figName)
        break
    elseif was_a_key && strcmp(get(fig, 'CurrentKey'), 'escape')
        disp(strcat({'k = '},num2str(k)))
        figName = strcat('whiskQuartiles_',expInfo.mouseName,'_',expInfo.expDate,'_cell_',num2str(plotCells(k)));
>>>>>>> master
        close(fig)
        break
    end
end
