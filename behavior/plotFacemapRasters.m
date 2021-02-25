function figName = plotFacemapRasters(expInfo, behavioralData,eyeData,whichROIs, whichTrials, k)


%% initialize experiment details
alignedFace = eyeData.eta.alignedFace;
eventWindow = eyeData.eta.eventWindow;
et = behavioralData.eventTimes;
wm = behavioralData.wheelMoves;


%% 
%define ROIs from facemap
plotROIs = whichROIs; 

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
rasterLabels = {'Whisk' 'Non-whisk trials'};
psthColors = [1 0 1; 0 1 1];


%% plot all facemap ROIs (all trials)separated by whisking trials vs non-whisk trials
rasterColors = colormap(winter); 
%rasterColors = [0 0 1; 1 0 0];
rasterLabels = {'Whisk' 'Non-whisk trials'};
%use blue and red for whisk vs non-whisk conditions 
%psthColors = [0 0 1; 1 0 0];
%use a colourmap for 4 diff lines garded by whisking 
psthColors = colormap(winter);

facemapfig = figure;
set(facemapfig, 'Position', [80 250 870 660]);
hold on;

if ~exist('k') == 1
    k = 1;
end

max_k = length(plotROIs);

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
            %[relativeTimes, sortIdx] = sortTrialTimes(whichTrials, et, wm);
            
            % sort the trials by whisking
            [sortIdx] = sortTrialByWhisk(whichTrials,eyeData);
           
            [spidx1, spidx2] = goToRasterSubplot(length(trialLists), total_length, cellfun(@length, trialLists{a})', a, iCond, 1, psth);
            ax = subplot(total_length,length(trialLists),[spidx1 spidx2]);
            f = imagesc(eventWindow,1:numTrials,alignedFace{a}(whichTrials(sortIdx),:,plotROIs(k)));
            colormap(ax,rasterColor(rasterColors(iCond,:)));
            hold on;
            %             line([0 0],[1 length(whichTrials)],'LineStyle','-','Color',[0 0 0]);
            %mt = plot(relativeTimes(:,2,a),1:numTrials,'bo');
            %st = plot(relativeTimes(:,1,a),1:numTrials,'ko');
            %rt = plot(relativeTimes(:,3,a),1:numTrials,'ko');
            %set([st mt rt], 'MarkerSize',1,'MarkerFaceColor','k','MarkerEdgeColor','none');
            box off
            ln = line([0 0],[-1 10],'LineStyle','-','Color',[0 0 0],'linewidth',1);
            uistack(ln,'bottom');
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
    
%{
        %plot 1 mean psth at the top
        spidxA = sub2ind([length(trialLists) total_length], a, 1);
        spidxB = sub2ind([length(trialLists) total_length], a, psth);
        subplot(total_length,length(trialLists),[spidxA spidxB])
        [meanPSTH, semPSTH, rasters] = computePSTHs(alignedFace{a}(:,:,plotROIs(k)),trialLists{a});
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
        %}
        
        %plot psth with 4 lines of diff whisk level
        %first, divide trials into 4 quesrtiles 
        
        q=0.25*(length(trialLists{1,1}));
        quartile{1}=trialLists{1,1}(1:q);
        quartile{2}=trialLists{1,1}(q+1:2*q);
        quartile{3}=trialLists{1,1}(2*q+1:3*q-1);
        quartile{4}=trialLists{1,1}(end-q:end);

        for f =1:length(quartile)
            spidxA = sub2ind([length(trialLists) total_length], a, 1);
            spidxB = sub2ind([length(trialLists) total_length], a, psth);
            subplot(total_length,length(trialLists),[spidxA spidxB])
            [meanPSTH, semPSTH, rasters] = computePSTHs(alignedFace{a}(:,:,plotROIs(k)),trialLists{a});
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
            caxis([min(iMin) max(iMax)*.1]);
            if iCond == length(trialLists{a})
                if a == 1
                    xlabel('Time from stimulus (s)')
                elseif a ==2
                    xlabel('Time from movement (s)')
                   
                elseif a == 3
                    xlabel('Time from outcome (s)')
                end
            %title(strcat('ROI',num2str(plotROIs(k))))
            end
        end
        spidxA = sub2ind([length(trialLists) total_length], a, 1);
        spidxB = sub2ind([length(trialLists) total_length], a, psth);    
        subplot(total_length,length(trialLists),[spidxA spidxB])
        box off
        ylim([min(yMin) max(yMax)]);
        if a == 1 
            ylabel('Arousal')
        end

    end    

    
    was_a_key = waitforbuttonpress;
    
    if was_a_key && strcmp(get(facemapfig, 'CurrentKey'), 'leftarrow')
      k = max(1, k - 1);
    elseif was_a_key && strcmp(get(facemapfig, 'CurrentKey'), 'rightarrow')
      k = min(max_k, k + 1);
    elseif was_a_key && strcmp(get(facemapfig, 'CurrentKey'), 'return')
        disp(k)
         saveName = strcat('C:\Users\Ella Svahn\Documents\eyedata\LEW031\Rasters\',expInfo.mouseName,'_',expInfo.expDate,'_ROI_',num2str(plotROIs(k)));
        print(gcf,'-dpng',saveName)
        break
    end
end
%%
%see if whisk ROI (2) correlate with ROI 5
[r, p] = corrcoef(eyeData.eta.alignedFace{1}(:,:,2),eyeData.eta.alignedFace{1}(:,:,5));