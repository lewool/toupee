function figName = plotFacemapRasters(expInfo, behavioralData,eyeData,whichROIs, whichTrials, k)


%% initialize experiment details
alignedFace = eyeData.eta.alignedFace;
eventWindow = eyeData.eta.eventWindow;
et = behavioralData.eventTimes;
wm = behavioralData.wheelMoves;


%% 
%define ROIs from facemap
plotROIs = whichROIs; 

%% set up some plot values

%compute the total size of the figure based on trials
totalTrials = sum(cellfun(@length, whichTrials));
totalRasters = max(cellfun(@numel, whichTrials));
border = 1;
total_raster = totalTrials + border*(totalRasters+1);
psth = round(total_raster*.5);
total_length = total_raster + psth;

rasterColors = [1 0 1; 0 1 1];
rasterLabels = {'Whisk' 'Non-whisk trials'};
psthColors = [1 0 1; 0 1 1];

%% plot all facemap ROIs (all trials)separated by whisking trials vs non-whisk trials
rasterColors = [1 0 0; 0 0 1];
rasterLabels = {'Whisk' 'Non-whisk trials'};
psthColors = [1 0 0; 0 0 1];

facemapfig = figure;
set(facemapfig, 'Position', [80 250 870 660]);
hold on;

if ~exist('k') == 1
    k = 1;
end

max_k = length(plotROIs);

while k <= max_k

    %clear subplots
    for a = 1:length(whichTrials)
        for iCond = 1:length(whichTrials)
            [spidx1, spidx2] = goToRasterSubplot(length(whichTrials), total_length, cellfun(@length, whichTrials)', a, iCond, 1, psth);
            subplot(total_length,length(whichTrials),[spidx1 spidx2])
            cla;
        end
        
        spidxA = sub2ind([(length(whichTrials)) total_length], a, 1);
        spidxB = sub2ind([(length(whichTrials)) total_length], a, psth);    
        subplot(total_length,length(whichTrials),[spidxA spidxB])
        cla;
    end
    
    iMin = [];
    iMax = [];
    yMin = [];
    yMax = [];

    for a = 1:3
        for iCond = 1:length(whichTrials)
            % extract the relevant trials for that condition
            %= whichTrials{iCond};
            numTrials = length(whichTrials);
            
            % sort the trials by event time
            %[relativeTimes, sortIdx] = sortTrialTimes(whichTrials, et, wm);
            
            % sort the trials by whisking
            [relativeTimes, sortIdx] = sortTrialByWhisk(whichTrials,eyeData,et, wm);
            
            %plot 2 lines in PSTH for the 25% of trials wiht the most and
            %least whisking
            %firstQuarter = whichTrials(1:(length(whichTrials)/4));
            %lastQuarter = whichTrials(end-(length(whichTrials)/4:end));

        
            %[spidx1, spidx2] = goToRasterSubplot((length(whichTrials), total_length, cellfun(@length, whichTrials)', a, iCond, 1, psth);
            ax = subplot(total_length,length(whichTrials),[spidx1 spidx2]);
            f = imagesc(eventWindow,1:numTrials,alignedFace{a}(whichTrials(sortIdx),:,plotROIs(k)));
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
            if iCond < length(whichTrials)
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
        spidxA = sub2ind([length(whichTrials) total_length], a, 1);
        spidxB = sub2ind([length(whichTrials) total_length], a, psth);
        subplot(total_length,length(whichTrials),[spidxA spidxB])
        [meanPSTH, semPSTH, rasters] = computePSTHs(alignedFace{a}(:,:,plotROIs(k)),whichTrials);
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
        for iCond = 1:length(whichTrials)
            [spidx1, spidx2] = goToRasterSubplot(length(whichTrials), total_length, cellfun(@length, whichTrials{a})', a, iCond, 1, psth);
            subplot(total_length,length(whichTrials),[spidx1 spidx2])
            caxis([min(iMin) max(iMax)*.3]);
            if iCond == length(whichTrials)
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
        spidxA = sub2ind([length(whichTrials) total_length], a, 1);
        spidxB = sub2ind([length(whichTrials) total_length], a, psth);    
        subplot(total_length,length(whichTrials),[spidxA spidxB])
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
disp(strcat('Correlation ROI2 & 5: r=',r))






