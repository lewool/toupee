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
nt = numel(expInfo.block.events.endTrialValues);
yOverride = 1;

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

ETAs = [1 4 2 3];

for iETA = 1:length(ETAs)
    if contains(pickTrials{iETA},'side')
        if contains(pickTrials{iETA},'_direction') || contains(pickTrials{iETA},'_block') || contains(pickTrials{iETA},'_outcome')
            rasterColors{iETA} = [0 .4 1; 0 .4 1; 0 0 0; 0 0 0; 1 0 0; 1 0 0];
            rasterLabels{iETA} = {'L' 'L' '0' '0' 'R' 'R'};
            psthColors{iETA} = [0 .4 1; .75 .75 .75; 1 0 0];
        else
            rasterColors{iETA} = [0 0 1; 0 0 0; 1 0 0];
            rasterLabels{iETA} = {'Left' 'Zero' 'Right'};
            psthColors{iETA} = [0 .4 1; .75 .75 .75; 1 0 0];
        end
    elseif contains(pickTrials{iETA},'outcome')
        if contains(pickTrials{iETA},'_direction') || contains(pickTrials{iETA},'_block')
            rasterColors{iETA} = [.1 .7 .1; .1 .7 .1; .75 0 0; .75 0 0];
            rasterLabels{iETA} = {'Cor' 'Cor' 'Inc' 'Inc'};
            psthColors{iETA} = [.1 .7 .1; .75 0 0];
        else
            rasterColors{iETA} = [.5 1 .5; 1 0 0];
            rasterLabels{iETA} = {'Cor' 'Inc'};
            psthColors{iETA} = [.1 .7 .1; .75 0 0];
        end
    elseif contains(pickTrials{iETA},'reward')
        if contains(pickTrials{iETA},'_direction') 
            psthColors{iETA} = [0 .5 0; .5 1 .5; .7 0.1 0];
            rasterLabels{iETA} = {'Hi L' 'Hi R' 'Lo L' 'Lo R' 'No L' 'No R'};
            rasterColors{iETA} = [0 .5 0; 0 .5 0; .5 1 .5; .5 1 .5; .7 0.1 0; .7 0.1 0];
        else
            rasterColors{iETA} = [0 .5 0; .5 1 .5; .7 .1 0];
            rasterLabels{iETA} = {'Hi' 'Lo' 'No'};
            psthColors{iETA} = [0 .5 0; .5 1 .5; .7 .1 0];
        end
    elseif contains(pickTrials{iETA},'direction') && ~contains(pickTrials{iETA},'_direction')
        if contains(pickTrials{iETA},'_block')
%             rasterColors{iETA} = [.5 0 1; .5 0 1; 1 0 .5; 1 0 .5];
            rasterColors{iETA} = [0 0 0; 0 0 0; 0 0 0; 0 0 0];
            rasterLabels{iETA} = {'L' 'L' 'R' 'R'};
%             psthColors{iETA} = [.5 0 1; 1 0 .5];
            psthColors{iETA} = [0 0 0; 0 0 0];
        else
            rasterColors{iETA} = [0 0 0; 0 0 0];
            rasterLabels{iETA} = {'L' 'R'};
            psthColors{iETA} = [0 0 0;0 0 0]; %[0.3 0 1; 1 0.2 0.8]
        end
    elseif contains(pickTrials{iETA},'block') && ~contains(pickTrials{iETA},'_block')
        if contains(pickTrials{iETA},'_direction') 
            psthColors{iETA} = [0.1 .7 .1; 1 .6 0];
            rasterLabels{iETA} = {'LL' 'LR' 'RL' 'RR'};
            rasterColors{iETA} = [0.1 .7 .1; 0.1 .7 .1; 1 .6 0; 1 .6 0];
        else
            rasterColors{iETA} = [0.1 .7 .1; 1 .6 0];
            rasterLabels{iETA} = {'L' 'R'};
            psthColors{iETA} = [0.1 .7 .1; 1 .6 0];
        end
             
    elseif contains(pickTrials{iETA},'contrast')
        contrasts = getUniqueContrasts(expInfo);
        zeroIdx = find(contrasts == 0);
        walkup = length(contrasts) - zeroIdx;
        walkback = zeroIdx - 1;
        allColors = [0 0 .5; 0 0 1; 0 .4 1; .6 .8 1; .75 .75 .75; .8 .45 .45; 1 0 0; .5 0 0; .25 0 0];
        zeroGray = find(allColors(:,1) == .75);
        rasterColors{iETA} = allColors(zeroGray-walkback:zeroGray + walkup,:);
        rasterLabels{iETA} = {'100' '50' '12' '5' '0' '5' '12' '50' '100'};
        psthColors{iETA} = allColors(zeroGray-walkback:zeroGray + walkup,:);
    end
end

%% behavioral stuff to plot
try
    allTrials = cat(2, trialStruct.contrast{:});
catch
    allTrials = cat(2, trialStruct.contrast_direction{:});
end

TTG = et(2).daqTime(allTrials) - et(1).daqTime(allTrials);
medianTTG = median(TTG);
lbTTG = prctile(TTG,2.5);
ubTTG = prctile(TTG,97.5);
bins = 0:.02:1;
densityTTG = histcounts(TTG,bins)/max(histcounts(TTG,bins));

TTM = wm.epochs(5).onsetTimes(allTrials) - et(2).daqTime(allTrials);
medianTTM = median(TTM);
lbTTM = prctile(TTM,2.5);
ubTTM = prctile(TTM,97.5);
bins = 0:.02:1;
densityTTM = histcounts(TTM,bins)/max(histcounts(TTM,bins));

TTF = et(5).daqTime(allTrials) - wm.epochs(5).onsetTimes(allTrials);
medianTTF = median(TTF);
lbTTF = prctile(TTF,2.5);
ubTTF = prctile(TTF,97.5);
bins = 0:.02:1;
densityTTF = histcounts(TTF,bins)/max(histcounts(TTF,bins));

allTrials(allTrials == nt) = [];
TTS = et(1).daqTime(allTrials+1) - wm.epochs(5).onsetTimes(allTrials);
medianTTS = median(TTS);
lbTTS = prctile(TTS,2.5);
ubTTS = prctile(TTS,97.5);
bins = 0:.02:1;
densityTTS = histcounts(TTS,bins)/max(histcounts(TTS,bins));

% figure;
% hold on;
% for t = 1:length(densityTTF)
%     fill([bins(t) bins(t+1) bins(t+1) bins(t)],[0 0 1 1],'k','EdgeColor','none','FaceAlpha',densityTTF(t));
% end

%% plot (all trials)
fig = figure;
set(fig, 'Position', [80 250 870 660]);
set(fig, 'Position', [460 800 524 664]);
set(fig, 'Position', [460 860 950 605]);
hold on;

if ~exist('k') == 1
    k = 1;
end

max_k = length(plotCells);

while k <= max_k

    %clear subplots
    for a = 1:length(trialLists)
        for iCond = 1:length(trialLists{a})
            [spidx1, spidx2] = goToRasterSubplot(length(trialLists), total_length, cellfun(@length, trialLists{a})', a, iCond, border, psth);
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

    for a = 1:length(ETAs)
        for iCond = 1:length(trialLists{a})
            
            % extract the relevant trials for that condition
            whichTrials = trialLists{a}{iCond};
            numTrials = size(whichTrials,2);
            
            % sort the trials by event time
            [relativeTimes, sortIdx] = sortTrialTimes(whichTrials, et, wm);
            
            [spidx1, spidx2] = goToRasterSubplot(length(trialLists), total_length, cellfun(@length, trialLists{a})', a, iCond, 1, psth);
            ax = subplot(total_length,length(trialLists),[spidx1 spidx2]);
            f = imagesc(eventWindow,1:numTrials,alignedResps{ETAs(a)}(whichTrials(sortIdx),:,plotCells(k)));
            colormap(ax,rasterColor(rasterColors{a}(iCond,:)));
            hold on;
            line([0 0],[1 length(whichTrials)],'LineStyle',':','Color',[0 0 0]);
            if a == 1
                lv = medianTTG;
                 l1 = line([lv lv],[1 length(whichTrials)],'LineStyle',':','Color',[0 0 0],'linewidth',.5);
            elseif a == 2
                lv = medianTTM;
                 l1 = line([lv lv],[1 length(whichTrials)],'LineStyle',':','Color',[0 0 0],'linewidth',.5);
            elseif a == 3
                lv = medianTTF;
                 l1 = line([lv lv],[1 length(whichTrials)],'LineStyle',':','Color',[0 0 0],'linewidth',.5);
            end
           
%             mt = plot(relativeTimes(:,2,a),1:numTrials,'bo');
%             st = plot(relativeTimes(:,1,a),1:numTrials,'ko');
%             rt = plot(relativeTimes(:,3,a),1:numTrials,'ko');
%             set([st mt rt], 'MarkerSize',1,'MarkerFaceColor','k','MarkerEdgeColor','none');
            box off
            set(gca,'ytick',[]);
            set(gca,'tickdir','out')
            ylabel(rasterLabels{a}(iCond))
            
            iMax(end+1) = max(max(f.CData));
            iMin(end+1) = min(min(f.CData));
           if a > 1
                xlim([-.5 2]);
            else
                xlim([-.5 2]);
            end
            set(gca,'xtick',[0 1])
            if iCond < length(trialLists{a})
                set(gca,'xtick',[]);
            end
        end
    

        %plot psth at the top
        spidxA = sub2ind([length(trialLists) total_length], a, 1);
        spidxB = sub2ind([length(trialLists) total_length], a, psth);
        subplot(total_length,length(trialLists),[spidxA spidxB])
        [meanPSTH, semPSTH, rasters] = computePSTHs(alignedResps{ETAs(a)}(:,:,plotCells(k)),trialArrays{a});
        if strcmp(pickTrials{a},'direction') && ~strcmp(pickTrials{a},'_direction')
            meanPSTH = meanPSTH';
            semPSTH = semPSTH';
            rasters = rasters';
        end
        for d = 1:size(meanPSTH,2)
            if d == 1
                ls = '-';
            else
                ls = ':';
            end
            plotPSTHs(eventWindow, cell2mat(meanPSTH(:,d)), cell2mat(semPSTH(:,d)), psthColors{a},ls);
            if a > 1
                xlim([-.5 2]);
            else
                xlim([-.5 2]);
            end
            yMin(end+1) = min(min(cell2mat(meanPSTH)-cell2mat(semPSTH)));
            yMax(end+1) = max([.01 max(max(cell2mat(meanPSTH)+cell2mat(semPSTH)))]);
        end
        
        ln = line([0 0],[-1 10],'LineStyle',':','Color',[0 0 0],'linewidth',.5);
        uistack(ln,'bottom');
        if a == 1
            lv = medianTTG;
            l1 = line([lv lv],[-1 10],'LineStyle',':','Color',[0 0 0],'linewidth',.5);
        uistack(l1,'bottom');
        elseif a == 2
            lv = medianTTM;
            l1 = line([lv lv],[-1 10],'LineStyle',':','Color',[0 0 0],'linewidth',.5);
        uistack(l1,'bottom');
        elseif a == 3
            lv = medianTTF;
            l1 = line([lv lv],[-1 10],'LineStyle',':','Color',[0 0 0],'linewidth',.5);
        uistack(l1,'bottom');
        end
        
        
        set(gca,'xtick',[]);
        set(gca,'TickDir','out');
    end
    
    %unify axis limits across plots and add labels as appropriate
    for a = 1:length(ETAs)
        for iCond = 1:length(trialLists{a})
            [spidx1, spidx2] = goToRasterSubplot(length(trialLists), total_length, cellfun(@length, trialLists{a})', a, iCond, 1, psth);
            subplot(total_length,length(trialLists),[spidx1 spidx2])
            caxis([min(iMin) max(iMax)*.5]);
            if iCond == length(trialLists{a})
                if a == 1
                    xlabel('Time – stimOn (s)')
                elseif a == 2
                    xlabel('Time – cueOn (s)')
                elseif a == 3
                    xlabel('Time – moveOn (s)')
                elseif a == 4
                    xlabel('Time – rewOn (s)')
                end
            end
        end
        spidxA = sub2ind([length(trialLists) total_length], a, 1);
        spidxB = sub2ind([length(trialLists) total_length], a, psth);    
        subplot(total_length,length(trialLists),[spidxA spidxB])
        box off
        if ~exist('yOverride')
            ymx = ceil(max(yMax)*10)/10;
        else
            ymx = yOverride;
        end
        ymn = 0;
        ymd = ymx/2;
        ylim([ymn ymx]);
        set(gca,'ytick',[ymd ymx]);
        if a == 1 
            ylabel('Activity')
        end
        if a == 1
            for t = 1:length(densityTTG)
                fill([bins(t) bins(t+1) bins(t+1) bins(t)],[ymx*.95 ymx*.95 ymx ymx],'k','EdgeColor','none','FaceAlpha',densityTTG(t));
            end
            [~, mi] = max(densityTTG);
%             txt = strcat(num2str(medianTTF*1000),{' ms'});
%             sprintf('Line #1\nThe second line.'
            text(bins(mi),ymx*1.03,strcat(num2str(medianTTG)),'HorizontalAlignment', 'center','FontSize',8);
         end
         if a == 2
            for t = 1:length(densityTTM)
                fill([bins(t) bins(t+1) bins(t+1) bins(t)],[ymx*.95 ymx*.95 ymx ymx],'k','EdgeColor','none','FaceAlpha',densityTTM(t));
            end
            [~, mi] = max(densityTTM);
%             txt = strcat(num2str(medianTTF*1000),{' ms'});
%             sprintf('Line #1\nThe second line.'
            text(bins(mi),ymx*1.03,strcat(num2str(medianTTM)),'HorizontalAlignment', 'center','FontSize',8);
        end
        if a == 3
            for t = 1:length(densityTTF)
                fill([bins(t) bins(t+1) bins(t+1) bins(t)],[ymx*.95 ymx*.95 ymx ymx],'k','EdgeColor','none','FaceAlpha',densityTTF(t));
            end
            [~, mi] = max(densityTTF);
%             txt = strcat(num2str(medianTTF*1000),{' ms'});
%             sprintf('Line #1\nThe second line.'
            text(bins(mi),ymx*1.03,strcat(num2str(medianTTF)),'HorizontalAlignment', 'center','FontSize',8);
        end
%         if a == 3
%             for t = 1:length(densityTTF)
%                 fill([-bins(t) -bins(t+1) -bins(t+1) -bins(t)],[max(yMax)*.95 max(yMax)*.95 max(yMax) max(yMax)],'k','EdgeColor','none','FaceAlpha',densityTTF(t));
%             end
%         end

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
    elseif was_a_key && strcmp(get(fig, 'CurrentKey'), 'escape')
        disp(k)
        close(gcf)
        break
    end
end
