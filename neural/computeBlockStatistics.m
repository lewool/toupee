for m = 4:length(mouseList)
    
    %% IMPORT DATA
    mouseName = char(mouseList{m});
    expDate = char(expList{m}{1});
    expNum = expList{m}{2};
    hemisphere = hemList(m);
    [expRef, ~] = data.constructExpRef(mouseName,expDate,expNum);
    fprintf('fetching %s...\n',expRef)
    [behavioralData, expInfo, neuralData] = data.loadDataset(mouseName, expDate, expNum);
    [baselineResps, stimResps, pmovResps, movResps, rewResps, preCueResps] = getEpochResps(neuralData.eta);
    
    %% SET PARAMS
    trimLength = 10;
    nt = length(behavioralData.eventTimes(1).daqTime);
    np = 1000;
    nc = size(baselineResps,2);
    
    [~, whichTrials] = selectCondition(expInfo, getUniqueContrasts(expInfo), behavioralData, ...
        initTrialConditions('responseType','all','movementTime','all','specificRTs',[0.8 5]));

    %% EXTRACT TRIAL VARIABLES AND GENERATE PSEUDOBLOCKS
    
    % correct for hemisphere
    if hemisphere < 0
        trueStimuli = expInfo.block.events.contrastValues(1:nt);
        trueCorrectChoice = expInfo.block.events.correctResponseValues(1:nt);
        trueBlocks = expInfo.block.events.highRewardSideValues(1:nt);
        trueActualChoice = expInfo.block.events.responseValues(1:nt);
    else
        trueStimuli = -expInfo.block.events.contrastValues(1:nt);
        trueCorrectChoice = -expInfo.block.events.correctResponseValues(1:nt);
        trueBlocks = -expInfo.block.events.highRewardSideValues(1:nt);
        trueActualChoice = -expInfo.block.events.responseValues(1:nt);
    end
    trueFeedback = expInfo.block.events.feedbackValues(1:nt);
    
    % remove the first trials after the block switch by overwriting nans
    idxTrim = false(1,nt);
    switchPoints = find(diff(expInfo.block.events.highRewardSideValues(1:nt)) ~= 0);
    if numel(switchPoints) == 1
        idxTrim(trimLength+1:switchPoints) = true;
        idxTrim(switchPoints+trimLength+1:nt) = true;
    else
        for s = 1:length(switchPoints)
            if s == 1
                idxTrim(trimLength+1:switchPoints(s)) = true;
            elseif s < length(switchPoints)
                idxTrim(switchPoints(s-1)+trimLength+1:switchPoints(s)) = true;
            elseif s == length(switchPoints)
                idxTrim(switchPoints(s-1)+trimLength+1:switchPoints(s)) = true;
                idxTrim(switchPoints(s)+trimLength+1:nt) = true;
            end
        end
    end
    trueBlocks(~idxTrim) = nan;
    
    % assign the 0% stimuli as either 'left' or 'right' depending on where
    % the reward was preassigned (contra or ipsi)
    trueStimuli(trueStimuli == 0) = eps;
    trueStimuli(abs(trueStimuli) < .05) = ...
        trueStimuli(abs(trueStimuli) < .05).* trueCorrectChoice(abs(trueStimuli) < .05);

    % low rewards are possible on sign-mismatched block and stimulus
    % high rewards are possible on sign-matched block and stimulus
    % 2 = high, 1 = low, 0 = error/none, nan = omitted trim trials
    trueValue(trueBlocks.*sign(trueStimuli) == -1) = 1;
    trueValue(trueBlocks.*sign(trueStimuli) == 1) = 2;
    trueValue(trueFeedback==0) = 0;
    trueValue(~idxTrim) = nan;

    % generate pseudoblocks
    blockStart = 'fixed';
    for p = 1:np
        if strcmp(blockStart,'fixed')
            firstSide = trueBlocks(trimLength+1);
        elseif strcmp(blockStart,'rand')
            firstSide = randsample([-1, 1],1,true);
        end
        b=nan(1,nt);
        switches = cumsum(125+randi(100,1,20));
        for s = 1:length(switches)
            if s == 1
                b((1+trimLength):switches(s)-1) = firstSide;
            elseif mod(s,2) == 1
                b((switches(s-1)+trimLength):switches(s)-1) = firstSide;
            elseif mod(s,2) == 0
                b((switches(s-1)+trimLength):switches(s)-1) = -firstSide;
            end
        end
        pseudoBlocks(:,p) = b(1:nt);
    end

    % assign new high/low trials based on the pseudoblocks + true stim
    for p = 1:np
        pseudoValue(pseudoBlocks(:,p).*sign(trueStimuli)' == -1,p) = 1;
        pseudoValue(pseudoBlocks(:,p).*sign(trueStimuli)' == 1,p) = 2;
        pseudoValue(trueFeedback==0,p) = 0;
        pseudoValue(isnan(pseudoBlocks(:,p)),p) = nan;
    end
    
    %% SET UP TRIALS TO COMPARE
    
    %1. contra vs ipsi blocks, all trials
    trialLabels{1} = 'block, all trials';
    trialSets{1,1} = intersect(find(trueBlocks > 0),whichTrials);
    trialSets{1,2} = intersect(find(trueBlocks < 0),whichTrials);
    
    %2. contra vs ipsi blocks, hcc trials
    trialLabels{2} = 'block, hcc trials';
    trialSets{2,1} = intersect(find(trueStimuli > .2 & trueBlocks > 0),whichTrials);
    trialSets{2,2} = intersect(find(trueStimuli > .2 & trueBlocks < 0),whichTrials);
    
    %3. contra vs ipsi blocks, lcc trials
    trialLabels{3} = 'block, lcc trials';
    trialSets{3,1} = intersect(find(trueStimuli < .2 & trueStimuli > .01 & trueBlocks > 0),whichTrials);
    trialSets{3,2} = intersect(find(trueStimuli < .2 & trueStimuli > .01 & trueBlocks < 0),whichTrials);
    
    %4. contra vs ipsi blocks, hci trials
    trialLabels{4} = 'block, hci trials';
    trialSets{4,1} = intersect(find(trueStimuli < -.2 & trueBlocks > 0),whichTrials);
    trialSets{4,2} = intersect(find(trueStimuli < -.2 & trueBlocks < 0),whichTrials);
    
    %5. contra vs ipsi blocks, lci trials
    trialLabels{5} = 'block, lci trials';
    trialSets{5,1} = intersect(find(trueStimuli > -.2 & trueStimuli < -.01 & trueBlocks > 0),whichTrials);
    trialSets{5,2} = intersect(find(trueStimuli > -.2 & trueStimuli < -.01 & trueBlocks < 0),whichTrials);
    
    %6. contra vs ipsi blocks, 0% trials
    trialLabels{6} = 'block, 0% trials';
    trialSets{6,1} = intersect(find(abs(trueStimuli) < .01 & trueBlocks > 0),whichTrials);
    trialSets{6,2} = intersect(find(abs(trueStimuli) < .01 & trueBlocks < 0),whichTrials);
    
    %7. high vs low value, all hc trials
    trialLabels{7} = 'value, hc trials';
    trialSets{7,1} = intersect(find(abs(trueStimuli) > .1 & trueValue == 2),whichTrials);
    trialSets{7,2} = intersect(find(abs(trueStimuli) > .1 & trueValue == 1),whichTrials);
    
    %8. high vs low value, all lc trials
    trialLabels{8} = 'value, lc trials';
    trialSets{8,1} = intersect(find(abs(trueStimuli) < .2 & abs(trueStimuli) > .01 & trueValue == 2),whichTrials);
    trialSets{8,2} = intersect(find(abs(trueStimuli) < .2 & abs(trueStimuli) > .01 & trueValue == 1),whichTrials);
    
    %9. high vs low value, 0% trials
    trialLabels{9} = 'value, 0% trials';
    trialSets{9,1} = intersect(find(abs(trueStimuli) < .01 & trueValue == 2),whichTrials);
    trialSets{9,2} = intersect(find(abs(trueStimuli) < .01 & trueValue == 1),whichTrials);
    
    trialSets_pseudo = cell(length(trialLabels),2,np);
    for p = 1:np
        %contra vs ipsi block, all trials
        trialSets_pseudo{1,1,p} = intersect(find(pseudoBlocks(:,p) > 0),whichTrials);
        trialSets_pseudo{1,2,p} = intersect(find(pseudoBlocks(:,p) < 0),whichTrials);
        
        %contra vs ipsi block, hc contralateral stim
        trialSets_pseudo{2,1,p} = intersect(find(trueStimuli > .2 & pseudoBlocks(:,p)' > 0),whichTrials);
        trialSets_pseudo{2,2,p} = intersect(find(trueStimuli > .2 & pseudoBlocks(:,p)' < 0),whichTrials);
        
        %contra vs ipsi block, lc contralateral stim
        trialSets_pseudo{3,1,p} = intersect(find(trueStimuli < .2 & trueStimuli > .01 & pseudoBlocks(:,p)' > 0),whichTrials);
        trialSets_pseudo{3,2,p} = intersect(find(trueStimuli < .2 & trueStimuli > .01 & pseudoBlocks(:,p)' < 0),whichTrials);
        
        %contra vs ipsi block, hc ipsilateral stim
        trialSets_pseudo{4,1,p} = intersect(find(trueStimuli < -.2 & pseudoBlocks(:,p)' > 0),whichTrials);
        trialSets_pseudo{4,2,p} = intersect(find(trueStimuli < -.2 & pseudoBlocks(:,p)' < 0),whichTrials);
        
        %contra vs ipsi block, lc ipsilateral stim
        trialSets_pseudo{5,1,p} = intersect(find(trueStimuli > -.2 & trueStimuli < -.01 & pseudoBlocks(:,p)' > 0),whichTrials);
        trialSets_pseudo{5,2,p} = intersect(find(trueStimuli > -.2 & trueStimuli < -.01 & pseudoBlocks(:,p)' < 0),whichTrials);
        
        %contra vs ipsi, 0% stim
        trialSets_pseudo{6,1,p} = intersect(find(abs(trueStimuli) < .01 & pseudoBlocks(:,p)' > 0),whichTrials);
        trialSets_pseudo{6,2,p} = intersect(find(abs(trueStimuli) < .01 & pseudoBlocks(:,p)' < 0),whichTrials);
        
        % high vs low value, all high-contrast stim
        trialSets_pseudo{7,1,p} = intersect(find(abs(trueStimuli) > .2 & pseudoValue(:,p)' == 2),whichTrials);
        trialSets_pseudo{7,2,p} = intersect(find(abs(trueStimuli) > .2 & pseudoValue(:,p)' == 1),whichTrials);
        
        %high vs low value, all low-contrast stim
        trialSets_pseudo{8,1,p} = intersect(find(abs(trueStimuli) < .2 & abs(trueStimuli) > .01 & pseudoValue(:,p)' == 2),whichTrials);
        trialSets_pseudo{8,2,p} = intersect(find(trueStimuli < .2 & abs(trueStimuli) > .01 & pseudoValue(:,p)' == 1),whichTrials);
        
        %high vs low value, 0% stim
        trialSets_pseudo{9,1,p} = intersect(find(abs(trueStimuli) < .01 & pseudoValue(:,p)' == 2),whichTrials);
        trialSets_pseudo{9,2,p} = intersect(find(abs(trueStimuli) < .01 & pseudoValue(:,p)' == 1),whichTrials);
    end
    
    %% COMPUTE 'BLOCK METRIC' FOR EACH CELL, TRUE & PSEUDO
    allResps(:,:,1) = baselineResps;
    allResps(:,:,2) = stimResps;
    allResps(:,:,3) = movResps;
    allResps(:,:,4) = rewResps;

    respLabels = {'baselineResps' 'stimResps' 'movResps' 'rewResps'};
    r=0;
    fprintf(1,'computing metric...\n')
    
    for r = 1:size(allResps,3)
        fprintf(1,'\b%3.0f',r);
        whichResps = allResps(:,:,r);
        metrics_true = nan(nc,length(trialLabels));
        for t = 1:length(trialLabels)
            for c = 1:nc
                respsA = whichResps(trialSets{t,1},c);
                respsB = whichResps(trialSets{t,2},c);
                metrics_true(c,t) = (nanmean(respsA) - nanmean(respsB)) / (nanmean(respsA) + nanmean(respsB) + eps);
            end
        end

        metrics_pseudo = nan(nc,length(trialLabels),np);
        for t = 1:length(trialLabels)
            for c = 1:nc
                for p = 1:np
                    respsA = whichResps(trialSets_pseudo{t,1,p},c);
                    respsB = whichResps(trialSets_pseudo{t,2,p},c);
                    metrics_pseudo(c,t,p) = (nanmean(respsA) - nanmean(respsB)) / (nanmean(respsA) + nanmean(respsB) + eps);
                end
            end
        end     
    
        %% COMPUTE TEST STATISTICS, TRUE AND PSEUDO

        thresh = 0.2;
        tsLabels = {'Mean(abs)' 'Proportion > 0' 'Proportion > thresh' 'Proportion (abs) > thresh'};
        testStatistics_true = [...
            nanmean(abs(metrics_true));
            sum(metrics_true > 0)./nc;
            sum(metrics_true > thresh)./nc;
            sum(abs(metrics_true) > thresh)./nc;
            ];

        testStatistics_pseudo(:,:,1) = squeeze(nanmean(abs(metrics_pseudo)))';
        testStatistics_pseudo(:,:,2) = squeeze(sum(metrics_pseudo > 0) ./ nc)';
        testStatistics_pseudo(:,:,3) = squeeze(sum(metrics_pseudo > thresh) ./ nc)';
        testStatistics_pseudo(:,:,4) = squeeze(sum(abs(metrics_pseudo) > thresh) ./ nc)';
        
        %capture the sigrank for each metric/ts against the pseudo distribution
        for t = 1:length(trialLabels)
            for ts = 1:length(tsLabels)
                [~, idx] = min(sort(abs(testStatistics_pseudo(:,t,ts) - testStatistics_true(ts,t))));
                sigrank(ts,t,r,m) = (sum(testStatistics_pseudo(:,t,ts) < testStatistics_true(ts,t)))/np;
            end
        end
    
        %% PLOT

%         sdim = [length(tsLabels)+1 length(trialLabels)];
% 
%         figure;
%         set(gcf,'position',[450 555 2350 1070]);
%         for t = 1:length(trialLabels)
%             subplot(sdim(1),sdim(2),t)
%             hold on;
%             h = histogram(metrics_true(:,t),41,'FaceColor',[.5 0 0]);
%             maxx = max(h.BinEdges);
%             maxy = max(h.Values);
%             box off
%             set(gca,'tickdir','out')
%             xlim([-maxx maxx])
%             ylim([0 maxy*1.1])
%             title(trialLabels{t})
%             if t == 1
%                 ylabel('No. neurons')
%                 bigTitle = text(-maxx,maxy*1.6,strcat(expRef,{': '},respLabels{r}),'Color',[0 0 0],'FontSize',15,'FontWeight','bold','HorizontalAlignment','left','Interpreter','none');
%             end
%         end
% 
%         for t = 1:length(trialLabels)
%             for ts = 1:size(testStatistics_true,1)
%                 subplot(sdim(1),sdim(2),length(trialLabels)*ts+t)
%                 hold on;
%                 h = histogram(testStatistics_pseudo(:,t,ts),41,'FaceColor',[.5 .5 .5]);
%                 if testStatistics_true(ts,t) > prctile(testStatistics_pseudo(:,t,ts),97.5) || testStatistics_true(ts,t) < prctile(testStatistics_pseudo(:,t,ts),2.5)
%                     line([testStatistics_true(ts,t) testStatistics_true(ts,t)],[0 400],'LineStyle','-','LineWidth',3,'Color','r');
%                 else
%                     line([testStatistics_true(ts,t) testStatistics_true(ts,t)],[0 400],'LineStyle','-','LineWidth',3,'Color',[.5 0 0]);
%                 end
%                 maxx = max(h.BinEdges);
%                 maxy = max(h.Values);
%                 box off
%                 set(gca,'tickdir','out')
%                 ylim([0 maxy*1.1])
%                 if t == 1
%                     ylabel(tsLabels{ts})
%                 end
%             end
%         end
% 
%         printfig(gcf,char(strcat(expRef,{' block statistics - '},respLabels{r})))
%         close all
    end
    %%
    clearvars -except mouseList expList hemList m sigrank
    
end
    
    

%%

for t = 1:9
    for ts = 1:4
        baselineTests{ts,t} = squeeze(sigrank(ts,t,1,:));
    end
end

%%

for t = 1:9
    for ts = 1:4
        stimTests{ts,t} = squeeze(sigrank(ts,t,2,:));
    end
end

%%

for t = 1:9
    for ts = 1:4
        movTests{ts,t} = squeeze(sigrank(ts,t,3,:));
    end
end

%%

for t = 1:9
    for ts = 1:4
        rewTests{ts,t} = squeeze(sigrank(ts,t,4,:));
    end
end

%% plot
figure;
set(gcf,'position',[1250 825 1270 800])
for t = 1:9
    for ts = 1:4
        subplot(4,9,9*(ts-1)+t)
        scatter(rand(1,34),baselineTests{ts,t},'MarkerFaceAlpha',.5,'MarkerFaceColor','k','MarkerEdgeColor','none')
        ylim([0 1])
        xlim([-.5 1.5])
        xticks({})
    end
end

%%
figure;
set(gcf,'position',[1250 825 1270 800])
for t = 1:9
    for ts = 1:4
        subplot(4,9,9*(ts-1)+t)
        scatter(rand(1,34),stimTests{ts,t},'MarkerFaceAlpha',.5,'MarkerFaceColor','k','MarkerEdgeColor','none')
        ylim([0 1])
        xlim([-.5 1.5])
        xticks({})
    end
end

%%
figure;
set(gcf,'position',[1250 825 1270 800])
for t = 1:9
    for ts = 1:4
        subplot(4,9,9*(ts-1)+t)
        scatter(rand(1,34),movTests{ts,t},'MarkerFaceAlpha',.5,'MarkerFaceColor','k','MarkerEdgeColor','none')
        ylim([0 1])
        xlim([-.5 1.5])
        xticks({})
    end
end

%%
figure;
set(gcf,'position',[1250 825 1270 800])
for t = 1:9
    for ts = 1:4
        subplot(4,9,9*(ts-1)+t)
        scatter(rand(1,34),rewTests{ts,t},'MarkerFaceAlpha',.5,'MarkerFaceColor','k','MarkerEdgeColor','none')
        ylim([0 1])
        xlim([-.5 1.5])
        xticks({})
    end
end
    
    
    
    
    
  