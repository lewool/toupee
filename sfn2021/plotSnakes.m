%% select a subset of trial types to visualize
contrasts = getUniqueContrasts(expInfo);
[~, allTrials] = selectCondition(expInfo, contrasts, behavioralData, ...
    initTrialConditions('movementTime','late'));

[~, mL] = selectCondition(expInfo, contrasts, behavioralData, ...
    initTrialConditions('movementDir','cw','movementTime','late'));
[~, mR] = selectCondition(expInfo, contrasts, behavioralData, ...
    initTrialConditions('movementDir','ccw','movementTime','late'));

[~, sL] = selectCondition(expInfo, contrasts(contrasts<0), behavioralData, ...
    initTrialConditions('movementTime','late'));
[~, sR] = selectCondition(expInfo, contrasts(contrasts>0), behavioralData, ...
    initTrialConditions('movementTime','late'));

%% STIMULUS: same thing, but aligned to stimulus

% compute each cell's mean event-aligned activity (even trials vs odd trials)
stimResponses = neuralData.eta.alignedResps{1}(allTrials,:,:);
evenStimTrials = stimResponses(2:2:end,:,:);
oddStimTrials = stimResponses(1:2:end,:,:);
odd_sL = neuralData.eta.alignedResps{1}(intersect(allTrials(logical(mod(allTrials,2))), sL),:,:);
odd_sR = neuralData.eta.alignedResps{1}(intersect(allTrials(logical(mod(allTrials,2))), sR),:,:);

stimEven = squeeze(nanmean(evenStimTrials,1))';
stimOdd = squeeze(nanmean(oddStimTrials,1))';
stimOdd_sL = squeeze(nanmean(odd_sL,1))';
stimOdd_sR = squeeze(nanmean(odd_sR,1))';

% %normalize each cell's activity to its own min/max (0 - 1)
% normStimEven = (stimEven - min(stimEven,[],2)) ./ max(stimEven - min(stimEven,[],2),[],2);
% normStimOdd = (stimOdd - min(stimOdd,[],2)) ./ max(stimOdd - min(stimOdd,[],2),[],2);
% normStimOdd_sL = (stimOdd_sL - min(stimOdd_sL,[],2)) ./ max(stimOdd_sL - min(stimOdd_sL,[],2),[],2);
% normStimOdd_sR = (stimOdd_sR - min(stimOdd_sR,[],2)) ./ max(stimOdd_sR - min(stimOdd_sR,[],2),[],2);

%normalize each cell's activity to its own min/max (0 - 1)
stimEven_baseline = mean(baselineResps(2:2:end,:),1);
stimOdd_baseline = mean(baselineResps(1:2:end,:),1);
odd_sL_baseline = mean(baselineResps(intersect(allTrials(logical(mod(allTrials,2))), sL),:),1);
odd_sR_baseline = mean(baselineResps(intersect(allTrials(logical(mod(allTrials,2))), sR),:),1);

%normalize each cell's activity to its own min/max (0 - 1)
normStimEven = stimEven - stimEven_baseline';
normStimOdd = stimOdd - stimOdd_baseline';
normStimOdd_sL = stimOdd_sL - odd_sL_baseline';
normStimOdd_sR = stimOdd_sR - odd_sR_baseline';

% sort the activity
%timestamp of the max of each cell's trace
[~, maxIdx] = max(normStimEven(:,15:41),[],2);

%sort the timestamps by early to late
[~,sortIdx] = sort(maxIdx);

%%
gs = 4;
g = .4;

figure;
set(gcf,'position',[680   103   560   995]);
ax1 = subplot(1,3,1);
imagesc(...
    neuralData.eta.eventWindow, ...
    1:length(neuralData.eta.alignedResps{1}), ...
    smoothdata(normStimOdd_sL(sortIdx,:),2,'gaussian',gs));
xlim([-0.5 2]);
caxis([-0.12 .5])
box off
set(gca,'tickdir','out')
line([0 0],[1 length(neuralData.eta.alignedResps{1})],'LineStyle','--','Color','k');
colormap(ax1, flipud(gray));

ax2 = subplot(1,3,2);
imagesc(...
    neuralData.eta.eventWindow, ...
    1:length(neuralData.eta.alignedResps{1}), ...
    smoothdata(normStimOdd_sR(sortIdx,:),2,'gaussian',gs));
xlim([-0.5 2]);
caxis([-0.12 .5])
box off
set(gca,'tickdir','out')
line([0 0],[1 length(neuralData.eta.alignedResps{1})],'LineStyle','--','Color','k');
colormap(ax2, flipud(gray));
 
ax3 = subplot(1,3,3);
imagesc(...
    neuralData.eta.eventWindow, ...
    1:length(neuralData.eta.alignedResps{1}), ...
    smoothdata(normStimOdd_sR(sortIdx,:) - normStimOdd_sL(sortIdx,:),2,'gaussian',gs));
xlim([-0.5 2]);
caxis([-.5 .5 ])
box off
set(gca,'tickdir','out')
line([0 0],[1 length(neuralData.eta.alignedResps{1})],'LineStyle','--','Color','k');
colormap(ax3, BlueWhiteRed(100,g))
% colormap(ax3, flipud(gray));

printfig(gcf,'snakes test')
% close all

%% MOVEMENT

[baselineResps, stimResps, pmovResps, movResps, rewResps] = getEpochResps(neuralData.eta);

% compute each cell's mean event-aligned activity (even trials vs odd trials)
movResponses = neuralData.eta.alignedResps{2}(allTrials,:,:);
evenMoveTrials = movResponses(2:2:end,:,:);
oddMoveTrials = movResponses(1:2:end,:,:);
odd_mL = neuralData.eta.alignedResps{2}(intersect(allTrials(logical(mod(allTrials,2))), mL),:,:);
odd_mR = neuralData.eta.alignedResps{2}(intersect(allTrials(logical(mod(allTrials,2))), mR),:,:);

moveEven = squeeze(nanmean(evenMoveTrials,1))';
moveOdd = squeeze(nanmean(oddMoveTrials,1))';
moveOdd_mL = squeeze(nanmean(odd_mL,1))';
moveOdd_mR = squeeze(nanmean(odd_mR,1))';

% %normalize each cell's activity to its own min/max (0 - 1)
% normMoveEven = (moveEven - min(moveEven,[],2)) ./ max(moveEven - min(moveEven,[],2),[],2);
% normMoveOdd = (moveOdd - min(moveOdd,[],2)) ./ max(moveOdd - min(moveOdd,[],2),[],2);
% normMoveOdd_mL = (moveOdd_mL - min(moveOdd_mL,[],2)) ./ max(moveOdd_mL - min(moveOdd_mL,[],2),[],2);
% normMoveOdd_mR = (moveOdd_mR - min(moveOdd_mR,[],2)) ./ max(moveOdd_mR - min(moveOdd_mR,[],2),[],2);

%normalize each cell's activity to its own min/max (0 - 1)
moveEven_baseline = mean(baselineResps(2:2:end,:),1);
moveOdd_baseline = mean(baselineResps(1:2:end,:),1);
odd_mL_baseline = mean(baselineResps(intersect(allTrials(logical(mod(allTrials,2))), mL),:),1);
odd_mR_baseline = mean(baselineResps(intersect(allTrials(logical(mod(allTrials,2))), mR),:),1);

normMoveEven = moveEven - moveEven_baseline';
normMoveOdd = moveOdd - moveOdd_baseline';
normMoveOdd_mL = moveOdd_mL - odd_mL_baseline';
normMoveOdd_mR = moveOdd_mR - odd_mR_baseline';

% % sort the activity
% %timestamp of the max of each cell's trace
% [~, maxIdx] = max(normStimEven,[],2);
% 
% %sort the timestamps by early to late
% [~,sortIdx] = sort(maxIdx);


%%
gs = 4;
g = .4;

figure;
set(gcf,'position',[680   103   560   995]);
ax1 = subplot(1,3,1);
imagesc(...
    neuralData.eta.eventWindow, ...
    1:length(neuralData.eta.alignedResps{2}), ...
    smoothdata(normMoveOdd_mL(sortIdx,:),2,'gaussian',gs));
xlim([-1 2]);
caxis([-.12 .6])
box off
set(gca,'tickdir','out')
line([0 0],[1 length(neuralData.eta.alignedResps{2})],'LineStyle','--','Color','k');
colormap(ax1, flipud(gray));

ax2 = subplot(1,3,2);
imagesc(...
    neuralData.eta.eventWindow, ...
    1:length(neuralData.eta.alignedResps{2}), ...
    smoothdata(normMoveOdd_mR(sortIdx,:),2,'gaussian',gs));
xlim([-1 2]);
caxis([-.12 .6])
box off
set(gca,'tickdir','out')
line([0 0],[1 length(neuralData.eta.alignedResps{2})],'LineStyle','--','Color','k');
colormap(ax2, flipud(gray));

ax3 = subplot(1,3,3);
imagesc(...
    neuralData.eta.eventWindow, ...
    1:length(neuralData.eta.alignedResps{1}), ...
    smoothdata(normMoveOdd_mR(sortIdx,:) - normMoveOdd_mL(sortIdx,:),2,'gaussian',gs));
xlim([-1 2]);
caxis([-.6 .6 ])
box off
set(gca,'tickdir','out')
line([0 0],[1 length(neuralData.eta.alignedResps{1})],'LineStyle','--','Color','k');
colormap(ax3, BlueWhiteRed(100,g))

printfig(gcf,'snakes test')
% close all