%% make a list of mice/experiments you want to analyze

mouseList = { ...
    {'LEW031'}...
    };

expList = { ...
    {'2020-02-13',1,[1]}};

%% load all the experiments into expInfo

expInfo = initExpInfo(mouseList,expList);

%% process the usual data
% the script knows to loop over all the experiments you listed above
% this will take a while but the command line will print progress

[expInfo, neuralData, behavioralData] = processExperiment(expInfo);

%% select a subset of trial types to visualize
contrasts = getUniqueContrasts(expInfo);
[~, whichTrials] = selectCondition(expInfo, contrasts, behavioralData, ...
    initTrialConditions('movementTime','late'));

[~, mL] = selectCondition(expInfo, contrasts, behavioralData, ...
    initTrialConditions('movementDir','cw','movementTime','late'));
[~, mR] = selectCondition(expInfo, contrasts, behavioralData, ...
    initTrialConditions('movementDir','ccw','movementTime','late'));

[~, sL] = selectCondition(expInfo, contrasts(contrasts<0), behavioralData, ...
    initTrialConditions('movementTime','late'));
[~, sR] = selectCondition(expInfo, contrasts(contrasts>0), behavioralData, ...
    initTrialConditions('movementTime','late'));

%% 1. MOVEMENT: get cell-by-cell mean activity (aligned to first movement)
% here, we compute each cell's mean event-aligned activity 
% and sort them by time of max activity. we do this in a cross-validated way
% by splitting the trials in half arbitrarily (evens vs odds), and use one
% half to sort the times (train, evens), and the other half to display
% (test, odds) using the sorting index we determined from the evens

% compute each cell's mean event-aligned activity (even trials vs odd trials)
moveResponses = neuralData.eta.alignedResps{2}(whichTrials,:,:);
evenMoveTrials = moveResponses(2:2:end,:,:);
oddMoveTrials = moveResponses(1:2:end,:,:);

moveEven = squeeze(nanmean(evenMoveTrials,1))';
moveOdd = squeeze(nanmean(oddMoveTrials,1))';

%normalize each cell's activity to its own min/max (0 - 1). (we can't
%really say much about relative activity between 2p ROIs anyway, so we
%might as well optimize each so we can see the activity easily in a plot)
normMoveEven = (moveEven - min(moveEven,[],2)) ./ max(moveEven - min(moveEven,[],2),[],2);
normMoveOdd = (moveOdd - min(moveOdd,[],2)) ./ max(moveOdd - min(moveOdd,[],2),[],2);

% sort the activity
%timestamp of the max of each cell's trace
[~, maxIdx] = max(normMoveEven,[],2);

%sort the timestamps by early to late
[~,sortIdx] = sort(maxIdx);
 
%% 2. STIMULUS: same thing, but aligned to stimulus

% compute each cell's mean event-aligned activity (even trials vs odd trials)
stimResponses = neuralData.eta.alignedResps{1}(whichTrials,:,:);
evenStimTrials = stimResponses(2:2:end,:,:);
oddStimTrials = stimResponses(1:2:end,:,:);

stimEven = squeeze(nanmean(evenStimTrials,1))';
stimOdd = squeeze(nanmean(oddStimTrials,1))';

%normalize each cell's activity to its own min/max (0 - 1)
normStimEven = (stimEven - min(stimEven,[],2)) ./ max(stimEven - min(stimEven,[],2),[],2);
normStimOdd = (stimOdd - min(stimOdd,[],2)) ./ max(stimOdd - min(stimOdd,[],2),[],2);

% FYI we will use the same sorting as we found for the movement, so no need to
% recompute here

%% 3. OUTCOME: same thing, but aligned to feedback

% compute each cell's mean event-aligned activity (even trials vs odd trials)
rewResponses = neuralData.eta.alignedResps{3}(whichTrials,:,:);
evenRewTrials = rewResponses(2:2:end,:,:);
oddRewTrials = rewResponses(1:2:end,:,:);

rewEven = squeeze(nanmean(evenRewTrials,1))';
rewOdd = squeeze(nanmean(oddRewTrials,1))';

%normalize each cell's activity to its own min/max (0 - 1)
normRewEven = (rewEven - min(rewEven,[],2)) ./ max(rewEven - min(rewEven,[],2),[],2);
normRewOdd = (rewOdd - min(rewOdd,[],2)) ./ max(rewOdd - min(rewOdd,[],2),[],2);

% FYI we will use the same sorting as we found for the movement, so no need to
% recompute here

% normalize across the three alignments
smr(:,:,1) = normStimOdd;
smr(:,:,2) = normMoveOdd;
smr(:,:,3) = normRewOdd;

for c = 1:length(smr)
    cellSlice = squeeze(smr(c,:,:));
    sliceMin = min(min(cellSlice));
    normSlice = cellSlice - sliceMin;
    sliceMax = max(max(normSlice));
    smr_norm (c,:,:) = normSlice / sliceMax;
end

%% 4. PLOT
figure;
try 
    cm = BlueWhiteRed(100,.3);
    cm = flipud(gray);
catch
    cm = flipud(gray);
end

subplot(1,9,1)
line([1 1],[length(normMoveOdd)-500 length(normMoveOdd)],'LineWidth',3,'Color','k')
xlim([0 1])
ylim([1 length(normMoveOdd)])
h=text(.8,length(normMoveOdd)-500,'  500 cells');
set(h,'Rotation',90);
axis off

subplot(1,9,[2 3])
imagesc(...
    neuralData.eta.eventWindow,...
    1:length(normStimOdd),...
    smoothdata(smr_norm(sortIdx,:,1),2,'gaussian',5))
colormap(cm)
stimLine = line([0 0],[1 length(normMoveOdd)]);
set(stimLine, 'Color','k','LineStyle','--');
xlim([-0.5 1.5]);
box off;
set(gca,'tickdir','out')
xlabel('Time (s)')
title('Stimulus')
set(gca,'ytick',[])
set(gca,'xtick',[-1 0 1])

subplot(1,9,[5 6]);
imagesc(...
    neuralData.eta.eventWindow,...
    1:length(normMoveOdd),...
    smoothdata(smr_norm(sortIdx,:,2),2,'gaussian',5))
colormap(cm)
stimLine = line([0 0],[1 length(normMoveOdd)]);
set(stimLine, 'Color','k','LineStyle','--');
xlim([-1 1]);
box off;
set(gca,'tickdir','out')
title('Movement')
set(gca,'ytick',[])
set(gca,'xtick',[-1 0 1])

subplot(1,9,[8 9])
imagesc(...
    neuralData.eta.eventWindow,...
    1:length(normRewOdd),...
    smoothdata(smr_norm(sortIdx,:,3),2,'gaussian',5))
colormap(cm)
stimLine = line([0 0],[1 length(normMoveOdd)]);
set(stimLine, 'Color','k','LineStyle','--');
xlim([-.5 1.5]);
box off;
set(gca,'tickdir','out')
title('Outcome')
set(gca,'ytick',[])
set(gca,'xtick',[-1 0 1])

%% 5. Left vs right direction preference

% compute each cell's mean event-aligned activity (even trials vs odd trials)
mR_responses = neuralData.eta.alignedResps{2}(mR,:,:);
mL_responses = neuralData.eta.alignedResps{2}(mL,:,:);

mL_even = squeeze(nanmean(mL_responses(2:2:end,:,:),1))';
mL_odd = squeeze(nanmean(mL_responses(1:2:end,:,:),1))';
mR_even = squeeze(nanmean(mR_responses(2:2:end,:,:),1))';
mR_odd = squeeze(nanmean(mR_responses(1:2:end,:,:),1))';

R_evens(:,:,1) = mL_even;
R_evens(:,:,2) = mR_even;
R_odds(:,:,1) = mL_odd;
R_odds(:,:,2) = mR_odd;

%normalize the responses across the mL and mR responses (so you preserve
%relative activity between the two conditions even while normalizing to 1)
for c = 1:length(R_evens)
    cellSlice = squeeze(R_evens(c,:,:));
    sliceMin = min(min(cellSlice));
    normSlice = cellSlice - sliceMin;
    sliceMax = max(max(normSlice));
    R_evens_norm (c,:,:) = normSlice / sliceMax;
end

for c = 1:length(R_odds)
    cellSlice = squeeze(R_odds(c,:,:));
    sliceMin = min(min(cellSlice));
    normSlice = cellSlice - sliceMin;
    sliceMax = max(max(normSlice));
    R_odds_norm (c,:,:) = normSlice / sliceMax;
end

assessPeak = 17:30;
respR = max(mean(R_evens_norm(:,assessPeak,2),3),[],2);
respL = max(mean(R_evens_norm(:,assessPeak,1),3),[],2);
dirPref = (respR - respL);
[~,dirIdx] = sort(dirPref);
cmin = min(min(min(R_evens_norm))) + 0.025;
cmax = max(max(max(R_evens_norm))) - 0.025;

%% PLOT R - L responsiveness
cm = flipud(gray);

figure;
subplot(1,9,1);
line([1 1],[length(normMoveOdd)-500 length(normMoveOdd)],'LineWidth',3,'Color','k')
xlim([0 1])
ylim([1 length(normMoveOdd)])
h=text(.8,length(normMoveOdd)-500,'  500 cells');
set(h,'Rotation',90);
axis off

subplot(1,9,[2 3]);
imagesc(...
    neuralData.eta.eventWindow,...
    1:length(R_odds_norm),...
    smoothdata(R_odds_norm(dirIdx,:,1),2,'gaussian',5))
hold on;
stimLine = line([0 0],[1 length(R_odds_norm)]);
set(stimLine, 'Color','k','LineStyle','--');
xlim([-1 1]);
box off;
set(gca,'tickdir','out')
title('Movement (L)')
set(gca,'ytick',[])
set(gca,'xtick',[-1 0 1])
xlabel('Time (s)')

subplot(1,9,[5 6]);
imagesc(...
    neuralData.eta.eventWindow,...
    1:length(R_odds_norm),...
    smoothdata(R_odds_norm(dirIdx,:,2),2,'gaussian',5))
caxis([cmin cmax]);
stimLine = line([0 0],[1 length(R_odds_norm)]);
set(stimLine, 'Color','k','LineStyle','--');
colormap(cm)
xlim([-1 1]);
box off;
set(gca,'tickdir','out')
title('Movement (R)')
set(gca,'ytick',[])
set(gca,'xtick',[-1 0 1])
%%
figure;
imagesc(...
    neuralData.eta.eventWindow,...
    1:length(R_odds_norm),...
    smoothdata(...
        R_odds_norm(dirIdx,:,2)-R_odds_norm(dirIdx,:,1),...
        2,'gaussian',5));



