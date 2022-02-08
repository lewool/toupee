for m = 1:length(mouseList)
    
mouseName = char(mouseList{m});
expDate = char(expList{m}{1});
expNum = expList{m}{2};
hemisphere = hemList(m);
[expRef, ~] = data.constructExpRef(mouseName,expDate,expNum);
fprintf('fetching %s...\n',expRef)
[behavioralData, expInfo, neuralData] = data.loadDataset(mouseName, expDate, expNum);
[neuralData] = getSignificantActivity(expInfo, behavioralData, neuralData,0);
eventWindow = neuralData.eta.eventWindow;

%%
try
trialTypes = getTrialTypes(expInfo, behavioralData, 'early');
[baselineResps, stimResps, pmovResps, movResps, rewResps, preCueResps] = getEpochResps(neuralData.eta);
    
whichETA = 2;

sL_trials = trialTypes.singleVar.side{1};
sR_trials = trialTypes.singleVar.side{3};
mL_trials = trialTypes.singleVar.direction{1};
mR_trials = trialTypes.singleVar.direction{2};
hit_trials = trialTypes.singleVar.outcome{1};
miss_trials = trialTypes.singleVar.outcome{2};
mLhit_trials = trialTypes.intVar.cb3D.direction_block{1,1};
mRhit_trials = trialTypes.intVar.cb3D.direction_block{2,1};
mLmiss_trials = trialTypes.intVar.cb3D.direction_block{1,2};
mRmiss_trials = trialTypes.intVar.cb3D.direction_block{2,2};
sLhit_trials = trialTypes.intVar.cb3D.side_block{1,1};
sRhit_trials = trialTypes.intVar.cb3D.side_block{3,1};
sLmiss_trials = trialTypes.intVar.cb3D.side_block{1,2};
sRmiss_trials = trialTypes.intVar.cb3D.side_block{3,2};

%% movement
whichCells = find(neuralData.stats.bfcH(:,5) | neuralData.stats.bfcH(:,6));

leftMovCells = find(neuralData.stats.bfcH(:,5));
rightMovCells = find(neuralData.stats.bfcH(:,6));


mL_resps = squeeze(mean(neuralData.eta.alignedResps{whichETA}(mL_trials,:,whichCells)))' - mean(baselineResps(mL_trials,whichCells),1)';
mR_resps = squeeze(mean(neuralData.eta.alignedResps{whichETA}(mR_trials,:,whichCells)))' - mean(baselineResps(mR_trials,whichCells),1)';

%diff_resps = contraTrials - ipsiTrials (so contra responses are positive)
if hemisphere > 0
    diff_resps = mL_resps - mR_resps;
else
    diff_resps = mR_resps - mL_resps;
end
[maxx] = max(diff_resps(:,30:35),[],2);
[~,sortIdx] = sort(maxx,'descend');

figure;
set(gcf,'position',[694 666 1122 960]);
subplot(3,3,[2 5]);
imagesc(neuralData.eta.eventWindow, [1 size(diff_resps,1)],diff_resps(sortIdx,:))
colormap(colormapThruWhite([0 .75 0],[.75 0 .75],100,.6))
caxis([-1 1]);
hold on;
line([0 0],[1 numel(whichCells)],'LineStyle','--','Color',[.5 .5 .5]);
line([0.8 0.8],[1 numel(whichCells)],'LineStyle','--','Color',[.5 .5 .5]);
xlim([-.5 2]);
box off
axis off
title('\DeltaResponse (contra – ipsi choice)')

if hemisphere > 0
    rangeContraCells = [1 numel(leftMovCells)];
    rangeIpsiCells = [numel(leftMovCells)+1 numel([leftMovCells; rightMovCells])];
else
    rangeContraCells = [1 numel(rightMovCells)];
    rangeIpsiCells = [numel(rightMovCells)+1 numel([leftMovCells; rightMovCells])];
end

contraLabel = line([-.5 -.5],rangeContraCells,'LineWidth',7,'Color',[0 .4 1]);
ipsiLabel = line([-.5 -.5],rangeIpsiCells,'LineWidth',7,'Color','r');
h=text(-.8, mean(rangeContraCells),strcat({'contraMov cells (n = '},num2str(diff(rangeContraCells)+1),')'));
set(h,'Rotation',90,'HorizontalAlignment', 'center');
h=text(-.8, mean(rangeIpsiCells),strcat({'ipsiMov cells (n = '},num2str(diff(rangeIpsiCells)+1),')'));
set(h,'Rotation',90,'HorizontalAlignment', 'center');

resps_sLhitTrials_lCells = squeeze(mean(neuralData.eta.alignedResps{whichETA}(sLhit_trials,:,leftMovCells)))' - mean(baselineResps(sLhit_trials,leftMovCells),1)';
resps_sRhitTrials_lCells = squeeze(mean(neuralData.eta.alignedResps{whichETA}(sRhit_trials,:,leftMovCells)))' - mean(baselineResps(sRhit_trials,leftMovCells),1)';
resps_sLhitTrials_rCells = squeeze(mean(neuralData.eta.alignedResps{whichETA}(sLhit_trials,:,rightMovCells)))' - mean(baselineResps(sLhit_trials,rightMovCells),1)';
resps_sRhitTrials_rCells = squeeze(mean(neuralData.eta.alignedResps{whichETA}(sRhit_trials,:,rightMovCells)))' - mean(baselineResps(sRhit_trials,rightMovCells),1)';

resps_sLmissTrials_lCells = squeeze(mean(neuralData.eta.alignedResps{whichETA}(sLmiss_trials,:,leftMovCells)))' - mean(baselineResps(sLmiss_trials,leftMovCells),1)';
resps_sRmissTrials_lCells = squeeze(mean(neuralData.eta.alignedResps{whichETA}(sRmiss_trials,:,leftMovCells)))' - mean(baselineResps(sRmiss_trials,leftMovCells),1)';
resps_sLmissTrials_rCells = squeeze(mean(neuralData.eta.alignedResps{whichETA}(sLmiss_trials,:,rightMovCells)))' - mean(baselineResps(sLmiss_trials,rightMovCells),1)';
resps_sRmissTrials_rCells = squeeze(mean(neuralData.eta.alignedResps{whichETA}(sRmiss_trials,:,rightMovCells)))' - mean(baselineResps(sRmiss_trials,rightMovCells),1)';

if hemisphere > 0
    %deltaResp = contraTrials - ipsiTrials (so contra responses are positive)
    deltaHitResps_contraCells = mean(resps_sLhitTrials_lCells - resps_sRhitTrials_lCells,1);
    deltaHitResps_ipsiCells = mean(resps_sLhitTrials_rCells - resps_sRhitTrials_rCells,1);
    deltaMissResps_contraCells = mean(resps_sLmissTrials_lCells - resps_sRmissTrials_lCells,1);
    deltaMissResps_ipsiCells = mean(resps_sLmissTrials_rCells - resps_sRmissTrials_rCells,1);
    numContraMovCells(m) = numel(leftMovCells);
    numIpsiMovCells(m) = numel(rightMovCells);
else
    %deltaResp = contraTrials - ipsiTrials (so contra responses are positive)
    deltaHitResps_contraCells = mean(resps_sRhitTrials_rCells - resps_sLhitTrials_rCells,1);
    deltaHitResps_ipsiCells = mean(resps_sRhitTrials_lCells - resps_sLhitTrials_lCells,1);
    deltaMissResps_contraCells = mean(resps_sRmissTrials_rCells - resps_sLmissTrials_rCells,1);
    deltaMissResps_ipsiCells = mean(resps_sRmissTrials_lCells - resps_sLmissTrials_lCells,1);
    numContraMovCells(m) = numel(rightMovCells);
    numIpsiMovCells(m) = numel(leftMovCells);
end

% resps_mLhitTrials_lCells = squeeze(mean(neuralData.eta.alignedResps{whichETA}(mLhit_trials,:,leftMovCells)))' - mean(baselineResps(mLhit_trials,leftMovCells),1)';
% resps_mRhitTrials_lCells = squeeze(mean(neuralData.eta.alignedResps{whichETA}(mRhit_trials,:,leftMovCells)))' - mean(baselineResps(mRhit_trials,leftMovCells),1)';
% resps_mLhitTrials_rCells = squeeze(mean(neuralData.eta.alignedResps{whichETA}(mLhit_trials,:,rightMovCells)))' - mean(baselineResps(mLhit_trials,rightMovCells),1)';
% resps_mRhitTrials_rCells = squeeze(mean(neuralData.eta.alignedResps{whichETA}(mRhit_trials,:,rightMovCells)))' - mean(baselineResps(mRhit_trials,rightMovCells),1)';
% 
% resps_mLmissTrials_lCells = squeeze(mean(neuralData.eta.alignedResps{whichETA}(mLmiss_trials,:,leftMovCells)))' - mean(baselineResps(mLmiss_trials,leftMovCells),1)';
% resps_mRmissTrials_lCells = squeeze(mean(neuralData.eta.alignedResps{whichETA}(mRmiss_trials,:,leftMovCells)))' - mean(baselineResps(mRmiss_trials,leftMovCells),1)';
% resps_mLmissTrials_rCells = squeeze(mean(neuralData.eta.alignedResps{whichETA}(mLmiss_trials,:,rightMovCells)))' - mean(baselineResps(mLmiss_trials,rightMovCells),1)';
% resps_mRmissTrials_rCells = squeeze(mean(neuralData.eta.alignedResps{whichETA}(mRmiss_trials,:,rightMovCells)))' - mean(baselineResps(mRmiss_trials,rightMovCells),1)';
% 
% if hemisphere > 0
%     %deltaResp = contraTrials - ipsiTrials (so contra responses are positive)
%     deltaHitResps_contraCells = mean(resps_mLhitTrials_lCells - resps_mRhitTrials_lCells,1);
%     deltaHitResps_ipsiCells = mean(resps_mLhitTrials_rCells - resps_mRhitTrials_rCells,1);
%     deltaMissResps_contraCells = mean(resps_mLmissTrials_lCells - resps_mRmissTrials_lCells,1);
%     deltaMissResps_ipsiCells = mean(resps_mLmissTrials_rCells - resps_mRmissTrials_rCells,1);
% else
%     %deltaResp = contraTrials - ipsiTrials (so contra responses are positive)
%     deltaHitResps_contraCells = mean(resps_mRhitTrials_rCells - resps_mLhitTrials_rCells,1);
%     deltaHitResps_ipsiCells = mean(resps_mRhitTrials_lCells - resps_mLhitTrials_lCells,1);
%     deltaMissResps_contraCells = mean(resps_mRmissTrials_rCells - resps_mLmissTrials_rCells,1);
%     deltaMissResps_ipsiCells = mean(resps_mRmissTrials_lCells - resps_mLmissTrials_lCells,1);
% end

minY = 1.1*min([deltaHitResps_contraCells deltaHitResps_ipsiCells deltaMissResps_contraCells deltaMissResps_ipsiCells]);
maxY = 1.1*max([deltaHitResps_contraCells deltaHitResps_ipsiCells deltaMissResps_contraCells deltaMissResps_ipsiCells]);
absMax = max([abs(minY) abs(maxY)]);

subplot(3,3,[8]);
pHc = plot(neuralData.eta.eventWindow, smooth(deltaHitResps_contraCells,'lowess'));
hold on;
pHi = plot(neuralData.eta.eventWindow, smooth(deltaHitResps_ipsiCells,'lowess'));
pMc = plot(neuralData.eta.eventWindow, smooth(deltaMissResps_contraCells,'lowess'));
pMi = plot(neuralData.eta.eventWindow, smooth(deltaMissResps_ipsiCells,'lowess'));

set([pHc pHi pMc pMi],'LineWidth',2);
set([pHc pMc],'Color',[0 .4 1]);
set([pHi pMi],'Color',[1 0 0]);
set([pHc pHi],'LineStyle','-');
set([pMc pMi],'LineStyle',':');

zeroLine = line([-.5 2],[0 0],'LineStyle','-','Color',[.5 .5 .5]);
stimOnLine = line([0 0],[-absMax absMax],'LineStyle','--','Color',[.5 .5 .5]);
earliestMoveLine = line([0.8 0.8],[-absMax absMax],'LineStyle','--','Color',[.5 .5 .5]);
% line([-.5 0],[-absMax -absMax],'LineWidth',5,'Color',[0.75 .75 .75]);
% line([0 .8],[-absMax -absMax],'LineWidth',5,'Color',[0 .75 .75]);
% line([.8 2],[-absMax -absMax],'LineWidth',5,'Color',[0 .5 .5]);
h=text(-.25,-absMax+absMax*.1,'\it pre');
set(h,'Color',[.5 .5 .5],'HorizontalAlignment', 'center');
h=text(0.4,-absMax+absMax*.1,'\it delay');
set(h,'Color',[.5 .5 .5],'HorizontalAlignment', 'center');
h=text(1.4,-absMax+absMax*.1,'\it response');
set(h,'Color',[.5 .5 .5],'HorizontalAlignment', 'center');


ylim([-absMax absMax]);
xlim([-.5 2]);
box off
set(gca,'tickdir','out')
ylabel('\DeltaResponse (contra – ipsi trials)')
xlabel('Time from stimulus (s)')

hitResps_contraM(m,:) = deltaHitResps_contraCells;
missResps_contraM(m,:) = deltaMissResps_contraCells;
hitResps_ipsiM(m,:) = deltaHitResps_ipsiCells;
missResps_ipsiM(m,:) = deltaMissResps_ipsiCells;

%% outcome
whichCells = find(neuralData.stats.bfcH(:,7) | neuralData.stats.bfcH(:,9));
hitCells = find(neuralData.stats.bfcH(:,7));
missCells = find(neuralData.stats.bfcH(:,9));

hit_resps = squeeze(mean(neuralData.eta.alignedResps{whichETA}(hit_trials,:,whichCells)))' - mean(baselineResps(hit_trials,whichCells),1)';
miss_resps = squeeze(mean(neuralData.eta.alignedResps{whichETA}(miss_trials,:,whichCells)))' - mean(baselineResps(miss_trials,whichCells),1)';
diff_resps = hit_resps - miss_resps;
[maxx] = max(diff_resps(:,30:40),[],2);
[~,sortIdx] = sort(maxx,'descend');

subplot(3,3,[3 6]);
imagesc(neuralData.eta.eventWindow, [1 size(diff_resps,1)],diff_resps(sortIdx,:))
colormap(colormapThruWhite([0 .75 0],[.75 0 .75],100,.6))
caxis([-1 1]);
hold on;
line([0 0],[1 numel(whichCells)],'LineStyle','--','Color',[.5 .5 .5]);
line([0.8 0.8],[1 numel(whichCells)],'LineStyle','--','Color',[.5 .5 .5]);
xlim([-0.5 2]);
box off
axis off
title('\DeltaResponse (correct – error)')

rangeHitCells = [1 numel(hitCells)];
rangeMissCells = [numel(hitCells)+1 numel([hitCells; missCells])];

contraLabel = line([-.5 -.5],rangeHitCells,'LineStyle','-','LineWidth',7,'Color',[0 .75 .1]);
ipsiLabel = line([-.5 -.5],rangeMissCells,'LineStyle','-','LineWidth',7,'Color',[.75 0 0]);
h=text(-.8, mean(rangeHitCells),strcat({'Correct cells (n = '},num2str(diff(rangeHitCells)+1),')'));
set(h,'Rotation',90,'HorizontalAlignment', 'center');
h=text(-.8, mean(rangeMissCells),strcat({'Error cells (n = '},num2str(diff(rangeMissCells)+1),')'));
set(h,'Rotation',90,'HorizontalAlignment', 'center');

resps_sLhitTrials_hCells = squeeze(mean(neuralData.eta.alignedResps{whichETA}(sLhit_trials,:,hitCells)))' - mean(baselineResps(sLhit_trials,hitCells),1)';
resps_sRhitTrials_hCells = squeeze(mean(neuralData.eta.alignedResps{whichETA}(sRhit_trials,:,hitCells)))' - mean(baselineResps(sRhit_trials,hitCells),1)';
resps_sLhitTrials_mCells = squeeze(mean(neuralData.eta.alignedResps{whichETA}(sLhit_trials,:,missCells)))' - mean(baselineResps(sLhit_trials,missCells),1)';
resps_sRhitTrials_mCells = squeeze(mean(neuralData.eta.alignedResps{whichETA}(sRhit_trials,:,missCells)))' - mean(baselineResps(sRhit_trials,missCells),1)';

resps_sLmissTrials_hCells = squeeze(mean(neuralData.eta.alignedResps{whichETA}(sLmiss_trials,:,hitCells)))' - mean(baselineResps(sLmiss_trials,hitCells),1)';
resps_sRmissTrials_hCells = squeeze(mean(neuralData.eta.alignedResps{whichETA}(sRmiss_trials,:,hitCells)))' - mean(baselineResps(sRmiss_trials,hitCells),1)';
resps_sLmissTrials_mCells = squeeze(mean(neuralData.eta.alignedResps{whichETA}(sLmiss_trials,:,missCells)))' - mean(baselineResps(sLmiss_trials,missCells),1)';
resps_sRmissTrials_mCells = squeeze(mean(neuralData.eta.alignedResps{whichETA}(sRmiss_trials,:,missCells)))' - mean(baselineResps(sRmiss_trials,missCells),1)';

if hemisphere > 0
    %deltaResp = contraTrials - ipsiTrials (so contra responses are positive)
    deltaHitResps_hitCells = mean(resps_sLhitTrials_hCells - resps_sRhitTrials_hCells,1);
    deltaHitResps_missCells = mean(resps_sLhitTrials_mCells - resps_sRhitTrials_mCells,1);
    deltaMissResps_hitCells = mean(resps_sLmissTrials_hCells - resps_sRmissTrials_hCells,1);
    deltaMissResps_missCells = mean(resps_sLmissTrials_mCells - resps_sRmissTrials_mCells,1);
    
else
    %deltaResp = contraTrials - ipsiTrials (so contra responses are positive)
    deltaHitResps_hitCells = mean(resps_sRhitTrials_hCells - resps_sLhitTrials_hCells,1);
    deltaHitResps_missCells = mean(resps_sRhitTrials_mCells - resps_sLhitTrials_mCells,1);
    deltaMissResps_hitCells = mean(resps_sRmissTrials_hCells - resps_sLmissTrials_hCells,1);
    deltaMissResps_missCells = mean(resps_sRmissTrials_mCells - resps_sLmissTrials_mCells,1);
end

numHitCells(m) = numel(hitCells);
numMissCells(m) = numel(missCells);

% resps_mLhitTrials_hCells = squeeze(mean(neuralData.eta.alignedResps{whichETA}(mLhit_trials,:,hitCells)))' - mean(baselineResps(mLhit_trials,hitCells),1)';
% resps_mRhitTrials_hCells = squeeze(mean(neuralData.eta.alignedResps{whichETA}(mRhit_trials,:,hitCells)))' - mean(baselineResps(mRhit_trials,hitCells),1)';
% resps_mLhitTrials_mCells = squeeze(mean(neuralData.eta.alignedResps{whichETA}(mLhit_trials,:,missCells)))' - mean(baselineResps(mLhit_trials,missCells),1)';
% resps_mRhitTrials_mCells = squeeze(mean(neuralData.eta.alignedResps{whichETA}(mRhit_trials,:,missCells)))' - mean(baselineResps(mRhit_trials,missCells),1)';
% 
% resps_mLmissTrials_hCells = squeeze(mean(neuralData.eta.alignedResps{whichETA}(mLmiss_trials,:,hitCells)))' - mean(baselineResps(mLmiss_trials,hitCells),1)';
% resps_mRmissTrials_hCells = squeeze(mean(neuralData.eta.alignedResps{whichETA}(mRmiss_trials,:,hitCells)))' - mean(baselineResps(mRmiss_trials,hitCells),1)';
% resps_mLmissTrials_mCells = squeeze(mean(neuralData.eta.alignedResps{whichETA}(mLmiss_trials,:,missCells)))' - mean(baselineResps(mLmiss_trials,missCells),1)';
% resps_mRmissTrials_mCells = squeeze(mean(neuralData.eta.alignedResps{whichETA}(mRmiss_trials,:,missCells)))' - mean(baselineResps(mRmiss_trials,missCells),1)';
% 
% if hemisphere > 0
%     %deltaResp = contraTrials - ipsiTrials (so contra responses are positive)
%     deltaContraResps_hitCells = mean(resps_mLhitTrials_hCells - resps_mLmissTrials_hCells,1);
%     deltaIpsiResps_hitCells = mean(resps_mRhitTrials_hCells - resps_mRmissTrials_hCells,1);
%     deltaContraResps_missCells = mean(resps_mLhitTrials_mCells - resps_mLmissTrials_mCells,1);
%     deltaIpsiResps_missCells = mean(resps_mRhitTrials_mCells - resps_mRmissTrials_mCells,1);
% else
%     %deltaResp = contraTrials - ipsiTrials (so contra responses are positive)
%     deltaContraResps_hitCells = mean(resps_mRhitTrials_hCells - resps_mRmissTrials_hCells,1);
%     deltaIpsiResps_hitCells = mean(resps_mLhitTrials_hCells - resps_mLmissTrials_hCells,1);
%     deltaContraResps_missCells = mean(resps_mRhitTrials_mCells - resps_mRmissTrials_mCells,1);
%     deltaIpsiResps_missCells = mean(resps_mLhitTrials_mCells - resps_mLmissTrials_mCells,1);
% end

% minY = 1.1*min([deltaContraResps_hitCells deltaIpsiResps_hitCells deltaContraResps_missCells deltaIpsiResps_missCells]);
% maxY = 1.1*max([deltaContraResps_hitCells deltaIpsiResps_hitCells deltaContraResps_missCells deltaIpsiResps_missCells]);

minY = 1.1*min([deltaHitResps_hitCells deltaHitResps_missCells deltaMissResps_hitCells deltaMissResps_missCells]);
maxY = 1.1*max([deltaHitResps_hitCells deltaHitResps_missCells deltaMissResps_hitCells deltaMissResps_missCells]);
absMax = max([abs(minY) abs(maxY)]);

subplot(3,3,[9]);
pHh = plot(neuralData.eta.eventWindow, smooth(deltaHitResps_hitCells,'lowess'));
hold on;
pHm = plot(neuralData.eta.eventWindow, smooth(deltaHitResps_missCells,'lowess'));
pMh = plot(neuralData.eta.eventWindow, smooth(deltaMissResps_hitCells,'lowess'));
pMm = plot(neuralData.eta.eventWindow, smooth(deltaMissResps_missCells,'lowess'));

set([pHh pHm pMh pMm],'LineWidth',2);
set([pHh pMh],'Color',[0 .5 .1]);
set([pHm pMm],'Color',[.75 0 0]);
set([pHh pHm],'LineStyle','-');
set([pMh pMm],'LineStyle',':');

zeroLine = line([-.5 2],[0 0],'LineStyle','-','Color',[.5 .5 .5]);
stimOnLine = line([0 0],[-absMax absMax],'LineStyle','--','Color',[.5 .5 .5]);
earliestMoveLine = line([0.8 0.8],[-absMax absMax],'LineStyle','--','Color',[.5 .5 .5]);
% line([-.5 0],[-absMax -absMax],'LineWidth',5,'Color',[0.75 .75 .75]);
% line([0 .8],[-absMax -absMax],'LineWidth',5,'Color',[0 .75 .75]);
% line([.8 2],[-absMax -absMax],'LineWidth',5,'Color',[0 .5 .5]);
h=text(-.25,-absMax+absMax*.1,'\it pre');
set(h,'Color',[.5 .5 .5],'HorizontalAlignment', 'center');
h=text(0.4,-absMax+absMax*.1,'\it delay');
set(h,'Color',[.5 .5 .5],'HorizontalAlignment', 'center');
h=text(1.4,-absMax+absMax*.1,'\it response');
set(h,'Color',[.5 .5 .5],'HorizontalAlignment', 'center');


ylim([-absMax absMax]);
xlim([-.5 2]);
box off
set(gca,'tickdir','out')
ylabel('\DeltaResponse (contra – ipsi trials)')
xlabel('Time from stimulus (s)')

% contraResps_hit(m,:) = deltaContraResps_hitCells;
% ipsiResps_hit(m,:) = deltaIpsiResps_hitCells;
% contraResps_miss(m,:) = deltaContraResps_missCells;
% ipsiResps_miss(m,:) = deltaIpsiResps_missCells;

hitResps_hit(m,:) = deltaHitResps_hitCells;
hitResps_miss(m,:) = deltaHitResps_missCells;
missResps_hit(m,:) = deltaMissResps_hitCells;
missResps_miss(m,:) = deltaMissResps_missCells;

%% stim
whichCells = find(neuralData.stats.bfcH(:,2) | neuralData.stats.bfcH(:,3));
leftStimCells = find(neuralData.stats.bfcH(:,2));
rightStimCells = find(neuralData.stats.bfcH(:,3));


sL_resps = squeeze(mean(neuralData.eta.alignedResps{whichETA}(sL_trials,:,whichCells)))' - mean(baselineResps(sL_trials,whichCells),1)';
sR_resps = squeeze(mean(neuralData.eta.alignedResps{whichETA}(sR_trials,:,whichCells)))' - mean(baselineResps(sR_trials,whichCells),1)';
if hemisphere > 0
    diff_resps = sL_resps - sR_resps;
else
    diff_resps = sR_resps - sL_resps;
end
[maxx] = max(diff_resps(:,20:25),[],2);
[~,sortIdx] = sort(maxx,'descend');

subplot(3,3,[1 4]);
imagesc(neuralData.eta.eventWindow, [1 size(diff_resps,1)],diff_resps(sortIdx,:))
colormap(colormapThruWhite([0 .75 0],[.75 0 .75],100,.6))
caxis([-2 2]);
hold on;
line([0 0],[1 numel(whichCells)],'LineStyle','--','Color',[.5 .5 .5]);
line([0.8 0.8],[1 numel(whichCells)],'LineStyle','--','Color',[.5 .5 .5]);
xlim([-0.5 2]);
box off
axis off
title('\DeltaResponse (contra – ipsi stimulus)')

if hemisphere > 0
    rangeContraCells = [1 numel(leftStimCells)];
    rangeIpsiCells = [numel(leftStimCells)+1 numel([leftStimCells; rightStimCells])];
else
    rangeContraCells = [1 numel(rightStimCells)];
    rangeIpsiCells = [numel(rightStimCells)+1 numel([leftStimCells; rightStimCells])];
end

contraLabel = line([-.5 -.5],rangeContraCells,'LineWidth',7,'Color',[0 .4 1]);
ipsiLabel = line([-.5 -.5],rangeIpsiCells,'LineWidth',7,'Color','r');
h=text(-.8, mean(rangeContraCells),strcat({'contraStim cells (n = '},num2str(diff(rangeContraCells)+1),')'));
set(h,'Rotation',90,'HorizontalAlignment', 'center');
h=text(-.8, mean(rangeIpsiCells),strcat({'ipsiStim cells (n = '},num2str(diff(rangeIpsiCells)+1),')'));
set(h,'Rotation',90,'HorizontalAlignment', 'center');

try
    resps_sLhitTrials_lCells = squeeze(mean(neuralData.eta.alignedResps{whichETA}(sLhit_trials,:,leftStimCells),1))' - mean(baselineResps(mLhit_trials,leftStimCells),1)';
catch
    resps_sLhitTrials_lCells = nan(1,size(neuralData.eta.alignedResps{whichETA},2));
end

try
    resps_sRhitTrials_lCells = squeeze(mean(neuralData.eta.alignedResps{whichETA}(sRhit_trials,:,leftStimCells),1))' - mean(baselineResps(mRhit_trials,leftStimCells),1)';
catch
    resps_sRhitTrials_lCells = nan(1,size(neuralData.eta.alignedResps{whichETA},2));
end

try
    resps_sLhitTrials_rCells = squeeze(mean(neuralData.eta.alignedResps{whichETA}(sLhit_trials,:,rightStimCells),1))' - mean(baselineResps(mLhit_trials,rightStimCells),1)';
catch
    resps_sLhitTrials_rCells = nan(1,size(neuralData.eta.alignedResps{whichETA},2));
end

try
    resps_sRhitTrials_rCells = squeeze(mean(neuralData.eta.alignedResps{whichETA}(sRhit_trials,:,rightStimCells),1))' - mean(baselineResps(mRhit_trials,rightStimCells),1)';
catch
    resps_sRhitTrials_rCells = nan(1,size(neuralData.eta.alignedResps{whichETA},2));
end

try
    resps_sLmissTrials_lCells = squeeze(mean(neuralData.eta.alignedResps{whichETA}(sLmiss_trials,:,leftStimCells),1))' - mean(baselineResps(mLmiss_trials,leftStimCells),1)';
catch
    resps_sLmissTrials_lCells = nan(1,size(neuralData.eta.alignedResps{whichETA},2));
end

try
    resps_sRmissTrials_lCells = squeeze(mean(neuralData.eta.alignedResps{whichETA}(sRmiss_trials,:,leftStimCells),1))' - mean(baselineResps(mRmiss_trials,leftStimCells),1)';
catch
    resps_sRmissTrials_lCells = nan(1,size(neuralData.eta.alignedResps{whichETA},2));
end

try
    resps_sLmissTrials_rCells = squeeze(mean(neuralData.eta.alignedResps{whichETA}(sLmiss_trials,:,rightStimCells),1))' - mean(baselineResps(mLmiss_trials,rightStimCells),1)';
catch
    resps_sLmissTrials_rCells = nan(1,size(neuralData.eta.alignedResps{whichETA},2));
end

try
    resps_sRmissTrials_rCells = squeeze(mean(neuralData.eta.alignedResps{whichETA}(sRmiss_trials,:,rightStimCells),1))' - mean(baselineResps(mRmiss_trials,rightStimCells),1)';
catch
    resps_sRmissTrials_rCells = nan(1,size(neuralData.eta.alignedResps{whichETA},2));
end

if size(resps_sLhitTrials_lCells,2) == 1
    resps_sLhitTrials_lCells = resps_sLhitTrials_lCells';
end
if size(resps_sRhitTrials_lCells,2) == 1
    resps_sRhitTrials_lCells = resps_sRhitTrials_lCells';
end
if size(resps_sLhitTrials_rCells,2) == 1
    resps_sLhitTrials_rCells = resps_sLhitTrials_rCells';
end
if size(resps_sRhitTrials_rCells,2) == 1
    resps_sRhitTrials_rCells = resps_sRhitTrials_rCells';
end
if size(resps_sLmissTrials_lCells,2) == 1
    resps_sLmissTrials_lCells = resps_sLmissTrials_lCells';
end
if size(resps_sRmissTrials_lCells,2) == 1
    resps_sRmissTrials_lCells = resps_sRmissTrials_lCells';
end
if size(resps_sLmissTrials_rCells,2) == 1
    resps_sLmissTrials_rCells = resps_sLmissTrials_rCells';
end
if size(resps_sRmissTrials_rCells,2) == 1
    resps_sRmissTrials_rCells = resps_sRmissTrials_rCells';
end


if hemisphere > 0
    %deltaResp = contraTrials - ipsiTrials (so contra responses are positive)
    deltaHitResps_contraCells = mean(resps_sLhitTrials_lCells - resps_sRhitTrials_lCells,1);
    deltaHitResps_ipsiCells = mean(resps_sLhitTrials_rCells - resps_sRhitTrials_rCells,1);
    deltaMissResps_contraCells = mean(resps_sLmissTrials_lCells - resps_sRmissTrials_lCells,1);
    deltaMissResps_ipsiCells = mean(resps_sLmissTrials_rCells - resps_sRmissTrials_rCells,1);
    numContraStimCells(m) = numel(leftStimCells);
    numIpsiStimCells(m) = numel(rightStimCells);
else
    %deltaResp = contraTrials - ipsiTrials (so contra responses are positive)
    deltaHitResps_contraCells = mean(resps_sRhitTrials_rCells - resps_sLhitTrials_rCells,1);
    deltaHitResps_ipsiCells = mean(resps_sRhitTrials_lCells - resps_sLhitTrials_lCells,1);
    deltaMissResps_contraCells = mean(resps_sRmissTrials_rCells - resps_sLmissTrials_rCells,1);
    deltaMissResps_ipsiCells = mean(resps_sRmissTrials_lCells - resps_sLmissTrials_lCells,1);
    numContraStimCells(m) = numel(rightStimCells);
    numIpsiStimCells(m) = numel(leftStimCells);
end

minY = 1.1*min([deltaHitResps_contraCells deltaHitResps_ipsiCells deltaMissResps_contraCells deltaMissResps_ipsiCells]);
maxY = 1.1*max([deltaHitResps_contraCells deltaHitResps_ipsiCells deltaMissResps_contraCells deltaMissResps_ipsiCells]);
absMax = max([abs(minY) abs(maxY)]);

subplot(3,3,[7]);
pHcs = plot(neuralData.eta.eventWindow, smooth(deltaHitResps_contraCells,'lowess'));
hold on;
pHis = plot(neuralData.eta.eventWindow, smooth(deltaHitResps_ipsiCells,'lowess'));
pMcs = plot(neuralData.eta.eventWindow, smooth(deltaMissResps_contraCells,'lowess'));
pMis = plot(neuralData.eta.eventWindow, smooth(deltaMissResps_ipsiCells,'lowess'));

set([pHcs pHis pMcs pMis],'LineWidth',2);
set([pHcs pMcs],'Color',[0 .4 1]);
set([pHis pMis],'Color',[1 0 0]);
set([pHcs pHis],'LineStyle','-');
set([pMcs pMis],'LineStyle',':');

zeroLine = line([-.5 2],[0 0],'LineStyle','-','Color',[.5 .5 .5]);
stimOnLine = line([0 0],[-absMax absMax],'LineStyle','--','Color',[.5 .5 .5]);
earliestMoveLine = line([0.8 0.8],[-absMax absMax],'LineStyle','--','Color',[.5 .5 .5]);
% line([-.5 0],[-absMax -absMax],'LineWidth',5,'Color',[0.75 .75 .75]);
% line([0 .8],[-absMax -absMax],'LineWidth',5,'Color',[0 .75 .75]);
% line([.8 2],[-absMax -absMax],'LineWidth',5,'Color',[0 .5 .5]);
h=text(-.25,-absMax+absMax*.1,'\it pre');
set(h,'Color',[.5 .5 .5],'HorizontalAlignment', 'center');
h=text(0.4,-absMax+absMax*.1,'\it delay');
set(h,'Color',[.5 .5 .5],'HorizontalAlignment', 'center');
h=text(1.4,-absMax+absMax*.1,'\it response');
set(h,'Color',[.5 .5 .5],'HorizontalAlignment', 'center');


ylim([-absMax absMax]);
xlim([-.5 2]);
box off
set(gca,'tickdir','out')
ylabel('\DeltaResponse (contra – ipsi trials)')
xlabel('Time from stimulus (s)')

hitResps_contraS(m,:) = deltaHitResps_contraCells;
missResps_contraS(m,:) = deltaMissResps_contraCells;
hitResps_ipsiS(m,:) = deltaHitResps_ipsiCells;
missResps_ipsiS(m,:) = deltaMissResps_ipsiCells;

%%
catch
    hitResps_contraM(m,:) = nan(1,size(neuralData.eta.alignedResps{whichETA},2));
missResps_contraM(m,:) = nan(1,size(neuralData.eta.alignedResps{whichETA},2));
hitResps_ipsiM(m,:) = nan(1,size(neuralData.eta.alignedResps{whichETA},2));
missResps_ipsiM(m,:) = nan(1,size(neuralData.eta.alignedResps{whichETA},2));

    hitResps_hit(m,:) = nan(1,size(neuralData.eta.alignedResps{whichETA},2));
hitResps_miss(m,:) = nan(1,size(neuralData.eta.alignedResps{whichETA},2));
missResps_hit(m,:) = nan(1,size(neuralData.eta.alignedResps{whichETA},2));
missResps_miss(m,:) = nan(1,size(neuralData.eta.alignedResps{whichETA},2));

    hitResps_contraS(m,:) = nan(1,size(neuralData.eta.alignedResps{whichETA},2));
missResps_contraS(m,:) = nan(1,size(neuralData.eta.alignedResps{whichETA},2));
hitResps_ipsiS(m,:) = nan(1,size(neuralData.eta.alignedResps{whichETA},2));
missResps_ipsiS(m,:) = nan(1,size(neuralData.eta.alignedResps{whichETA},2));


end
clearvars -except mouseList expList hemList ...
    hitResps_contraM hitResps_ipsiM missResps_contraM missResps_ipsiM ...
    hitResps_hit hitResps_miss missResps_hit missResps_miss ...
    hitResps_contraS missResps_contraS hitResps_ipsiS missResps_ipsiS ...
    numContraMovCells numIpsiMovCells numHitCells numMissCells numContraStimCells numIpsiStimCells
%     contraResps_hit ipsiResps_hit contraResps_miss ipsiResps_miss ...

end

%%
eventWindow = linspace(-2,2,41);
absMax = 0.5;

figure;
subplot(1,3,1);
hold on;
plotPSTHs(eventWindow, nanmean(hitResps_contraS), nanstd(hitResps_contraS)/sqrt(size(missResps_ipsiS,1)),[0 0.4 1],'-')
plotPSTHs(eventWindow, nanmean(hitResps_ipsiS), nanstd(hitResps_ipsiS)/sqrt(size(missResps_ipsiS,1)),[1 0 0],'-')
plotPSTHs(eventWindow, nanmean(missResps_contraS), nanstd(missResps_contraS)/sqrt(size(missResps_ipsiS,1)),[0 0.4 1],':')
plotPSTHs(eventWindow, nanmean(missResps_ipsiS), nanstd(missResps_ipsiS)/sqrt(size(missResps_ipsiS,1)),[1 0 0],':')
xlim([-0.5 2]);
xlabel('Time from stimulus onset (s)')
ylabel('\DeltaResponse (contra – ipsi trials)')
zeroLine = line([-.5 2],[0 0],'LineStyle','-','Color',[.5 .5 .5]);
stimOnLine = line([0 0],[-absMax absMax],'LineStyle','--','Color',[.5 .5 .5]);
earliestMoveLine = line([0.8 0.8],[-absMax absMax],'LineStyle','--','Color',[.5 .5 .5]);
box off 
set(gca,'tickdir','out')

subplot(1,3,2);
hold on;
plotPSTHs(eventWindow, nanmean(hitResps_contraM), nanstd(hitResps_contraM)/sqrt(size(missResps_ipsiS,1)),[0 0.4 1],'-')
plotPSTHs(eventWindow, nanmean(hitResps_ipsiM), nanstd(hitResps_ipsiM)/sqrt(size(missResps_ipsiS,1)),[1 0 0],'-')
plotPSTHs(eventWindow, nanmean(missResps_contraM), nanstd(missResps_contraM)/sqrt(size(missResps_ipsiS,1)),[0 0.4 1],':')
plotPSTHs(eventWindow, nanmean(missResps_ipsiM), nanstd(missResps_ipsiM)/sqrt(size(missResps_ipsiS,1)),[1 0 0],':')
xlim([-0.5 2]);
xlabel('Time from stimulus onset (s)')
ylabel('\DeltaResponse (contra – ipsi trials)')
zeroLine = line([-.5 2],[0 0],'LineStyle','-','Color',[.5 .5 .5]);
stimOnLine = line([0 0],[-absMax absMax],'LineStyle','--','Color',[.5 .5 .5]);
earliestMoveLine = line([0.8 0.8],[-absMax absMax],'LineStyle','--','Color',[.5 .5 .5]);
box off 
set(gca,'tickdir','out')

subplot(1,3,3);
hold on;
plotPSTHs(eventWindow, nanmean(hitResps_hit), nanstd(hitResps_hit)/sqrt(size(missResps_ipsiS,1)),[0 .5 .1],'-')
plotPSTHs(eventWindow, nanmean(hitResps_miss), nanstd(hitResps_miss)/sqrt(size(missResps_ipsiS,1)),[.75 0 0],'-')
plotPSTHs(eventWindow, nanmean(missResps_hit), nanstd(missResps_hit)/sqrt(size(missResps_ipsiS,1)),[0 .5 .1],':')
plotPSTHs(eventWindow, nanmean(missResps_miss), nanstd(missResps_miss)/sqrt(size(missResps_ipsiS,1)),[.75 0 0],':')
xlim([-0.5 2]);
xlabel('Time from stimulus onset (s)')
ylabel('\DeltaResponse (contra – ipsi trials)')
zeroLine = line([-.5 2],[0 0],'LineStyle','-','Color',[.5 .5 .5]);
stimOnLine = line([0 0],[-absMax absMax],'LineStyle','--','Color',[.5 .5 .5]);
earliestMoveLine = line([0.8 0.8],[-absMax absMax],'LineStyle','--','Color',[.5 .5 .5]);
box off 
set(gca,'tickdir','out')

%% move
height = 0.39;
y = 0.5;
figure;
subplot(1,3,1)
xlim([.5 4.5])
ylim([-y y]); 
t1 = 15;
t2 = 20;
jitter = ((rand(1,24)-.5)/15);
% jitter = 0;
xPositions1 = jitter+ones(1,24);
xPositions2 = jitter+ones(1,24)*2.25;
xPositions3 = jitter+ones(1,24)*2.75;
xPositions4 = jitter+ones(1,24)*4;
yValues1 = mean(hitResps_contraM(:,t1:t2),2)';
yValues2 = mean(missResps_contraM(:,t1:t2),2)';
yValues3 = mean(hitResps_ipsiM(:,t1:t2),2)';
yValues4 = mean(missResps_ipsiM(:,t1:t2),2)';
lc = line([xPositions1; xPositions2],[yValues1; yValues2],'Color',[.5 .5 .5]);
li = line([xPositions3; xPositions4],[yValues3; yValues4],'Color',[.5 .5 .5]);
hold on;
plot(xPositions1,yValues1,'o','MarkerFaceColor',[0 .4 1],'MarkerEdgeColor','w')
plot(xPositions2,yValues2,'o','MarkerFaceColor',[0 .4 1],'MarkerEdgeColor','w')
plot(xPositions3,yValues3,'o','MarkerFaceColor',[1 0 0],'MarkerEdgeColor','w')
plot(xPositions4,yValues4,'o','MarkerFaceColor',[1 0 0],'MarkerEdgeColor','w')

if signrank(yValues1,yValues2) <= 0.0001
    h = text(1.7,height+.01,'****');
    set(h,'Color','k','HorizontalAlignment', 'center','FontSize',8);
    line([1 2.25],[height height],'Color','k')
elseif signrank(yValues1,yValues2) <= 0.001
    h = text(1.7,height+.01,'***');
    set(h,'Color','k','HorizontalAlignment', 'center','FontSize',8);
    line([1 2.25],[height height],'Color','k')
elseif signrank(yValues1,yValues2) <= 0.01
    h = text(1.7,height+.01,'**');
    set(h,'Color','k','HorizontalAlignment', 'center','FontSize',8);
    line([1 2.25],[height height],'Color','k')
elseif signrank(yValues1,yValues2) <= 0.05
    h = text(1.7,height+.01,'*');
    set(h,'Color','k','HorizontalAlignment', 'center','FontSize',8);
    line([1 2.25],[height height],'Color','k')
elseif signrank(yValues1,yValues2) > 0.05
    h = text(1.7,height+.01,'ns');
    set(h,'Color','k','VerticalAlignment','bottom','HorizontalAlignment', 'center','FontSize',8);
    line([1 2.25],[height height],'Color','k')
end

if signrank(yValues3,yValues4) <= 0.0001
    h = text(3.375,height+.01,'****');
    set(h,'Color','k','HorizontalAlignment', 'center','FontSize',8);
    line([2.75 4],[height height],'Color','k')
elseif signrank(yValues3, yValues4) <= 0.001
    h = text(3.375,height+.01,'****');
    set(h,'Color','k','HorizontalAlignment', 'center','FontSize',8);
    line([2.75 4],[height height],'Color','k')
elseif signrank(yValues3,yValues4) <= 0.01
    h = text(3.375,height+.01,'**');
    set(h,'Color','k','HorizontalAlignment', 'center','FontSize',8);
    line([2.75 4],[height height],'Color','k')
elseif signrank(yValues3,yValues4) <= 0.05
    h = text(3.375,height+.01,'*');
    set(h,'Color','k','HorizontalAlignment', 'center','FontSize',8);
    line([2.75 4],[height height],'Color','k')
elseif signrank(yValues3,yValues4) > 0.05
    h = text(3.375,height+.01,'ns');
    set(h,'Color','k','VerticalAlignment','bottom','HorizontalAlignment', 'center','FontSize',8);
    line([2.75 4],[height height],'Color','k')
end

for g = 1:24
    lc(g).Color(4) = .25;
    li(g).Color(4) = .25;
end
set(gca,'tickdir','out')
title('Pre-trial epoch (-500–0 ms)')
set(gca, 'XTick', [1 2.25 2.75 4])
set(gca, 'XTickLabels', {'corr.','err.','corr.','err.'})
xlabel('Contra cells                Ipsi cells')
ylabel('\DeltaResponse (contra – ipsi trials)')

subplot(1,3,2)
xlim([.5 4.5])
ylim([-y y]); 
t1 = 21;
t2 = 26;
jitter = ((rand(1,24)-.5)/15);
% jitter = 0;
xPositions1 = jitter+ones(1,24);
xPositions2 = jitter+ones(1,24)*2.25;
xPositions3 = jitter+ones(1,24)*2.75;
xPositions4 = jitter+ones(1,24)*4;
yValues1 = mean(hitResps_contraM(:,t1:t2),2)';
yValues2 = mean(missResps_contraM(:,t1:t2),2)';
yValues3 = mean(hitResps_ipsiM(:,t1:t2),2)';
yValues4 = mean(missResps_ipsiM(:,t1:t2),2)';
lc = line([xPositions1; xPositions2],[yValues1; yValues2],'Color',[.5 .5 .5]);
li = line([xPositions3; xPositions4],[yValues3; yValues4],'Color',[.5 .5 .5]);
hold on;
plot(xPositions1,yValues1,'o','MarkerFaceColor',[0 .4 1],'MarkerEdgeColor','w')
plot(xPositions2,yValues2,'o','MarkerFaceColor',[0 .4 1],'MarkerEdgeColor','w')
plot(xPositions3,yValues3,'o','MarkerFaceColor',[1 0 0],'MarkerEdgeColor','w')
plot(xPositions4,yValues4,'o','MarkerFaceColor',[1 0 0],'MarkerEdgeColor','w')

if signrank(yValues1,yValues2) <= 0.0001
    h = text(1.7,height+.01,'****');
    set(h,'Color','k','HorizontalAlignment', 'center','FontSize',8);
    line([1 2.25],[height height],'Color','k')
elseif signrank(yValues1,yValues2) <= 0.001
    h = text(1.7,height+.01,'***');
    set(h,'Color','k','HorizontalAlignment', 'center','FontSize',8);
    line([1 2.25],[height height],'Color','k')
elseif signrank(yValues1,yValues2) <= 0.01
    h = text(1.7,height+.01,'**');
    set(h,'Color','k','HorizontalAlignment', 'center','FontSize',8);
    line([1 2.25],[height height],'Color','k')
elseif signrank(yValues1,yValues2) <= 0.05
    h = text(1.7,height+.01,'*');
    set(h,'Color','k','HorizontalAlignment', 'center','FontSize',8);
    line([1 2.25],[height height],'Color','k')
elseif signrank(yValues1,yValues2) > 0.05
    h = text(1.7,height+.01,'ns');
    set(h,'Color','k','VerticalAlignment','bottom','HorizontalAlignment', 'center','FontSize',8);
    line([1 2.25],[height height],'Color','k')
end

if signrank(yValues3,yValues4) <= 0.0001
    h = text(3.375,height+.01,'****');
    set(h,'Color','k','HorizontalAlignment', 'center','FontSize',8);
    line([2.75 4],[height height],'Color','k')
elseif signrank(yValues3, yValues4) <= 0.001
    h = text(3.375,height+.01,'****');
    set(h,'Color','k','HorizontalAlignment', 'center','FontSize',8);
    line([2.75 4],[height height],'Color','k')
elseif signrank(yValues3,yValues4) <= 0.01
    h = text(3.375,height+.01,'**');
    set(h,'Color','k','HorizontalAlignment', 'center','FontSize',8);
    line([2.75 4],[height height],'Color','k')
elseif signrank(yValues3,yValues4) <= 0.05
    h = text(3.375,height+.01,'*');
    set(h,'Color','k','HorizontalAlignment', 'center','FontSize',8);
    line([2.75 4],[height height],'Color','k')
elseif signrank(yValues3,yValues4) > 0.05
    h = text(3.375,height+.01,'ns');
    set(h,'Color','k','VerticalAlignment','bottom','HorizontalAlignment', 'center','FontSize',8);
    line([2.75 4],[height height],'Color','k')
end


for g = 1:24
    lc(g).Color(4) = .25;
    li(g).Color(4) = .25;
end
set(gca,'tickdir','out')
title('Delay epoch (0–500 ms)')
set(gca, 'XTick', [1 2.25 2.75 4])
set(gca, 'XTickLabels', {'corr.','err.','corr.','err.'})
xlabel('Contra cells                Ipsi cells')
ylabel('\DeltaResponse (contra – ipsi trials)')

subplot(1,3,3)
xlim([.5 4.5])
ylim([-y y]); 
t1 = 28;
t2 = 35;
jitter = ((rand(1,24)-.5)/15);
% jitter = 0;
xPositions1 = jitter+ones(1,24);
xPositions2 = jitter+ones(1,24)*2.25;
xPositions3 = jitter+ones(1,24)*2.75;
xPositions4 = jitter+ones(1,24)*4;
yValues1 = mean(hitResps_contraM(:,t1:t2),2)';
yValues2 = mean(missResps_contraM(:,t1:t2),2)';
yValues3 = mean(hitResps_ipsiM(:,t1:t2),2)';
yValues4 = mean(missResps_ipsiM(:,t1:t2),2)';
lc = line([xPositions1; xPositions2],[yValues1; yValues2],'Color',[.5 .5 .5]);
li = line([xPositions3; xPositions4],[yValues3; yValues4],'Color',[.5 .5 .5]);
hold on;
plot(xPositions1,yValues1,'o','MarkerFaceColor',[0 .4 1],'MarkerEdgeColor','w')
plot(xPositions2,yValues2,'o','MarkerFaceColor',[0 .4 1],'MarkerEdgeColor','w')
plot(xPositions3,yValues3,'o','MarkerFaceColor',[1 0 0],'MarkerEdgeColor','w')
plot(xPositions4,yValues4,'o','MarkerFaceColor',[1 0 0],'MarkerEdgeColor','w')

if signrank(yValues1,yValues2) <= 0.0001
    h = text(1.7,height+.01,'****');
    set(h,'Color','k','HorizontalAlignment', 'center','FontSize',8);
    line([1 2.25],[height height],'Color','k')
elseif signrank(yValues1,yValues2) <= 0.001
    h = text(1.7,height+.01,'***');
    set(h,'Color','k','HorizontalAlignment', 'center','FontSize',8);
    line([1 2.25],[height height],'Color','k')
elseif signrank(yValues1,yValues2) <= 0.01
    h = text(1.7,height+.01,'**');
    set(h,'Color','k','HorizontalAlignment', 'center','FontSize',8);
    line([1 2.25],[height height],'Color','k')
elseif signrank(yValues1,yValues2) <= 0.05
    h = text(1.7,height+.01,'*');
    set(h,'Color','k','HorizontalAlignment', 'center','FontSize',8);
    line([1 2.25],[height height],'Color','k')
elseif signrank(yValues1,yValues2) > 0.05
    h = text(1.7,height+.01,'ns');
    set(h,'Color','k','VerticalAlignment','bottom','HorizontalAlignment', 'center','FontSize',8);
    line([1 2.25],[height height],'Color','k')
end

if signrank(yValues3,yValues4) <= 0.0001
    h = text(3.375,height+.01,'****');
    set(h,'Color','k','HorizontalAlignment', 'center','FontSize',8);
    line([2.75 4],[height height],'Color','k')
elseif signrank(yValues3, yValues4) <= 0.001
    h = text(3.375,height+.01,'****');
    set(h,'Color','k','HorizontalAlignment', 'center','FontSize',8);
    line([2.75 4],[height height],'Color','k')
elseif signrank(yValues3,yValues4) <= 0.01
    h = text(3.375,height+.01,'**');
    set(h,'Color','k','HorizontalAlignment', 'center','FontSize',8);
    line([2.75 4],[height height],'Color','k')
elseif signrank(yValues3,yValues4) <= 0.05
    h = text(3.375,height+.01,'*');
    set(h,'Color','k','HorizontalAlignment', 'center','FontSize',8);
    line([2.75 4],[height height],'Color','k')
elseif signrank(yValues3,yValues4) > 0.05
    h = text(3.375,height+.01,'ns');
    set(h,'Color','k','VerticalAlignment','bottom','HorizontalAlignment', 'center','FontSize',8);
    line([2.75 4],[height height],'Color','k')
end
for g = 1:24
    lc(g).Color(4) = .25;
    li(g).Color(4) = .25;
end
set(gca,'tickdir','out')
title('Response epoch (800–1400 ms)')
set(gca, 'XTick', [1 2.25 2.75 4])
set(gca, 'XTickLabels', {'corr.','err.','corr.','err.'})
xlabel('Contra cells                Ipsi cells')
ylabel('\DeltaResponse (contra – ipsi trials)')




%% stim

figure;
subplot(1,3,1)
xlim([.5 4.5])
ylim([-y y]); 
t1 = 15;
t2 = 20;
jitter = ((rand(1,24)-.5)/15);
% jitter = 0;
xPositions1 = jitter+ones(1,24);
xPositions2 = jitter+ones(1,24)*2.25;
xPositions3 = jitter+ones(1,24)*2.75;
xPositions4 = jitter+ones(1,24)*4;
yValues1 = mean(hitResps_contraS(:,t1:t2),2)';
yValues2 = mean(missResps_contraS(:,t1:t2),2)';
yValues3 = mean(hitResps_ipsiS(:,t1:t2),2)';
yValues4 = mean(missResps_ipsiS(:,t1:t2),2)';
lc = line([xPositions1; xPositions2],[yValues1; yValues2],'Color',[.5 .5 .5]);
li = line([xPositions3; xPositions4],[yValues3; yValues4],'Color',[.5 .5 .5]);
hold on;
plot(xPositions1,yValues1,'o','MarkerFaceColor',[0 .4 1],'MarkerEdgeColor','w')
plot(xPositions2,yValues2,'o','MarkerFaceColor',[0 .4 1],'MarkerEdgeColor','w')
plot(xPositions3,yValues3,'o','MarkerFaceColor',[1 0 0],'MarkerEdgeColor','w')
plot(xPositions4,yValues4,'o','MarkerFaceColor',[1 0 0],'MarkerEdgeColor','w')

if signrank(yValues1,yValues2) <= 0.0001
    h = text(1.7,height+.01,'****');
    set(h,'Color','k','HorizontalAlignment', 'center','FontSize',8);
    line([1 2.25],[height height],'Color','k')
elseif signrank(yValues1,yValues2) <= 0.001
    h = text(1.7,height+.01,'***');
    set(h,'Color','k','HorizontalAlignment', 'center','FontSize',8);
    line([1 2.25],[height height],'Color','k')
elseif signrank(yValues1,yValues2) <= 0.01
    h = text(1.7,height+.01,'**');
    set(h,'Color','k','HorizontalAlignment', 'center','FontSize',8);
    line([1 2.25],[height height],'Color','k')
elseif signrank(yValues1,yValues2) <= 0.05
    h = text(1.7,height+.01,'*');
    set(h,'Color','k','HorizontalAlignment', 'center','FontSize',8);
    line([1 2.25],[height height],'Color','k')
elseif signrank(yValues1,yValues2) > 0.05
    h = text(1.7,height+.01,'ns');
    set(h,'Color','k','VerticalAlignment','bottom','HorizontalAlignment', 'center','FontSize',8);
    line([1 2.25],[height height],'Color','k')
end

if signrank(yValues3,yValues4) <= 0.0001
    h = text(3.375,height+.01,'****');
    set(h,'Color','k','HorizontalAlignment', 'center','FontSize',8);
    line([2.75 4],[height height],'Color','k')
elseif signrank(yValues3, yValues4) <= 0.001
    h = text(3.375,height+.01,'****');
    set(h,'Color','k','HorizontalAlignment', 'center','FontSize',8);
    line([2.75 4],[height height],'Color','k')
elseif signrank(yValues3,yValues4) <= 0.01
    h = text(3.375,height+.01,'**');
    set(h,'Color','k','HorizontalAlignment', 'center','FontSize',8);
    line([2.75 4],[height height],'Color','k')
elseif signrank(yValues3,yValues4) <= 0.05
    h = text(3.375,height+.01,'*');
    set(h,'Color','k','HorizontalAlignment', 'center','FontSize',8);
    line([2.75 4],[height height],'Color','k')
elseif signrank(yValues3,yValues4) > 0.05
    h = text(3.375,height+.01,'ns');
    set(h,'Color','k','VerticalAlignment','bottom','HorizontalAlignment', 'center','FontSize',8);
    line([2.75 4],[height height],'Color','k')
end

for g = 1:24
    lc(g).Color(4) = .25;
    li(g).Color(4) = .25;
end
set(gca,'tickdir','out')
title('Pre-trial epoch (-500–0 ms)')
set(gca, 'XTick', [1 2.25 2.75 4])
set(gca, 'XTickLabels', {'hit','miss','hit','miss'})
xlabel('Contra cells                Ipsi cells')
ylabel('\DeltaResponse (contra – ipsi trials)')

subplot(1,3,2)
xlim([.5 4.5])
ylim([-y y]); 
t1 = 21;
t2 = 26;
jitter = ((rand(1,24)-.5)/15);
% jitter = 0;
xPositions1 = jitter+ones(1,24);
xPositions2 = jitter+ones(1,24)*2.25;
xPositions3 = jitter+ones(1,24)*2.75;
xPositions4 = jitter+ones(1,24)*4;
yValues1 = mean(hitResps_contraS(:,t1:t2),2)';
yValues2 = mean(missResps_contraS(:,t1:t2),2)';
yValues3 = mean(hitResps_ipsiS(:,t1:t2),2)';
yValues4 = mean(missResps_ipsiS(:,t1:t2),2)';
lc = line([xPositions1; xPositions2],[yValues1; yValues2],'Color',[.5 .5 .5]);
li = line([xPositions3; xPositions4],[yValues3; yValues4],'Color',[.5 .5 .5]);
hold on;
plot(xPositions1,yValues1,'o','MarkerFaceColor',[0 .4 1],'MarkerEdgeColor','w')
plot(xPositions2,yValues2,'o','MarkerFaceColor',[0 .4 1],'MarkerEdgeColor','w')
plot(xPositions3,yValues3,'o','MarkerFaceColor',[1 0 0],'MarkerEdgeColor','w')
plot(xPositions4,yValues4,'o','MarkerFaceColor',[1 0 0],'MarkerEdgeColor','w')

if signrank(yValues1,yValues2) <= 0.0001
    h = text(1.7,height+.01,'****');
    set(h,'Color','k','HorizontalAlignment', 'center','FontSize',8);
    line([1 2.25],[height height],'Color','k')
elseif signrank(yValues1,yValues2) <= 0.001
    h = text(1.7,height+.01,'***');
    set(h,'Color','k','HorizontalAlignment', 'center','FontSize',8);
    line([1 2.25],[height height],'Color','k')
elseif signrank(yValues1,yValues2) <= 0.01
    h = text(1.7,height+.01,'**');
    set(h,'Color','k','HorizontalAlignment', 'center','FontSize',8);
    line([1 2.25],[height height],'Color','k')
elseif signrank(yValues1,yValues2) <= 0.05
    h = text(1.7,height+.01,'*');
    set(h,'Color','k','HorizontalAlignment', 'center','FontSize',8);
    line([1 2.25],[height height],'Color','k')
elseif signrank(yValues1,yValues2) > 0.05
    h = text(1.7,height+.01,'ns');
    set(h,'Color','k','VerticalAlignment','bottom','HorizontalAlignment', 'center','FontSize',8);
    line([1 2.25],[height height],'Color','k')
end

if signrank(yValues3,yValues4) <= 0.0001
    h = text(3.375,height+.01,'****');
    set(h,'Color','k','HorizontalAlignment', 'center','FontSize',8);
    line([2.75 4],[height height],'Color','k')
elseif signrank(yValues3, yValues4) <= 0.001
    h = text(3.375,height+.01,'****');
    set(h,'Color','k','HorizontalAlignment', 'center','FontSize',8);
    line([2.75 4],[height height],'Color','k')
elseif signrank(yValues3,yValues4) <= 0.01
    h = text(3.375,height+.01,'**');
    set(h,'Color','k','HorizontalAlignment', 'center','FontSize',8);
    line([2.75 4],[height height],'Color','k')
elseif signrank(yValues3,yValues4) <= 0.05
    h = text(3.375,height+.01,'*');
    set(h,'Color','k','HorizontalAlignment', 'center','FontSize',8);
    line([2.75 4],[height height],'Color','k')
elseif signrank(yValues3,yValues4) > 0.05
    h = text(3.375,height+.01,'ns');
    set(h,'Color','k','VerticalAlignment','bottom','HorizontalAlignment', 'center','FontSize',8);
    line([2.75 4],[height height],'Color','k')
end


for g = 1:24
    lc(g).Color(4) = .25;
    li(g).Color(4) = .25;
end
set(gca,'tickdir','out')
title('Delay epoch (0–500 ms)')
set(gca, 'XTick', [1 2.25 2.75 4])
set(gca, 'XTickLabels', {'corr.','err.','corr.','err.'})
xlabel('Contra cells                Ipsi cells')
ylabel('\DeltaResponse (contra – ipsi trials)')

subplot(1,3,3)
xlim([.5 4.5])
ylim([-y y]); 
t1 = 28;
t2 = 35;
jitter = ((rand(1,24)-.5)/15);
% jitter = 0;
xPositions1 = jitter+ones(1,24);
xPositions2 = jitter+ones(1,24)*2.25;
xPositions3 = jitter+ones(1,24)*2.75;
xPositions4 = jitter+ones(1,24)*4;
yValues1 = mean(hitResps_contraS(:,t1:t2),2)';
yValues2 = mean(missResps_contraS(:,t1:t2),2)';
yValues3 = mean(hitResps_ipsiS(:,t1:t2),2)';
yValues4 = mean(missResps_ipsiS(:,t1:t2),2)';
lc = line([xPositions1; xPositions2],[yValues1; yValues2],'Color',[.5 .5 .5]);
li = line([xPositions3; xPositions4],[yValues3; yValues4],'Color',[.5 .5 .5]);
hold on;
plot(xPositions1,yValues1,'o','MarkerFaceColor',[0 .4 1],'MarkerEdgeColor','w')
plot(xPositions2,yValues2,'o','MarkerFaceColor',[0 .4 1],'MarkerEdgeColor','w')
plot(xPositions3,yValues3,'o','MarkerFaceColor',[1 0 0],'MarkerEdgeColor','w')
plot(xPositions4,yValues4,'o','MarkerFaceColor',[1 0 0],'MarkerEdgeColor','w')

if signrank(yValues1,yValues2) <= 0.0001
    h = text(1.7,height+.01,'****');
    set(h,'Color','k','HorizontalAlignment', 'center','FontSize',8);
    line([1 2.25],[height height],'Color','k')
elseif signrank(yValues1,yValues2) <= 0.001
    h = text(1.7,height+.01,'***');
    set(h,'Color','k','HorizontalAlignment', 'center','FontSize',8);
    line([1 2.25],[height height],'Color','k')
elseif signrank(yValues1,yValues2) <= 0.01
    h = text(1.7,height+.01,'**');
    set(h,'Color','k','HorizontalAlignment', 'center','FontSize',8);
    line([1 2.25],[height height],'Color','k')
elseif signrank(yValues1,yValues2) <= 0.05
    h = text(1.7,height+.01,'*');
    set(h,'Color','k','HorizontalAlignment', 'center','FontSize',8);
    line([1 2.25],[height height],'Color','k')
elseif signrank(yValues1,yValues2) > 0.05
    h = text(1.7,height+.01,'ns');
    set(h,'Color','k','VerticalAlignment','bottom','HorizontalAlignment', 'center','FontSize',8);
    line([1 2.25],[height height],'Color','k')
end

if signrank(yValues3,yValues4) <= 0.0001
    h = text(3.375,height+.01,'****');
    set(h,'Color','k','HorizontalAlignment', 'center','FontSize',8);
    line([2.75 4],[height height],'Color','k')
elseif signrank(yValues3, yValues4) <= 0.001
    h = text(3.375,height+.01,'****');
    set(h,'Color','k','HorizontalAlignment', 'center','FontSize',8);
    line([2.75 4],[height height],'Color','k')
elseif signrank(yValues3,yValues4) <= 0.01
    h = text(3.375,height+.01,'**');
    set(h,'Color','k','HorizontalAlignment', 'center','FontSize',8);
    line([2.75 4],[height height],'Color','k')
elseif signrank(yValues3,yValues4) <= 0.05
    h = text(3.375,height+.01,'*');
    set(h,'Color','k','HorizontalAlignment', 'center','FontSize',8);
    line([2.75 4],[height height],'Color','k')
elseif signrank(yValues3,yValues4) > 0.05
    h = text(3.375,height+.01,'ns');
    set(h,'Color','k','VerticalAlignment','bottom','HorizontalAlignment', 'center','FontSize',8);
    line([2.75 4],[height height],'Color','k')
end
for g = 1:24
    lc(g).Color(4) = .25;
    li(g).Color(4) = .25;
end
set(gca,'tickdir','out')
title('Response epoch (800–1400 ms)')
set(gca, 'XTick', [1 2.25 2.75 4])
set(gca, 'XTickLabels', {'corr.','err.','corr.','err.'})
xlabel('Contra cells                Ipsi cells')
ylabel('\DeltaResponse (contra – ipsi trials)')

%%

figure;
subplot(1,3,1)
xlim([.5 4.5])
ylim([-y y]); 
t1 = 15;
t2 = 20;
jitter = ((rand(1,24)-.5)/15);
% jitter = 0;
xPositions1 = jitter+ones(1,24);
xPositions2 = jitter+ones(1,24)*2.25;
xPositions3 = jitter+ones(1,24)*2.75;
xPositions4 = jitter+ones(1,24)*4;
yValues1 = mean(hitResps_hit(:,t1:t2),2)';
yValues2 = mean(missResps_hit(:,t1:t2),2)';
yValues3 = mean(hitResps_miss(:,t1:t2),2)';
yValues4 = mean(missResps_miss(:,t1:t2),2)';
lc = line([xPositions1; xPositions2],[yValues1; yValues2],'Color',[.5 .5 .5]);
li = line([xPositions3; xPositions4],[yValues3; yValues4],'Color',[.5 .5 .5]);
hold on;
plot(xPositions1,yValues1,'o','MarkerFaceColor',[0 .5 .1],'MarkerEdgeColor','w')
plot(xPositions2,yValues2,'o','MarkerFaceColor',[0 .5 .1],'MarkerEdgeColor','w')
plot(xPositions3,yValues3,'o','MarkerFaceColor',[.75 0 0],'MarkerEdgeColor','w')
plot(xPositions4,yValues4,'o','MarkerFaceColor',[.75 0 0],'MarkerEdgeColor','w')

if signrank(yValues1,yValues2) <= 0.0001
    h = text(1.7,height+.01,'****');
    set(h,'Color','k','HorizontalAlignment', 'center','FontSize',8);
    line([1 2.25],[height height],'Color','k')
elseif signrank(yValues1,yValues2) <= 0.001
    h = text(1.7,height+.01,'***');
    set(h,'Color','k','HorizontalAlignment', 'center','FontSize',8);
    line([1 2.25],[height height],'Color','k')
elseif signrank(yValues1,yValues2) <= 0.01
    h = text(1.7,height+.01,'**');
    set(h,'Color','k','HorizontalAlignment', 'center','FontSize',8);
    line([1 2.25],[height height],'Color','k')
elseif signrank(yValues1,yValues2) <= 0.05
    h = text(1.7,height+.01,'*');
    set(h,'Color','k','HorizontalAlignment', 'center','FontSize',8);
    line([1 2.25],[height height],'Color','k')
elseif signrank(yValues1,yValues2) > 0.05
    h = text(1.7,height+.01,'ns');
    set(h,'Color','k','VerticalAlignment','bottom','HorizontalAlignment', 'center','FontSize',8);
    line([1 2.25],[height height],'Color','k')
end

if signrank(yValues3,yValues4) <= 0.0001
    h = text(3.375,height+.01,'****');
    set(h,'Color','k','HorizontalAlignment', 'center','FontSize',8);
    line([2.75 4],[height height],'Color','k')
elseif signrank(yValues3, yValues4) <= 0.001
    h = text(3.375,height+.01,'****');
    set(h,'Color','k','HorizontalAlignment', 'center','FontSize',8);
    line([2.75 4],[height height],'Color','k')
elseif signrank(yValues3,yValues4) <= 0.01
    h = text(3.375,height+.01,'**');
    set(h,'Color','k','HorizontalAlignment', 'center','FontSize',8);
    line([2.75 4],[height height],'Color','k')
elseif signrank(yValues3,yValues4) <= 0.05
    h = text(3.375,height+.01,'*');
    set(h,'Color','k','HorizontalAlignment', 'center','FontSize',8);
    line([2.75 4],[height height],'Color','k')
elseif signrank(yValues3,yValues4) > 0.05
    h = text(3.375,height+.01,'ns');
    set(h,'Color','k','VerticalAlignment','bottom','HorizontalAlignment', 'center','FontSize',8);
    line([2.75 4],[height height],'Color','k')
end

for g = 1:24
    lc(g).Color(4) = .25;
    li(g).Color(4) = .25;
end
set(gca,'tickdir','out')
title('Pre-trial epoch (-500–0 ms)')
set(gca, 'XTick', [1 2.25 2.75 4])
set(gca, 'XTickLabels', {'corr.','err.','corr.','err.'})
xlabel('Correct cells                Error cells')
ylabel('\DeltaResponse (contra – ipsi trials)')

subplot(1,3,2)
xlim([.5 4.5])
ylim([-y y]); 
t1 = 21;
t2 = 26;
jitter = ((rand(1,24)-.5)/15);
% jitter = 0;
xPositions1 = jitter+ones(1,24);
xPositions2 = jitter+ones(1,24)*2.25;
xPositions3 = jitter+ones(1,24)*2.75;
xPositions4 = jitter+ones(1,24)*4;
yValues1 = mean(hitResps_hit(:,t1:t2),2)';
yValues2 = mean(missResps_hit(:,t1:t2),2)';
yValues3 = mean(hitResps_miss(:,t1:t2),2)';
yValues4 = mean(missResps_miss(:,t1:t2),2)';
lc = line([xPositions1; xPositions2],[yValues1; yValues2],'Color',[.5 .5 .5]);
li = line([xPositions3; xPositions4],[yValues3; yValues4],'Color',[.5 .5 .5]);
hold on;
plot(xPositions1,yValues1,'o','MarkerFaceColor',[0 .5 .1],'MarkerEdgeColor','w')
plot(xPositions2,yValues2,'o','MarkerFaceColor',[0 .5 .1],'MarkerEdgeColor','w')
plot(xPositions3,yValues3,'o','MarkerFaceColor',[.75 0 0],'MarkerEdgeColor','w')
plot(xPositions4,yValues4,'o','MarkerFaceColor',[.75 0 0],'MarkerEdgeColor','w')

if signrank(yValues1,yValues2) <= 0.0001
    h = text(1.7,height+.01,'****');
    set(h,'Color','k','HorizontalAlignment', 'center','FontSize',8);
    line([1 2.25],[height height],'Color','k')
elseif signrank(yValues1,yValues2) <= 0.001
    h = text(1.7,height+.01,'***');
    set(h,'Color','k','HorizontalAlignment', 'center','FontSize',8);
    line([1 2.25],[height height],'Color','k')
elseif signrank(yValues1,yValues2) <= 0.01
    h = text(1.7,height+.01,'**');
    set(h,'Color','k','HorizontalAlignment', 'center','FontSize',8);
    line([1 2.25],[height height],'Color','k')
elseif signrank(yValues1,yValues2) <= 0.05
    h = text(1.7,height+.01,'*');
    set(h,'Color','k','HorizontalAlignment', 'center','FontSize',8);
    line([1 2.25],[height height],'Color','k')
elseif signrank(yValues1,yValues2) > 0.05
    h = text(1.7,height+.01,'ns');
    set(h,'Color','k','VerticalAlignment','bottom','HorizontalAlignment', 'center','FontSize',8);
    line([1 2.25],[height height],'Color','k')
end

if signrank(yValues3,yValues4) <= 0.0001
    h = text(3.375,height+.01,'****');
    set(h,'Color','k','HorizontalAlignment', 'center','FontSize',8);
    line([2.75 4],[height height],'Color','k')
elseif signrank(yValues3, yValues4) <= 0.001
    h = text(3.375,height+.01,'****');
    set(h,'Color','k','HorizontalAlignment', 'center','FontSize',8);
    line([2.75 4],[height height],'Color','k')
elseif signrank(yValues3,yValues4) <= 0.01
    h = text(3.375,height+.01,'**');
    set(h,'Color','k','HorizontalAlignment', 'center','FontSize',8);
    line([2.75 4],[height height],'Color','k')
elseif signrank(yValues3,yValues4) <= 0.05
    h = text(3.375,height+.01,'*');
    set(h,'Color','k','HorizontalAlignment', 'center','FontSize',8);
    line([2.75 4],[height height],'Color','k')
elseif signrank(yValues3,yValues4) > 0.05
    h = text(3.375,height+.01,'ns');
    set(h,'Color','k','VerticalAlignment','bottom','HorizontalAlignment', 'center','FontSize',8);
    line([2.75 4],[height height],'Color','k')
end


for g = 1:24
    lc(g).Color(4) = .25;
    li(g).Color(4) = .25;
end
set(gca,'tickdir','out')
title('Delay epoch (0–500 ms)')
set(gca, 'XTick', [1 2.25 2.75 4])
set(gca, 'XTickLabels', {'corr.','err.','corr.','err.'})
xlabel('Correct cells                Error cells')
ylabel('\DeltaResponse (contra – ipsi trials)')

subplot(1,3,3)
xlim([.5 4.5])
ylim([-y y]);
t1 = 28;
t2 = 35;
jitter = ((rand(1,24)-.5)/15);
% jitter = 0;
xPositions1 = jitter+ones(1,24);
xPositions2 = jitter+ones(1,24)*2.25;
xPositions3 = jitter+ones(1,24)*2.75;
xPositions4 = jitter+ones(1,24)*4;
yValues1 = mean(hitResps_hit(:,t1:t2),2)';
yValues2 = mean(missResps_hit(:,t1:t2),2)';
yValues3 = mean(hitResps_miss(:,t1:t2),2)';
yValues4 = mean(missResps_miss(:,t1:t2),2)';
lc = line([xPositions1; xPositions2],[yValues1; yValues2],'Color',[.5 .5 .5]);
li = line([xPositions3; xPositions4],[yValues3; yValues4],'Color',[.5 .5 .5]);
hold on;
plot(xPositions1,yValues1,'o','MarkerFaceColor',[0 .5 .1],'MarkerEdgeColor','w')
plot(xPositions2,yValues2,'o','MarkerFaceColor',[0 .5 .1],'MarkerEdgeColor','w')
plot(xPositions3,yValues3,'o','MarkerFaceColor',[.75 0 0],'MarkerEdgeColor','w')
plot(xPositions4,yValues4,'o','MarkerFaceColor',[.75 0 0],'MarkerEdgeColor','w')

if signrank(yValues1,yValues2) <= 0.0001
    h = text(1.7,height+.01,'****');
    set(h,'Color','k','HorizontalAlignment', 'center','FontSize',8);
    line([1 2.25],[height height],'Color','k')
elseif signrank(yValues1,yValues2) <= 0.001
    h = text(1.7,height+.01,'***');
    set(h,'Color','k','HorizontalAlignment', 'center','FontSize',8);
    line([1 2.25],[height height],'Color','k')
elseif signrank(yValues1,yValues2) <= 0.01
    h = text(1.7,height+.01,'**');
    set(h,'Color','k','HorizontalAlignment', 'center','FontSize',8);
    line([1 2.25],[height height],'Color','k')
elseif signrank(yValues1,yValues2) <= 0.05
    h = text(1.7,height+.01,'*');
    set(h,'Color','k','HorizontalAlignment', 'center','FontSize',8);
    line([1 2.25],[height height],'Color','k')
elseif signrank(yValues1,yValues2) > 0.05
    h = text(1.7,height+.01,'ns');
    set(h,'Color','k','VerticalAlignment','bottom','HorizontalAlignment', 'center','FontSize',8);
    line([1 2.25],[height height],'Color','k')
end

if signrank(yValues3,yValues4) <= 0.0001
    h = text(3.375,height+.01,'****');
    set(h,'Color','k','HorizontalAlignment', 'center','FontSize',8);
    line([2.75 4],[height height],'Color','k')
elseif signrank(yValues3, yValues4) <= 0.001
    h = text(3.375,height+.01,'****');
    set(h,'Color','k','HorizontalAlignment', 'center','FontSize',8);
    line([2.75 4],[height height],'Color','k')
elseif signrank(yValues3,yValues4) <= 0.01
    h = text(3.375,height+.01,'**');
    set(h,'Color','k','HorizontalAlignment', 'center','FontSize',8);
    line([2.75 4],[height height],'Color','k')
elseif signrank(yValues3,yValues4) <= 0.05
    h = text(3.375,height+.01,'*');
    set(h,'Color','k','HorizontalAlignment', 'center','FontSize',8);
    line([2.75 4],[height height],'Color','k')
elseif signrank(yValues3,yValues4) > 0.05
    h = text(3.375,height+.01,'ns');
    set(h,'Color','k','VerticalAlignment','bottom','HorizontalAlignment', 'center','FontSize',8);
    line([2.75 4],[height height],'Color','k')
end
for g = 1:24
    lc(g).Color(4) = .25;
    li(g).Color(4) = .25;
end
set(gca,'tickdir','out')
title('Response epoch (800–1400 ms)')
set(gca, 'XTick', [1 2.25 2.75 4])
set(gca, 'XTickLabels', {'corr.','err.','corr.','err.'})
xlabel('Correct cells                Error cells')
ylabel('\DeltaResponse (contra – ipsi trials)')

%%
ms = 40;
figure;
subplot(1,3,2)
scatter(numContraMovCells,numIpsiMovCells,ms,'MarkerEdgeColor','w','MarkerFaceColor','k')
line([0 1000],[0 1000],'LineStyle','--','Color',[.5 .5 .5]);
axis square
set(gca,'tickdir','out')
xlabel('No. contraversive move cells')
ylabel('No. ipsiversive move cells')

subplot(1,3,1)
scatter(numContraStimCells,numIpsiStimCells,ms,'MarkerEdgeColor','w','MarkerFaceColor','k')
line([0 400],[0 400],'LineStyle','--','Color',[.5 .5 .5]);
axis square
set(gca,'tickdir','out')
xlabel('No. contralateral stim cells')
ylabel('No. ipsilateral stim cells')

subplot(1,3,3)
scatter(numHitCells,numMissCells,ms,'MarkerEdgeColor','w','MarkerFaceColor','k')
line([0 1800],[0 1800],'LineStyle','--','Color',[.5 .5 .5]);
axis square
set(gca,'tickdir','out')
ylabel('No. error cells')
xlabel('No. reward cells')





