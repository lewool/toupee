function neuralData = getSingleCellStats(expInfo, behavioralData, neuralData)

labels = {'stimOnset' 'leftStim' 'rightStim' 'moveOnset' 'choice' 'outcome' 'block' 'value'};

 %%   
et = behavioralData;
contrasts = getUniqueContrasts(expInfo);
nt = length(et.eventTimes(1).daqTime);

trialCorrectChoice = expInfo.block.events.correctResponseValues(1:nt);
trueStimuli = expInfo.block.events.contrastValues(1:nt);
trueStimuli(trueStimuli == 0) = eps;
trueStimuli(abs(trueStimuli) < .05) = ...
    trueStimuli(abs(trueStimuli) < .05).* trialCorrectChoice(abs(trueStimuli) < .05);

trueChoices = et.wheelMoves.epochs(5).moveDir;
[~, leftStimTrials] = selectCondition(expInfo, contrasts(contrasts<0), et, ...
    initTrialConditions('movementTime','late','specificRTs',[.8 Inf]));
[~, rightStimTrials] = selectCondition(expInfo, contrasts(contrasts>0), et, ...
    initTrialConditions('movementTime','late','specificRTs',[.8 Inf]));
trueBlock = expInfo.block.events.highRewardSideValues(1:nt);
trueValue(trueBlock .* sign(trueStimuli) > 0) = 2;
trueValue(trueBlock .* sign(trueStimuli) < 0) = 1;


[~, stimTrials] = selectCondition(expInfo, contrasts, et, ...
        initTrialConditions('preStimMovement','quiescent','movementTime','late','specificRTs',[.8 Inf]));
[~, zeroStimTrials] = selectCondition(expInfo, 0, et, ...
        initTrialConditions('movementTime','late','specificRTs',[.8 Inf]));
[~, correctTrials] = selectCondition(expInfo, 0, et, initTrialConditions('responseType','correct'));
[~, incorrectTrials] = selectCondition(expInfo, 0, et, initTrialConditions('responseType','incorrect'));
[~, zeroTrials] = selectCondition(expInfo, 0, et, initTrialConditions('movementTime','all'));
nonnanTrials = find(~isnan(trueChoices));

%get responses
[baselineResps, stimResps, pmovResps, movResps, rewResps, preCueResps] = getEpochResps(neuralData.eta);

%% stimulus

%stim onset (general)
for c = 1:length(movResps)
    [~, h] = signrank(stimResps(stimTrials,c),baselineResps(stimTrials,c));
    hTest(c,1) = double(h);
end

% stimulus side
for c = 1:size(baselineResps,2)
    [~, hL] = ranksum(stimResps(leftStimTrials,c),stimResps(zeroStimTrials,c),'tail','right');
    [~, hR] = ranksum(stimResps(rightStimTrials,c),stimResps(zeroStimTrials,c),'tail','right');
    hTest(c,2) = double(hL);
    hTest(c,3) = double(hR);
end

%% movement & choice

whichTrials = intersect(nonnanTrials,zeroTrials);

% movement
for c = 1:length(movResps)
    [~,h] = signrank(movResps(whichTrials,c),pmovResps(whichTrials,c));
    hTest(c,4) = double(h);
end

% choice
trueChoices_0 = trueChoices(whichTrials);
movResps_0 = movResps(whichTrials,:);

trimLength = 20;
for c = 1:length(movResps)
    clear shiftDiff
    trueIdx = trimLength+1;
    mov0_trunc = movResps_0(trimLength+1:end - trimLength,c);
    for l = 1:trimLength*2+1
        ss = l;
        es = length(trueChoices_0) - (trimLength*2-l+1);
        shiftDiff(l) = mean(mov0_trunc(trueChoices_0(ss:es) == 1)) - mean(mov0_trunc(trueChoices_0(ss:es) == -1));
    end
    
    trueVal = shiftDiff(trueIdx);
    pseudoVals = shiftDiff([1:trimLength,trimLength+2:trimLength*2+1]);
    [~,h] = ranksum(trueVal,pseudoVals,'method','exact','tail','both');
    
    if trueVal < median(pseudoVals)
        hTest(c,5) = double(-h);
    else
        hTest(c,5) = double(h);
    end
end

%% reward
for c = 1:size(baselineResps,2)
    [~,h] = ranksum(rewResps(correctTrials,c),rewResps(incorrectTrials,c),'tail','both');
    if mean(rewResps(correctTrials,c)) > mean(rewResps(incorrectTrials,c))
        hTest(c,6) = double(h);
    else
        hTest(c,6) = double(-h);
    end
end

%% block & value

whichTrials = intersect(stimTrials,nonnanTrials);
whichIdx = false(1,nt);
whichIdx(whichTrials) = true;
whichResps = pmovResps;

%block
np = 1000;
pseudoBlocks = generatePseudoBlocks(expInfo,np,'fixed');

rbTrials = intersect(whichTrials, find(trueBlock > 0));
lbTrials = intersect(whichTrials, find(trueBlock < 0));

for p = 1:np
    rbTrials_pseudo{p} = find(whichIdx.*(pseudoBlocks(p,:) > 0));
    lbTrials_pseudo{p} = find(whichIdx.*(pseudoBlocks(p,:) < 0));
end

for c = 1:length(baselineResps)
    trueDiff = mean(whichResps(rbTrials,c)) - mean(whichResps(lbTrials,c));
    for p = 1:np
        pseudoDiff(p) = mean(whichResps(rbTrials_pseudo{p},c)) - mean(whichResps(lbTrials_pseudo{p},c));
    end
    [~,h] = ranksum(trueDiff,pseudoDiff,'method','exact','tail','both');
    
    if trueDiff > 0
        hTest(c,7) = double(h);
    else
        hTest(c,7) = double(-h);
    end
        
end

% value


hvTrials = intersect(whichTrials, find(trueValue == 2));
lvTrials = intersect(whichTrials, find(trueValue == 1));

for p = 1:np
    pseudoValues(p, (pseudoBlocks(p,:) .* sign(trueStimuli) > 0)) = 2;
    pseudoValues(p, (pseudoBlocks(p,:) .* sign(trueStimuli) < 0)) = 1;
end
for p = 1:np
    hvTrials_pseudo{p} = find(whichIdx.*(pseudoValues(p,:) == 2));
    lvTrials_pseudo{p} = find(whichIdx.*(pseudoValues(p,:) == 1));
end

for c = 1:length(whichResps)
    trueDiff = mean(whichResps(hvTrials,c)) - mean(whichResps(lvTrials,c));
    for p = 1:np
        pseudoDiff(p) = mean(whichResps(hvTrials_pseudo{p},c)) - mean(whichResps(lvTrials_pseudo{p},c));
    end
    [~,h] = ranksum(trueDiff,pseudoDiff,'method','exact','tail','both');
    
    if trueDiff > 0
        hTest(c,8) = double(h);
    else
        hTest(c,8) = double(-h);
    end 
end

%% Q


% [Q, Q_pseudo] = computeQ(behavioralData, expInfo);
% 
% whichTrials = intersect(stimTrials,nonnanTrials);
% whichResps = pmovResps;
% 
% for c = 1:length(whichResps)
%     trueCorrQc = corr(Q.Qc(whichTrials),whichResps(whichTrials,c));
%     trueCorrQtotal = corr(Q.Qtotal(whichTrials),whichResps(whichTrials,c));
%     trueCorrQrel = corr(Q.Qrel(whichTrials),whichResps(whichTrials,c));
%     
%     for p = 1:np
%         pseudoCorrQc(p) = corr(Q_pseudo.Qc(whichTrials,p),whichResps(whichTrials,c));
%         pseudoCorrQtotal(p) = corr(Q_pseudo.Qtotal(whichTrials,p),whichResps(whichTrials,c));
%         pseudoCorrQrel(p) = corr(Q_pseudo.Qrel(whichTrials,p),whichResps(whichTrials,c));
%     end
%        
%     [~,hQc] = ranksum(trueCorrQc,pseudoCorrQc,'method','exact','tail','both');
%     if trueCorrQc > 0
%         hTest(c,9) = double(hQc);
%     else
%         hTest(c,9) = double(-hQc);
%     end 
%     
%     [~,hQtotal] = ranksum(trueCorrQtotal,pseudoCorrQtotal,'method','exact','tail','both');
%     if trueCorrQtotal > 0
%         hTest(c,10) = double(hQtotal);
%     else
%         hTest(c,10) = double(-hQtotal);
%     end
%     
%     [~,hQrel] = ranksum(trueCorrQrel,pseudoCorrQrel,'method','exact','tail','both');
%     if trueCorrQrel > 0
%         hTest(c,11) = double(hQrel);
%     else
%         hTest(c,11) = double(-hQrel);
%     end
% end

%%
neuralData.stats = [];
neuralData.stats.labels = labels;
neuralData.stats.alpha = 0.05;
neuralData.stats.hTest = hTest;








%%
