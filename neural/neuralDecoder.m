function [gof, gof_pseudo, mi, mip] = neuralDecoder(expInfo, behavioralData, neuralData, Ytype, ETA)

%% Predict a task variable from neural activity
% This takes a 1D vector of neural activity at a particular trial timepoint (X) and uses
% this in lasso logistic regression to try and predict a Y outcome variable
% on that trial. A model is fitted for every timepoint (n = 41)

% Data is split into thirds. 2/3 goes into the cvglmnet protocol to train & pick the
% lambda, then the other 1/3 is held back to evaluate the model.




%%
    
nt = length(behavioralData.eventTimes(1).daqTime);

%% generate Y outcome vector

% extract stimulus, choice, feedback, value, and block values

trueStimuli = expInfo.block.events.contrastValues(1:nt);
trialCorrectChoice = expInfo.block.events.correctResponseValues(1:nt);
trueBlocks = expInfo.block.events.highRewardSideValues(1:nt);
trueChoices = expInfo.block.events.responseValues(1:nt);
allFeedback = double(expInfo.block.events.feedbackValues(1:nt));
feedbackIdx = find(trueStimuli == 0);
trueFeedback = allFeedback(feedbackIdx);

% assign the 0% stimuli as either 'left' or 'right' depending on the
% preassigned correct choice (not the mouse's choice)
trueStimuli(trueStimuli == 0) = eps;
trueStimuli(abs(trueStimuli) < .05) = ...
    trueStimuli(abs(trueStimuli) < .05).* trialCorrectChoice(abs(trueStimuli) < .05);

trueSide = trueStimuli > 0;

% low rewards are possible on sign-mismatched block and stimulus
% high rewards are possible on sign-matched block and stimulus
% 1 = high, 0 = low
trueValue(trueBlocks.*sign(trueStimuli) == -1) = 0;
trueValue(trueBlocks.*sign(trueStimuli) == 1) = 1;

% fetch RT and maxVel

trueVel = abs(behavioralData.wheelMoves.epochs(5).peakVel(1:nt));
trueRT = behavioralData.wheelMoves.epochs(5).onsetTimes(1:nt) - behavioralData.eventTimes(1).daqTime(1:nt);

%% generate pseudo Y vectors
np = 100;

%generate randomly shuffled stims
pseudoStimuli = nan(np,nt);
for p = 1:np
    pseudoStimuli(p,:) = randsample(trueStimuli,length(trueStimuli));
end

%generate randomly shuffled stims
pseudoSide = nan(np,nt);
for p = 1:np
    pseudoSide(p,:) = randsample(trueSide,length(trueSide));
end

%generate linear-shifted choices
% trimLength = 20;
% for l = 1:trimLength*2+1
%     ss = l;
%     es = length(trueChoices) - (trimLength*2-l+1);
%     shifts(l,:) = trueChoices(ss:es);
% end
% pseudoChoices = shifts([1:trimLength,trimLength+2:trimLength*2+1],:);

%generate randomly shuffled feedbacks (only 0% contrast)
pseudoFeedback = nan(np,length(trueFeedback));
for p = 1:np
    pseudoFeedback(p,:) = randsample(trueFeedback,length(trueFeedback));
end

%generate linear-shifted RTs
trimLength = 20;
for l = 1:trimLength*2+1
    ss = l;
    es = length(trueRT) - (trimLength*2-l+1);
    rt_shifts(l,:) = trueRT(ss:es);
end
pseudoRT = rt_shifts([1:trimLength,trimLength+2:trimLength*2+1],:);

%generate linear-shifted maxVels
trimLength = 20;
for l = 1:trimLength*2+1
    ss = l;
    es = length(trueVel) - (trimLength*2-l+1);
    vel_shifts(l,:) = trueVel(ss:es);
end
pseudoVel = vel_shifts([1:trimLength,trimLength+2:trimLength*2+1],:);

% generate pseudoblocks
pseudoBlocks = nan(np,nt);
pseudoValue = nan(np,nt);
trimLength = 0;
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
    pseudoBlocks(p,:) = b(1:nt);
end

% generate pseudo high/low trials based on the pseudoblocks + true stim
for p = 1:np
    pseudoValue(p,pseudoBlocks(p,:).*sign(trueStimuli) == -1) = 0;
    pseudoValue(p,pseudoBlocks(p,:).*sign(trueStimuli) == 1) = 1;
%     pseudoValue(p,trueFeedback==0) = nan;
%     pseudoValue(p,isnan(pseudoBlocks(p,:))) = nan;
end


%% convert to contra/ipsi

if expInfo.hemisphere > 0
    trueStimuli = -trueStimuli;
    trueChoices = (trueChoices == -1);
    trueBlocks = (trueBlocks == -1);
    pseudoBlocks = pseudoBlocks == -1;
    pseudoStimuli = -pseudoStimuli;
    trueSide = trueSide == 0;
    pseudoSide = pseudoSide == 0;
else
    trueChoices = (trueChoices == 1);
    trueBlocks = (trueBlocks == 1);
    pseudoBlocks = pseudoBlocks == 1;
end

%% generate stimulus one-hot encoding matrix

contrasts = getUniqueContrasts(expInfo);
contrasts(contrasts==0) = [];

stimOneHot = zeros(length(trueStimuli),length(contrasts));
for s=1:length(trueStimuli)
    stimOneHot(s,:) = contrasts == trueStimuli(s);
end
stimOneHot = stimOneHot * 1;

choiceOneHot = zeros(length(trueChoices),1);
for c=1:length(trueChoices)
    if trueChoices(c) == 1
        choiceOneHot(c,1) = 1;
    elseif trueChoices(c) == 0
        choiceOneHot(c,1) = -1;
    end
end

% choiceOneHot = zeros(length(trueChoices),2);
% for c=1:length(trueChoices)
%     if trueChoices(c) == 1
%         choiceOneHot(c,1) = 1;
%     elseif trueChoices(c) == 0
%         choiceOneHot(c,2) = 1;
%     end
% end
%% fit type

timeRange = 6:2:36;

for t = 1:length(timeRange)
    switch Ytype
        case 'side'
            % generally filter trials (e.g., only patient ones)
            [~, whichTrials] = selectCondition(expInfo, getUniqueContrasts(expInfo), behavioralData, ...
                    initTrialConditions('repeatType','random','movementTime','late','specificRTs',[.8 3]));
            testTrials = whichTrials(1:3:end);
            trainTrials = whichTrials;
            trainTrials(ismember(whichTrials,testTrials)) = [];
            
            X(t).train = squeeze(neuralData.eta.alignedResps{ETA}(trainTrials,timeRange(t),:));
            X(t).test = squeeze(neuralData.eta.alignedResps{ETA}(testTrials,timeRange(t),:));
    
            Y(t).true.train = trueSide(trainTrials);
            Y(t).true.test = trueSide(testTrials);
            Y(t).true.dumb = ones(1,length(testTrials))*mean(trueSide(trainTrials));
            Y(t).pseudo.train = pseudoSide(:,trainTrials);
            Y(t).pseudo.test = pseudoSide(:,testTrials);
            Y(t).pseudo.dumb = ones(1,length(testTrials)).*mean(pseudoSide(:,trainTrials),2);
            family = 'binomial';
            
        case 'stimulus'
            % generally filter trials (e.g., only patient ones)
            [~, whichTrials] = selectCondition(expInfo, getUniqueContrasts(expInfo), behavioralData, ...
                    initTrialConditions('repeatType','random','movementTime','late','specificRTs',[.8 3]));
            testTrials = whichTrials(1:3:end);
            trainTrials = whichTrials;
            trainTrials(ismember(whichTrials,testTrials)) = [];
            
            X(t).train = squeeze(neuralData.eta.alignedResps{ETA}(trainTrials,timeRange(t),:));
            X(t).test = squeeze(neuralData.eta.alignedResps{ETA}(testTrials,timeRange(t),:));
    
            Y(t).true.train = trueStimuli(trainTrials);
            Y(t).true.test = trueStimuli(testTrials);
            Y(t).true.dumb = ones(1,length(testTrials))*mean(trueStimuli(trainTrials));
            Y(t).pseudo.train = pseudoStimuli(:,trainTrials);
            Y(t).pseudo.test = pseudoStimuli(:,testTrials);
            Y(t).pseudo.dumb = ones(1,length(testTrials)).*mean(pseudoStimuli(:,trainTrials),2);
            family = 'gaussian';
            
        case 'choice'
            [~, whichTrials] = selectCondition(expInfo, getUniqueContrasts(expInfo), behavioralData, ...
                    initTrialConditions('repeatType','random','movementTime','late','specificRTs',[.8 3]));
            whichChoices = trueChoices(whichTrials);
            %generate linear-shifted choices
            trimLength = 20;
            for l = 1:trimLength*2+1
                ss = l;
                es = length(whichChoices) - (trimLength*2-l+1);
                shifts(l,:) = whichChoices(ss:es);
            end
            pseudoChoices = shifts([1:trimLength,trimLength+2:trimLength*2+1],:);
            testTrials = whichTrials(trimLength+1:3:end-trimLength);
            trainTrials = whichTrials(trimLength+1:end-trimLength);
            trainTrials(ismember(whichTrials(trimLength+1:end-trimLength),testTrials)) = [];
            trimLength = 20;
            truncIdx = trimLength+1:nt-trimLength;
            ptestTrials = 1:3:length(pseudoChoices);
            ptrainTrials = 1:length(pseudoChoices);
            ptrainTrials(ptestTrials) = [];
            
            X(t).train = squeeze(neuralData.eta.alignedResps{ETA}(trainTrials,timeRange(t),:));
            X(t).test = squeeze(neuralData.eta.alignedResps{ETA}(testTrials,timeRange(t),:));
            
            Y(t).true.train = trueChoices(trainTrials);
            Y(t).true.test = trueChoices(testTrials);
            Y(t).true.dumb = ones(1,length(testTrials))*mean(trueChoices(trainTrials));
            Y(t).pseudo.train = pseudoChoices(:,ptrainTrials);
            Y(t).pseudo.test = pseudoChoices(:,ptestTrials);
            Y(t).pseudo.dumb = ones(1,length(ptestTrials)).*mean(pseudoChoices(:,ptrainTrials),2);
            family = 'binomial';
        
        case 'prevChoice'
            [~, whichTrials] = selectCondition(expInfo, getUniqueContrasts(expInfo), behavioralData, ...
                    initTrialConditions('repeatType','random','movementTime','late','specificRTs',[.8 3]));
            whichChoices = trueChoices(whichTrials);
            %generate linear-shifted choices
            trimLength = 20;
            for l = 1:trimLength*2+1
                ss = l;
                es = length(whichChoices) - (trimLength*2-l+1);
                shifts(l,:) = whichChoices(ss:es);
            end
            pseudoChoices = shifts([1:trimLength,trimLength+2:trimLength*2+1],:);
            testTrials = whichTrials(trimLength+1:3:end-trimLength);
            trainTrials = whichTrials(trimLength+1:end-trimLength);
            trainTrials(ismember(whichTrials(trimLength+1:end-trimLength),testTrials)) = [];
            trimLength = 20;
            truncIdx = trimLength+1:nt-trimLength;
            ptestTrials = 1:3:length(pseudoChoices);
            ptrainTrials = 1:length(pseudoChoices);
            ptrainTrials(ptestTrials) = [];
            
            X(t).train = squeeze(neuralData.eta.alignedResps{ETA}(trainTrials+1,timeRange(t),:));
            X(t).test = squeeze(neuralData.eta.alignedResps{ETA}(testTrials+1,timeRange(t),:));
            
            Y(t).true.train = trueChoices(trainTrials);
            Y(t).true.test = trueChoices(testTrials);
            Y(t).pseudo.train = pseudoChoices(:,ptrainTrials);
            Y(t).pseudo.test = pseudoChoices(:,ptestTrials);
            family = 'binomial';
        
        case 'choiceIndStim'
            [~, whichTrials] = selectCondition(expInfo, getUniqueContrasts(expInfo), behavioralData, ...
                    initTrialConditions('repeatType','random','movementTime','late','specificRTs',[.8 3]));
            whichChoices = trueChoices(whichTrials);
            whichStimOneHot = stimOneHot(whichTrials,:);
            %generate linear-shifted choices
            trimLength = 20;
            for l = 1:trimLength*2+1
                ss = l;
                es = length(whichChoices) - (trimLength*2-l+1);
                shifts(l,:) = whichChoices(ss:es);
            end
            
            %generate linear-shifted stim one-hot-encoding matrix
            for l = 1:trimLength*2+1
                ss = l;
                es = length(whichStimOneHot) - (trimLength*2-l+1);
                sohShifts(l,:,:) = whichStimOneHot(ss:es,:);
            end
            
            pseudoChoices = shifts([1:trimLength,trimLength+2:trimLength*2+1],:);
            pseudoStimOneHot = sohShifts([1:trimLength,trimLength+2:trimLength*2+1],:,:);
            testTrials = whichTrials(trimLength+1:3:end-trimLength);
            trainTrials = whichTrials(trimLength+1:end-trimLength);
            trainTrials(ismember(whichTrials(trimLength+1:end-trimLength),testTrials)) = [];

            truncIdx = trimLength+1:nt-trimLength;
            ptestTrials = 1:3:length(pseudoChoices);
            ptrainTrials = 1:length(pseudoChoices);
            ptrainTrials(ptestTrials) = [];
            
            X(t).train = squeeze(neuralData.eta.alignedResps{ETA}(trainTrials,timeRange(t),:));
            X(t).test = squeeze(neuralData.eta.alignedResps{ETA}(testTrials,timeRange(t),:));
            
            OneHot(t).true.train = stimOneHot(trainTrials,:);
            OneHot(t).true.test = stimOneHot(testTrials,:);
            OneHot(t).pseudo.train = pseudoStimOneHot(:,ptrainTrials,:);
            OneHot(t).pseudo.test = pseudoStimOneHot(:,ptestTrials,:);
            
            Y(t).true.train = trueChoices(trainTrials);
            Y(t).true.test = trueChoices(testTrials);
            Y(t).pseudo.train = pseudoChoices(:,ptrainTrials);
            Y(t).pseudo.test = pseudoChoices(:,ptestTrials);
            
            family = 'binomial';
            
        case 'stimIndChoice'
            [~, whichTrials] = selectCondition(expInfo, contrasts, behavioralData, ...
                    initTrialConditions('repeatType','random','movementTime','late','specificRTs',[.8 3]));
            whichStim = trueSide(whichTrials);
            whichChoiceOneHot = choiceOneHot(whichTrials,:);
            %generate linear-shifted stim
            trimLength = 20;
            for l = 1:trimLength*2+1
                ss = l;
                es = length(whichStim) - (trimLength*2-l+1);
                shifts(l,:) = whichStim(ss:es);
            end
            
            %generate linear-shifted choice one-hot-encoding matrix
            for l = 1:trimLength*2+1
                ss = l;
                es = length(whichChoiceOneHot) - (trimLength*2-l+1);
                cohShifts(l,:,:) = whichChoiceOneHot(ss:es,:);
            end
            
            pseudoSide = shifts([1:trimLength,trimLength+2:trimLength*2+1],:);
            pseudoChoiceOneHot = cohShifts([1:trimLength,trimLength+2:trimLength*2+1],:,:);
            testTrials = whichTrials(trimLength+1:3:end-trimLength);
            trainTrials = whichTrials(trimLength+1:end-trimLength);
            trainTrials(ismember(whichTrials(trimLength+1:end-trimLength),testTrials)) = [];

            truncIdx = trimLength+1:nt-trimLength;
            ptestTrials = 1:3:length(pseudoSide);
            ptrainTrials = 1:length(pseudoSide);
            ptrainTrials(ptestTrials) = [];
            
            X(t).train = squeeze(neuralData.eta.alignedResps{ETA}(trainTrials,timeRange(t),:));
            X(t).test = squeeze(neuralData.eta.alignedResps{ETA}(testTrials,timeRange(t),:));
            
            OneHot(t).true.train = choiceOneHot(trainTrials,:);
            OneHot(t).true.test = choiceOneHot(testTrials,:);
            OneHot(t).pseudo.train = pseudoChoiceOneHot(:,ptrainTrials,:);
            OneHot(t).pseudo.test = pseudoChoiceOneHot(:,ptestTrials,:);
            
            Y(t).true.train = trueSide(trainTrials);
            Y(t).true.test = trueSide(testTrials);
            Y(t).pseudo.train = pseudoSide(:,ptrainTrials);
            Y(t).pseudo.test = pseudoSide(:,ptestTrials);
            
            family = 'binomial';
            
        case 'allFeedback'
            [~, whichTrials] = selectCondition(expInfo, getUniqueContrasts(expInfo), behavioralData, ...
                    initTrialConditions('repeatType','random','movementTime','late','specificRTs',[.8 3]));
            whichFeedback = allFeedback(whichTrials);
            %generate linear-shifted choices
            trimLength = 20;
            for l = 1:trimLength*2+1
                ss = l;
                es = length(whichFeedback) - (trimLength*2-l+1);
                shifts(l,:) = whichFeedback(ss:es);
            end
            pseudoFeedback = shifts([1:trimLength,trimLength+2:trimLength*2+1],:);
            testTrials = whichTrials(trimLength+1:3:end-trimLength);
            trainTrials = whichTrials(trimLength+1:end-trimLength);
            trainTrials(ismember(whichTrials(trimLength+1:end-trimLength),testTrials)) = [];
            trimLength = 20;
            truncIdx = trimLength+1:nt-trimLength;
            ptestTrials = 1:3:length(pseudoFeedback);
            ptrainTrials = 1:length(pseudoFeedback);
            ptrainTrials(ptestTrials) = [];
            
            X(t).train = squeeze(neuralData.eta.alignedResps{ETA}(trainTrials,timeRange(t),:));
            X(t).test = squeeze(neuralData.eta.alignedResps{ETA}(testTrials,timeRange(t),:));
            
            Y(t).true.train = allFeedback(trainTrials);
            Y(t).true.test = allFeedback(testTrials);
            Y(t).true.dumb = ones(1,length(testTrials))*mean(allFeedback(trainTrials));
            Y(t).pseudo.train = pseudoFeedback(:,ptrainTrials);
            Y(t).pseudo.test = pseudoFeedback(:,ptestTrials);
            Y(t).pseudo.dumb = ones(1,length(ptestTrials)).*mean(pseudoFeedback(:,ptrainTrials),2);
            family = 'binomial';
            
        case 'velocity'
            [~, whichTrials] = selectCondition(expInfo, getUniqueContrasts(expInfo), behavioralData, ...
                    initTrialConditions('repeatType','random','movementTime','late','specificRTs',[.8 3]));
            testTrials = whichTrials(1:3:end);
            trainTrials = whichTrials;
            trainTrials(ismember(whichTrials,testTrials)) = [];
            trimLength = 20;
            truncIdx = trimLength+1:nt-trimLength;
            
            X(t).train = squeeze(neuralData.eta.alignedResps{ETA}(intersect(truncIdx,trainTrials),timeRange(t),:));
            X(t).test = squeeze(neuralData.eta.alignedResps{ETA}(intersect(truncIdx,testTrials),timeRange(t),:));
            
            Y(t).true.train = trueVel(intersect(truncIdx,trainTrials));
            Y(t).true.test = trueVel(intersect(truncIdx,testTrials));
            Y(t).pseudo.train = pseudoVel(:,intersect(truncIdx,trainTrials)-trimLength);
            Y(t).pseudo.test = pseudoVel(:,intersect(truncIdx,testTrials)-trimLength);
            family = 'gaussian';
        
        case 'RT'
            [~, whichTrials] = selectCondition(expInfo, getUniqueContrasts(expInfo), behavioralData, ...
                    initTrialConditions('repeatType','random','movementTime','all','specificRTs',[.1 3]));
            testTrials = whichTrials(1:3:end);
            trainTrials = whichTrials;
            trainTrials(ismember(whichTrials,testTrials)) = [];
            trimLength = 20;
            truncIdx = trimLength+1:nt-trimLength;
            
            X(t).train = squeeze(neuralData.eta.alignedResps{ETA}(intersect(truncIdx,trainTrials),timeRange(t),:));
            X(t).test = squeeze(neuralData.eta.alignedResps{ETA}(intersect(truncIdx,testTrials),timeRange(t),:));
            
            Y(t).true.train = trueRT(intersect(truncIdx,trainTrials));
            Y(t).true.test = trueRT(intersect(truncIdx,testTrials));
            Y(t).pseudo.train = pseudoRT(:,intersect(truncIdx,trainTrials)-trimLength);
            Y(t).pseudo.test = pseudoRT(:,intersect(truncIdx,testTrials)-trimLength);
            family = 'gaussian';
            
        case 'feedback'
            [~, whichTrials] = selectCondition(expInfo, getUniqueContrasts(expInfo), behavioralData, ...
                    initTrialConditions('repeatType','random','movementTime','late','specificRTs',[.8 3]));
            testTrials = whichTrials(1:3:end);
            trainTrials = whichTrials;
            trainTrials(ismember(whichTrials,testTrials)) = [];
            
            [~,iTrain] = intersect(feedbackIdx,trainTrials);
            [~,iTest] = intersect(feedbackIdx,testTrials);
            
            X(t).train = squeeze(neuralData.eta.alignedResps{ETA}(intersect(feedbackIdx,trainTrials),timeRange(t),:));
            X(t).test = squeeze(neuralData.eta.alignedResps{ETA}(intersect(feedbackIdx,testTrials),timeRange(t),:));
            
            Y(t).true.train = trueFeedback(iTrain);
            Y(t).true.test = trueFeedback(iTest);
            Y(t).pseudo.train = pseudoFeedback(:,iTrain);
            Y(t).pseudo.test = pseudoFeedback(:,iTest);
            family = 'binomial';
            
        case 'block'
            [~, whichTrials] = selectCondition(expInfo, getUniqueContrasts(expInfo), behavioralData, ...
                    initTrialConditions('repeatType','random','movementTime','late','specificRTs',[.8 3]));
            testTrials = whichTrials(1:3:end);
            trainTrials = whichTrials;
            trainTrials(ismember(whichTrials,testTrials)) = [];
            
            X(t).train = squeeze(neuralData.eta.alignedResps{ETA}(trainTrials,timeRange(t),:));
            X(t).test = squeeze(neuralData.eta.alignedResps{ETA}(testTrials,timeRange(t),:));
            
            Y(t).true.train = trueBlocks(trainTrials);
            Y(t).true.test = trueBlocks(testTrials);
            Y(t).true.dumb = ones(1,length(testTrials))*mean(trueBlocks(trainTrials));
            Y(t).pseudo.train = pseudoBlocks(:,trainTrials);
            Y(t).pseudo.test = pseudoBlocks(:,testTrials);
            Y(t).pseudo.dumb = ones(1,length(testTrials)).*mean(pseudoBlocks(:,trainTrials),2);
            family = 'binomial';
            
        case 'value'
            [~, whichTrials] = selectCondition(expInfo, getUniqueContrasts(expInfo), behavioralData, ...
                    initTrialConditions('repeatType','random','movementTime','late','specificRTs',[.8 3]));
            testTrials = whichTrials(1:3:end);
            trainTrials = whichTrials;
            trainTrials(ismember(whichTrials,testTrials)) = [];
            
            X(t).train = squeeze(neuralData.eta.alignedResps{ETA}(trainTrials,timeRange(t),:));
            X(t).test = squeeze(neuralData.eta.alignedResps{ETA}(testTrials,timeRange(t),:));
            
            Y(t).true.train = trueValue(trainTrials);
            Y(t).true.test = trueValue(testTrials);
            Y(t).true.dumb = ones(1,length(testTrials))*mean(trueValue(trainTrials));
            Y(t).pseudo.train = pseudoValue(:,trainTrials);
            Y(t).pseudo.test = pseudoValue(:,testTrials);
            Y(t).pseudo.dumb = ones(1,length(testTrials)).*mean(pseudoValue(:,trainTrials),2);
            family = 'binomial';
    end
end

%% fit options

options.alpha = 1;
options.nlambda = 20;
options.standardize = 'false';
% options.exclude;
nf = 3;

try
    parpool();
catch
end

% %% DEV: FIT TRUE
% 
% if strcmp(Ytype,'choiceIndStim') || strcmp(Ytype,'stimIndChoice') %append the one-hot to the neural matrix
%     prog = 0;
%     fprintf(1,'fitting true: %3d%%\n',prog);
%     %neurons are always subject to a penalty (ones), but the one-hot matrix is never
%     %subject to a penalty (zeros)
%     options.penalty_factor = [ones(1,size(X(t).train,2)) zeros(1,size(OneHot(t).true.train,2))];
%     for t = 1:length(timeRange)
%         training_set = [X(t).train OneHot(t).true.train];
%         test_set = [X(t).test OneHot(t).true.test];
%         
%         %for the smart model, we don't exclude any predictors - we want
%         %both neural and one-hot columns to contribute
%         options.exclude = [];
%         fit.true{t} = cvglmnet(training_set, Y(t).true.train', family, options, 'deviance',nf,[],true);
%         Y(t).true.hat = cvglmnetPredict(fit.true{t}, test_set,'lambda_min','response')';
%         
%         %for the dumb model, we exclude the neural data and see how we get
%         %with just the one-hot columns
%         options.exclude = [1:size(X(t).train,2)]';
%         fit.dumb{t} = glmnet(training_set, Y(t).true.train', family, options);
%         tmp_fits = glmnetPredict(fit.dumb{t}, test_set,[],'response');
%         Y(t).true.dumb = tmp_fits(:,end)';
%         prog = ( 100*(t/length(timeRange)));
%         fprintf(1,'\b\b\b\b%3.0f%%',prog);
%     end
% else
% end
% 
% 
% 
% if strcmp(Ytype,'choiceIndStim') || strcmp(Ytype,'stimIndChoice') %append the one-hot to the neural matrix
%     prog = 0;
%     fprintf(1,'fitting true: %3d%%\n',prog);
%     for t = 1:length(timeRange)
%         fit.true{t} = cvglmnet([X(t).train OneHot(t).true.train],Y(t).true.train,family,options,'deviance',nf,[],true);
%         fit.dumb{t} = cvglmnet([OneHot(t).true.train],Y(t).true.train,family,options,'deviance',nf,[],true);
%         Y(t).true.hat = cvglmnetPredict(fit.true{t},[X(t).test OneHot(t).true.test],'lambda_min','response')';
%         Y(t).true.dumb = cvglmnetPredict(fit.dumb{t},[OneHot(t).true.test],'lambda_min','response')';
%         prog = ( 100*(t/length(timeRange)));
%         fprintf(1,'\b\b\b\b%3.0f%%',prog);
%     end
%     fprintf('\n');
% else
% end
%  for p = 1:np
%     for t = 1:length(timeRange)
%         fit.pseudo{p,t} = cvglmnet([X(t).train squeeze(OneHot(t).pseudo.train(p,:,:))],Y(t).pseudo.train(p,:),family,options,'deviance',nf,[],true);
%         Y(t).pseudo.hat(p,:) = cvglmnetPredict(fit.pseudo{p,t},[X(t).test squeeze(OneHot(t).pseudo.test(p,:,:))],'lambda_min','response')';
%         Y(t).pseudo.dumb = repmat(Y(t).true.dumb,[np 1]);
%         prog = ( 100*(p/np) );
%         fprintf(1,'\b\b\b\b%3.0f%%',prog);
%     end
% end
%% fit true

% fit
if strcmp(Ytype,'choiceIndStim') || strcmp(Ytype,'stimIndChoice') %append the one-hot to the neural matrix
    prog = 0;
    fprintf(1,'fitting true: %3d%%\n',prog);
    %neurons are always subject to a penalty (ones), but the one-hot matrix is never
    %subject to a penalty (zeros)
    options.penalty_factor = [ones(1,size(X(t).train,2)) zeros(1,size(OneHot(t).true.train,2))];
    for t = 1:length(timeRange)
        training_set = [X(t).train OneHot(t).true.train];
        test_set = [X(t).test OneHot(t).true.test];
        
        %for the smart model, we don't exclude any predictors - we want
        %both neural and one-hot columns to contribute
        options.exclude = [];
        fit.true{t} = cvglmnet(training_set, Y(t).true.train', family, options, 'deviance',nf,[],true);
        Y(t).true.hat = cvglmnetPredict(fit.true{t}, test_set,'lambda_min','response')';
        
        %for the dumb model, we exclude the neural data and see how we get
        %with just the one-hot columns
        options.exclude = [1:size(X(t).train,2)]';
        options.lambda = [0.00000001 0.00000001 0.00000001];
        fit.dumb{t} = cvglmnet(training_set, Y(t).true.train', family, options, 'deviance',nf,[],true);
        tmp_fits = cvglmnetPredict(fit.dumb{t}, test_set,'lambda_min','response')';
        Y(t).true.dumb = tmp_fits(:,end)';
        prog = ( 100*(t/length(timeRange)));
        fprintf(1,'\b\b\b\b%3.0f%%',prog);
    end
else
    prog = 0;
    fprintf(1,'fitting true: %3d%%\n',prog);
    for t = 1:length(timeRange)
        fit.true{t} = cvglmnet(X(t).train,Y(t).true.train,family,options,'deviance',nf,[],true);
    %     fit.true{t} = cvglmnet(X(t).train,Y(t).true.train,'multinomial',options,'deviance',5,[],true);
        Y(t).true.hat = cvglmnetPredict(fit.true{t},X(t).test,'lambda_min','response')';
        prog = ( 100*(t/length(timeRange)));
        fprintf(1,'\b\b\b\b%3.0f%%',prog);
    end
    fprintf('\n');
end

%goodness of fit
%(gof = MSE for continuous (stimulus) or gof = LL for binomial)
gof = nan(1,t);
if strcmp(family,'gaussian') %|| strcmp(family,'binomial')
    for t = 1:length(timeRange)
        gof(t) = sum((Y(t).true.test - (Y(t).true.hat)).^2)/length(Y(t).true.test);
    end
else
    for t = 1:length(timeRange)
        term1 = Y(t).true.test.*log(Y(t).true.hat);
        term2 = (1-Y(t).true.test).*log(1-Y(t).true.hat);
        ll(t) = sum(term1 + term2);
        gof(t) = exp(ll(t)/size(Y(t).true.hat,2));
    end
    
    for t = 1:length(timeRange)
        term1 = Y(t).true.test.*log(Y(t).true.dumb);
        term2 = (1-Y(t).true.test).*log(1-Y(t).true.dumb);
        lld(t) = sum(term1 + term2);
    end
end

info_hat = ll./size(Y(1).true.hat,2)./log(2);
info_dumb = lld./size(Y(1).true.dumb,2)./log(2);
mi = info_hat - info_dumb;

%% fit pseudos

np = size(Y(1).pseudo.train,1);

%fit
if strcmp(Ytype,'choiceIndStim') || strcmp(Ytype,'stimIndChoice')
    prog = 0;
    fprintf(1,'fitting pseudos: %3d%%\n',prog);
    %neurons are always subject to a penalty (ones), but the one-hot matrix is never
    %subject to a penalty (zeros)
    options.penalty_factor = [ones(1,size(X(t).train,2)) zeros(1,size(OneHot(t).true.train,2))];
    for p = 1:np
        for t = 1:length(timeRange)
            try
                training_set = [X(t).train squeeze(OneHot(t).pseudo.train(p,:,:))];
                test_set = [X(t).test squeeze(OneHot(t).pseudo.test(p,:,:))];
            catch
                training_set = [X(t).train squeeze(OneHot(t).pseudo.train(p,:,:))'];
                test_set = [X(t).test squeeze(OneHot(t).pseudo.test(p,:,:))'];
            end
            %for the smart model, we don't exclude any predictors - we want
            %both neural and one-hot columns to contribute
            options.exclude = [];
            fit.pseudo{p,t} = cvglmnet(training_set,Y(t).pseudo.train(p,:),family,options,'deviance',nf,[],true);
            Y(t).pseudo.hat(p,:) = cvglmnetPredict(fit.pseudo{p,t},test_set,'lambda_min','response')';
            
            %for the dumb model, we exclude the neural data and see how we get
            %with just the one-hot columns
            options.exclude = [1:size(X(t).train,2)]';
            Y(t).pseudo.dumb = repmat(Y(t).true.dumb,[np 1]);
            prog = ( 100*(p/np) );
            fprintf(1,'\b\b\b\b%3.0f%%',prog);
        end
    end
    fprintf('\n');
else
    prog = 0;
    fprintf(1,'fitting pseudos: %3d%%\n',prog);
    for p = 1:np
        for t = 1:length(timeRange)
            fit.pseudo{p,t} = cvglmnet(X(t).train,Y(t).pseudo.train(p,:),family,options,'deviance',nf,[],true);
            Y(t).pseudo.hat(p,:) = cvglmnetPredict(fit.pseudo{p,t},X(t).test,'lambda_min','response')';
            prog = ( 100*(p/np) );
            fprintf(1,'\b\b\b\b%3.0f%%',prog);
        end
    end
    fprintf('\n');
end

% goodness of fit
gof_pseudo = nan(np,t);
if strcmp(family,'gaussian')
    for p = 1:np
        for t = 1:length(timeRange)
            gof_pseudo(p,t) = sum((Y(t).pseudo.test(p,:) - (Y(t).pseudo.hat(p,:))).^2)/length(Y(t).pseudo.test(p,:));
        end
    end
else
    for p = 1:np
        for t = 1:length(timeRange)
            term1 = Y(t).pseudo.test(p,:).*log(Y(t).pseudo.hat(p,:));
            term2 = (1-Y(t).pseudo.test(p,:)).*log(1-Y(t).pseudo.hat(p,:));
            ll(p,t) = sum(term1 + term2);
            gof_pseudo(p,t) = exp(ll(p,t)/size(Y(t).pseudo.hat(p,:),2));
        end
        
        for t = 1:length(timeRange)
            term1 = Y(t).pseudo.test(p,:).*log(Y(t).pseudo.dumb(p,:));
            term2 = (1-Y(t).pseudo.test(p,:)).*log(1-Y(t).pseudo.dumb(p,:));
            lld(p,t) = sum(term1 + term2);
        end
    end
end

info_phat = ll./size(Y(1).pseudo.hat,2)./log(2);
info_pdumb = lld./size(Y(1).pseudo.dumb,2)./log(2);
mip = info_phat - info_pdumb;

    
%% plot
% 
% ew = neuralData.eta.eventWindow(timeRange);
% 
% try %true vs pseudos
%     gof_UB = prctile(gof_pseudo,97.5);
%     gof_LB = prctile(gof_pseudo,2.5);
%     expRef = strcat(expInfo.mouseName,'_',expInfo.expDate);
%     maxy = max(max(max([gof_LB; gof; gof_UB])));
%     miny = min(min(min([gof_LB; gof; gof_UB])));
%     buffer = 0.1*abs((maxy-miny));
%     figure;
%     plotSignal(ew,gof,gof_UB,gof_LB,[0 0 0],'-');
%     prettyPlot(gca);
%     line([0 0],[miny-buffer maxy+buffer],'LineStyle','--','Color',[.5 .5 .5]);
%     xlim([-.5 2]);
%     ylim([miny-buffer maxy+buffer]);
%     xlabel('Time from go cue (s)')
%     ylabel('1/error')
%     title(strcat(expRef,{': Neural decoding of '},Ytype),'Interpreter','none')
% catch %just true
%     expRef = strcat(expInfo.mouseName,'_',expInfo.expDate);
%     maxy = max(max(max([gof])));
%     miny = min(min(min([gof])));
%     buffer = 0.1*abs((maxy-miny));
%     figure;
%     plot(neuralData.eta.eventWindow,gof,'Color',[0 0 0],'LineWidth',2,'LineStyle','-');
%     prettyPlot(gca);
%     line([0 0],[0-buffer maxy+buffer],'LineStyle','--','Color',[.5 .5 .5]);
%     xlim([-.5 2]);
%     ylim([0.4 1]);
%     xlabel('Time from stimulus onset (s)')
%     ylabel('ll')
%     title(strcat(expRef,{': Neural decoding of '},Ytype),'Interpreter','none')
% end


%% print and save

% printfig(gcf, char(strcat(expRef,{' neural decoding of '},Ytype)));

% sessDir = fullfile('G:\Workspaces',expInfo.mouseName,expInfo.expDate,num2str(expInfo.expNum));
% cd(sessDir)
% 
% clearvars -except family fit gof gof_UB gof_LB gof_pseudo options testTrials whichTrials trainTrials X Y Ytype
% 
% if strcmp(Ytype,'stimulus')
%     save stimulusDecoderFits.mat -v7.3
% elseif strcmp(Ytype,'choice')
%     save choiceDecoderFits.mat -v7.3
% elseif strcmp(Ytype,'block')
%     save blockDecoderFits.mat -v7.3
% elseif strcmp(Ytype,'feedback')
%     save feedbackDecoderFits.mat -v7.3
% elseif strcmp(Ytype,'value')
%     save valueDecoderFits.mat -v7.3
% end

%% WORKBENCH

% function fharrell
%     % this is from a blog post by Frank Harrell, which proposes the following R2 
%     % as a way to expressed explained outcome variance in probability models
%     % https://www.fharrell.com/post/addvalue
%     for t = 1:size(neuralData.eta.alignedResps{1},2)
%         %variance of the model, (Phat probability that Y = 1)
%         term1 = var(Y(t).true.hat); 
% 
%         %variance of 'the rest' (1-Phat)
%         term2 = sum((Y(t).true.hat).*(1-Y(t).true.hat))/size(Y(t).true.hat,2);
% 
%         %proportion
%         R2(t) = term1/(term1+term2);
%     end
% end
            
            
end            
            
            
            
            
            
