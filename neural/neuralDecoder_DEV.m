function [mi, mip] = neuralDecoder(expInfo, behavioralData, neuralData, Ytype, ETA)

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
            X(t).true = squeeze(neuralData.eta.alignedResps{ETA}(whichTrials,timeRange(t),:));
    
            Y(t).true.smart = trueSide(whichTrials);
            Y(t).true.dumb = ones(1,length(whichTrials))*mean(trueSide(whichTrials));
            Y(t).pseudo.smart = pseudoSide(:,whichTrials);
            Y(t).pseudo.dumb = ones(1,length(whichTrials)).*mean(pseudoSide(:,whichTrials),2);
            family = 'binomial';
            
        case 'stimulus'
            % generally filter trials (e.g., only patient ones)
            [~, whichTrials] = selectCondition(expInfo, getUniqueContrasts(expInfo), behavioralData, ...
                    initTrialConditions('repeatType','random','movementTime','late','specificRTs',[.8 3]));
            X(t).true = squeeze(neuralData.eta.alignedResps{ETA}(whichTrials,timeRange(t),:));
    
            Y(t).true.smart = trueStimuli(whichTrials);
            Y(t).true.dumb = ones(1,length(testTrials))*mean(trueStimuli(trainTrials));
            Y(t).pseudo.smart = pseudoStimuli(:,whichTrials);
            Y(t).pseudo.dumb = ones(1,length(whichTrials)).*mean(pseudoStimuli(:,whichTrials),2);
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
            trimLength = 20;
            truncIdx = trimLength+1:length(whichChoices)-trimLength;
            
            X(t).true = squeeze(neuralData.eta.alignedResps{ETA}(whichTrials(truncIdx),timeRange(t),:));
            
            Y(t).true.smart = trueChoices(whichTrials(truncIdx));
            Y(t).true.dumb = ones(1,length(truncIdx))*mean(trueChoices(whichTrials(truncIdx)));
            Y(t).pseudo.smart = pseudoChoices;
            Y(t).pseudo.dumb = ones(1,length(truncIdx)).*mean(pseudoChoices,2);
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
            trimLength = 20;
            truncIdx = trimLength+1:length(whichChoices)-trimLength;
            
            X(t).true = squeeze(neuralData.eta.alignedResps{ETA}(whichTrials(truncIdx),timeRange(t),:));
            
            Y(t).true.smart = trueChoices(whichTrials(truncIdx));
            Y(t).true.dumb = ones(1,length(truncIdx))*mean(trueChoices(whichTrials(truncIdx)));
            Y(t).pseudo.smart = pseudoChoices;
            Y(t).pseudo.dumb = ones(1,length(truncIdx)).*mean(pseudoChoices,2);
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
            

            truncIdx = trimLength+1:length(whichTrials)-trimLength;
            
            X(t).true = squeeze(neuralData.eta.alignedResps{ETA}(whichTrials(truncIdx),timeRange(t),:));
            
            OneHot(t).true = stimOneHot(whichTrials(truncIdx),:);
            OneHot(t).pseudo = pseudoStimOneHot;
            
            Y(t).true.smart = trueChoices(whichTrials(truncIdx));
            Y(t).pseudo.smart = pseudoChoices;            
            family = 'binomial';
            
        case 'stimIndChoice'
            [~, whichTrials] = selectCondition(expInfo, [0], behavioralData, ...
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

            truncIdx = trimLength+1:length(whichTrials)-trimLength;
            
            X(t).true = squeeze(neuralData.eta.alignedResps{ETA}(whichTrials(truncIdx),timeRange(t),:));
            
            OneHot(t).true = choiceOneHot(whichTrials(truncIdx),:);
            OneHot(t).pseudo = pseudoChoiceOneHot;
            
            Y(t).true.smart = trueSide(whichTrials(truncIdx));
            Y(t).pseudo.smart = pseudoSide;
            
            family = 'binomial';
            
        case 'allFeedback'
            [~, whichTrials] = selectCondition(expInfo, getUniqueContrasts(expInfo), behavioralData, ...
                    initTrialConditions('repeatType','random','movementTime','late','specificRTs',[.8 3]));
            whichFeedback = allFeedback(whichTrials);
            %generate linear-shifted trials
            trimLength = 20;
            for l = 1:trimLength*2+1
                ss = l;
                es = length(whichFeedback) - (trimLength*2-l+1);
                shifts(l,:) = whichFeedback(ss:es);
            end
            pseudoFeedback = shifts([1:trimLength,trimLength+2:trimLength*2+1],:);
            
            truncIdx = trimLength+1:length(whichTrials)-trimLength;
            
            X(t).true = squeeze(neuralData.eta.alignedResps{ETA}(whichTrials(truncIdx),timeRange(t),:));
            
            Y(t).true.smart = allFeedback(whichTrials(truncIdx));
            Y(t).true.dumb = ones(1,length(truncIdx))*mean(allFeedback(whichTrials(truncIdx)));
            Y(t).pseudo.smart = pseudoFeedback;
            Y(t).pseudo.dumb = ones(1,length(truncIdx)).*mean(pseudoFeedback,2);
            family = 'binomial';
            
        case 'velocity'
            [~, whichTrials] = selectCondition(expInfo, getUniqueContrasts(expInfo), behavioralData, ...
                    initTrialConditions('repeatType','random','movementTime','late','specificRTs',[.8 3]));
            whichVel = trueVel(whichTrials);
            %generate linear-shifted trials
            trimLength = 20;
            for l = 1:trimLength*2+1
                ss = l;
                es = length(whichVel) - (trimLength*2-l+1);
                shifts(l,:) = whichVel(ss:es);
            end
            truncIdx = trimLength+1:length(whichTrials)-trimLength;
            pseudoVel = shifts([1:trimLength,trimLength+2:trimLength*2+1],:);
            
            X(t).true = squeeze(neuralData.eta.alignedResps{ETA}(whichTrials(truncIdx),timeRange(t),:));
            
            Y(t).true.smart = trueVel(whichTrials(truncIdx));
            Y(t).true.dumb = ones(1,length(truncIdx))*mean(trueVel(whichTrials(truncIdx)));
            Y(t).pseudo.smart = pseudoVel;
            Y(t).pseudo.dumb = ones(1,length(truncIdx)).*mean(pseudoVel,2);
            family = 'gaussian';
        
        case 'RT'
            [~, whichTrials] = selectCondition(expInfo, getUniqueContrasts(expInfo), behavioralData, ...
                    initTrialConditions('repeatType','random','movementTime','all','specificRTs',[.1 3]));
            whichRT = trueRT(whichTrials);
            %generate linear-shifted trials
            trimLength = 20;
            for l = 1:trimLength*2+1
                ss = l;
                es = length(whichRT) - (trimLength*2-l+1);
                shifts(l,:) = whichRT(ss:es);
            end
            truncIdx = trimLength+1:length(whichTrials)-trimLength;
            pseudoRT = shifts([1:trimLength,trimLength+2:trimLength*2+1],:);
            
            X(t).true = squeeze(neuralData.eta.alignedResps{ETA}(whichTrials(truncIdx),timeRange(t),:));
            
            Y(t).true.smart = trueRT(whichTrials(truncIdx));
            Y(t).true.dumb = ones(1,length(truncIdx))*mean(trueRT(whichTrials(truncIdx)));
            Y(t).pseudo.smart = pseudoRT;
            Y(t).pseudo.dumb = ones(1,length(truncIdx)).*mean(pseudoRT,2);
            family = 'gaussian';
            
        case 'block'
            [~, whichTrials] = selectCondition(expInfo, getUniqueContrasts(expInfo), behavioralData, ...
                    initTrialConditions('repeatType','random','movementTime','late','specificRTs',[.8 3]));
            
            X(t).true = squeeze(neuralData.eta.alignedResps{ETA}(whichTrials,timeRange(t),:));
            
            Y(t).true.smart = trueBlocks(whichTrials);
            Y(t).true.dumb = ones(1,length(whichTrials))*mean(trueBlocks(whichTrials));
            Y(t).pseudo.smart = pseudoBlocks(:,whichTrials);
            Y(t).pseudo.dumb = ones(1,length(whichTrials)).*mean(pseudoBlocks(:,whichTrials),2);
            family = 'binomial';
            
        case 'value'
            [~, whichTrials] = selectCondition(expInfo, getUniqueContrasts(expInfo), behavioralData, ...
                    initTrialConditions('repeatType','random','movementTime','late','specificRTs',[.8 3]));
            testTrials = whichTrials(1:3:end);
            trainTrials = whichTrials;
            trainTrials(ismember(whichTrials,testTrials)) = [];
            
            X(t).true = squeeze(neuralData.eta.alignedResps{ETA}(whichTrials,timeRange(t),:));
            
            Y(t).true.smart = trueValue(whichTrials);
            Y(t).true.dumb = ones(1,length(whichTrials))*mean(trueValue(whichTrials));
            Y(t).pseudo.smart = pseudoValue(:,whichTrials);
            Y(t).pseudo.dumb = ones(1,length(whichTrials)).*mean(pseudoValue(:,whichTrials),2);
            family = 'binomial';
    end
end

%% fit options

options.alpha = 1;
options.nlambda = 20;
options.standardize = 'false';
nf = 3;

try
    parpool();
catch
end

%% fit true

% fit
if strcmp(Ytype,'choiceIndStim') || strcmp(Ytype,'stimIndChoice') %append the one-hot to the neural matrix
    clear options
    options.alpha = 1;
    options.nlambda = 20;
    options.standardize = 'false';
    nf = 3;
    
    prog = 0;
    fprintf(1,'fitting true: %3d%%\n',prog);

    %neurons are always subject to a penalty (ones), but the one-hot matrix is never
    %subject to a penalty (zeros)
    options.penalty_factor = [ones(1,size(X(t).true,2)) zeros(1,size(OneHot(t).true,2))];
    
    %for the smart model, we don't exclude any predictors - we want
    %both neural and one-hot columns to contribute
    options.exclude = [];
    for t = 1:length(timeRange)
        preds = [X(t).true OneHot(t).true];
        fit.true{t} = cvglmnet(preds, Y(t).true.smart', family, options, 'deviance',nf,[],true);
        info_hat(t) = fit.true{t}.cvm(fit.true{t}.lambda == fit.true{t}.lambda_min)  / (-2 * log(2));
    end
    
    %for the dumb model, we exclude the neural data and see how we get
    %with just the one-hot columns. and generate some fake lambdas since
    %they're not actually being used (this is to avoid a glmnet bug)
    options.exclude = [1:size(X(t).true,2)]';
    options.lambda = ones(1,options.nlambda) * 0.0000001;
    for t = 1:length(timeRange)
        preds = [X(t).true OneHot(t).true];
        fit.dumb{t} = cvglmnet(preds, Y(t).true.smart', family, options, 'deviance',nf,[],true);
        info_dumb(t) = fit.dumb{t}.cvm(1) / (-2 * log(2));
        prog = ( 100*(t/length(timeRange)));
        fprintf(1,'\b\b\b\b%3.0f%%',prog);
    end

else
    clear options
    options.alpha = 1;
    options.nlambda = 20;
    options.standardize = 'false';
    nf = 3;
    
    prog = 0;
    fprintf(1,'fitting true: %3d%%\n',prog);
    for t = 1:length(timeRange)
        fit.true{t} = cvglmnet(X(t).true,Y(t).true.smart,family,options,'deviance',nf,[],true);
        Y(t).true.hat = cvglmnetPredict(fit.true{t},X(t).true,'lambda_min','response')';
        prog = ( 100*(t/length(timeRange)));
        fprintf(1,'\b\b\b\b%3.0f%%',prog);
    end
    fprintf('\n');
    
    for t = 1:length(timeRange)
        term1 = Y(t).true.smart.*log(Y(t).true.hat);
        term2 = (1-Y(t).true.smart).*log(1-Y(t).true.hat);
        lls(t) = sum(term1 + term2);
        gof(t) = exp(lls(t)/size(Y(t).true.hat,2));
    end
    
    for t = 1:length(timeRange)
        term1 = Y(t).true.smart.*log(Y(t).true.dumb);
        term2 = (1-Y(t).true.smart).*log(1-Y(t).true.dumb);
        lld(t) = sum(term1 + term2);
    end
    info_hat = lls./size(Y(1).true.hat,2)./log(2);
    info_dumb = lld./size(Y(1).true.dumb,2)./log(2);
end

mi = info_hat - info_dumb;

%% fit pseudos

np = size(Y(1).pseudo.smart,1);

%fit
if strcmp(Ytype,'choiceIndStim') || strcmp(Ytype,'stimIndChoice')
    clear options
    options.alpha = 1;
    options.nlambda = 20;
    options.standardize = 'false';
    nf = 3;
    
    prog = 0;
    fprintf(1,'fitting pseudos: %3d%%\n',prog);
    
    %neurons are always subject to a penalty (ones), but the one-hot matrix is never
    %subject to a penalty (zeros)
    options.penalty_factor = [ones(1,size(X(t).true,2)) zeros(1,size(OneHot(t).pseudo,2))];
    for p = 1:np
        
        %for the smart model, we don't exclude any predictors - we want
        %both neural and one-hot columns to contribute
        options.exclude = [];
        for t = 1:length(timeRange)
            try
                preds = [X(t).true squeeze(OneHot(t).pseudo(p,:,:))];
            catch
                preds = [X(t).true squeeze(OneHot(t).pseudo(p,:,:))'];
            end
            fit.pseudo{p,t} = cvglmnet(preds,Y(t).pseudo.smart(p,:),family,options,'deviance',nf,[],true);
            info_phat(p,t) = fit.pseudo{p,t}.cvm(fit.pseudo{p,t}.lambda == fit.pseudo{p,t}.lambda_min)  / (-2 * log(2));    
        end
        prog = ( 100*(t/length(timeRange)));
        fprintf(1,'\b\b\b\b%3.0f%%',prog);
    end
    
    %for the dumb model, we exclude the neural data and see how we get
    %with just the one-hot columns. this is identical to the true
    %fitting since the one-hot and Y outcome don't shift against each
    %other
    info_pdumb = repmat(info_dumb,[np 1]);
    mip = info_phat-info_pdumb;
    fprintf('\n');
else
    
    clear options
    options.alpha = 1;
    options.nlambda = 20;
    options.standardize = 'false';
    nf = 3;
    
    prog = 0;
    fprintf(1,'fitting pseudos: %3d%%\n',prog);
    for p = 1:np
        for t = 1:length(timeRange)
            fit.pseudo{p,t} = cvglmnet(X(t).true,Y(t).pseudo.smart(p,:),family,options,'deviance',nf,[],true);
            Y(t).pseudo.hat(p,:) = cvglmnetPredict(fit.pseudo{p,t},X(t).true,'lambda_min','response')';
        end
        prog = ( 100*(p/np) );
        fprintf(1,'\b\b\b\b%3.0f%%',prog);
    end
    
    for p = 1:np
        for t = 1:length(timeRange)
            term1 = Y(t).pseudo.smart(p,:).*log(Y(t).pseudo.hat(p,:));
            term2 = (1-Y(t).pseudo.smart(p,:)).*log(1-Y(t).pseudo.hat(p,:));
            lls(p,t) = sum(term1 + term2);
        end
        
        for t = 1:length(timeRange)
            term1 = Y(t).pseudo.smart(p,:).*log(Y(t).pseudo.dumb(p,:));
            term2 = (1-Y(t).pseudo.smart(p,:)).*log(1-Y(t).pseudo.dumb(p,:));
            lld(p,t) = sum(term1 + term2);
        end
    end
    
    info_phat = lls./size(Y(1).pseudo.hat,2)./log(2);
    info_pdumb = lld./size(Y(1).pseudo.dumb,2)./log(2);
    mip = info_phat - info_pdumb;
    fprintf('\n');
end            
            
end            
            
            
            
            
            
