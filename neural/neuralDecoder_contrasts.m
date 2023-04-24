function [fits, mutual_info, allX, allY, allD] = neuralDecoder_contrasts(expInfo, behavioralData, neuralData, Ytype, ETA)

%% Predict a task variable from neural activity
% This takes a 2D vector of neural activity at a particular trial timepoint (X) and uses
% this in lasso logistic regression to try and predict a 1D Y outcome
% variable. A model is fitted for every timepoint (n = 16)

% try
%     parpool();
% catch
% end

%% Collect task data for Y vectors containing true data
nt = length(behavioralData.eventTimes(1).daqTime);

% extract stimulus values
trueStimuli = expInfo.block.events.contrastValues(1:nt);
trialCorrectChoice = expInfo.block.events.correctResponseValues(1:nt);
trueChoices = expInfo.block.events.responseValues(1:nt);

% assign the 0% stimuli as either 'left' or 'right' depending on the
% preassigned correct choice (not the mouse's choice)
trueStimuli(trueStimuli == 0) = eps;
trueStimuli(abs(trueStimuli) < .05) = ...
    trueStimuli(abs(trueStimuli) < .05).* trialCorrectChoice(abs(trueStimuli) < .05);
trueSide = trueStimuli > 0;

%% Generate Y vectors containing pseudo data
%shuffle method for independent vars like stim
%pseudo-generated from task code for block and value
%these go into Y.pseudo, below
np = 100;

%generate randomly shuffled stims
pseudoStimuli = nan(np,nt);
for p = 1:np
    randIndices(p,:) = randsample([1:length(trueStimuli)],length(trueStimuli));
    pseudoStimuli(p,:) = trueStimuli(randIndices(p,:));
    pseudoSide(p,:) = trueSide(randIndices(p,:));
end

%% Convert to contra/ipsi

if expInfo.hemisphere > 0
    trueStimuli = -trueStimuli;
    trueChoices = (trueChoices == -1);
    pseudoStimuli = -pseudoStimuli;
    trueSide = trueSide == 0;
    pseudoSide = pseudoSide == 0;
else
    trueChoices = (trueChoices == 1);
end

%% Generate one-hot encoders
%this becomes the design matrix D, below

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


%% %% Generate other Y vectors containing pseudo data
%use linear-shift for non-independent vars with unknown distribution
%these go into Y.pseudo, below

[~, whichTrials] = selectCondition(expInfo, getUniqueContrasts(expInfo), behavioralData, ...
    initTrialConditions('repeatType','random','movementTime','late','specificRTs',[.8 3]));
trimLength = 20;
truncIdx = trimLength+1:length(whichTrials)-trimLength;

for l = 1:trimLength*2+1
    ss = l;
    es = length(whichTrials) - (trimLength*2-l+1);
    shiftChoices(l,:) = trueChoices(whichTrials(ss:es));
    shiftSide(l,:) = trueSide(whichTrials(ss:es));
    sohShifts(l,:,:) = stimOneHot(whichTrials(ss:es),:);
    cohShifts(l,:,:) = choiceOneHot(whichTrials(ss:es),:);
end

shiftChoices(trimLength+1,:) = [];
shiftSide(trimLength+1,:) = [];
sohShifts(trimLength+1,:,:) = [];
cohShifts(trimLength+1,:,:) = [];


%% Set up X, Y, and D arrays for fitting

%subsample time range to save time
st = -1.5;
et = 3;
[~,si] = min(abs(neuralData.eta.eventWindow - st));
[~,ei] = min(abs(neuralData.eta.eventWindow - et));
timeRange = si:2:ei;  
% timeRange = 6:2:36;

contrasts = unique(abs(getUniqueContrasts(expInfo)));
% contrasts = {absContrasts(1) absContrasts(2:3), absContrasts(4:5)};

for c = 1:length(contrasts)
    clear X Y fit
    fprintf('Loading %d...',contrasts(c))
    

    allTrials = whichTrials(truncIdx);
    [~, testTrials] = selectCondition(expInfo, [-contrasts(c) contrasts(c)], behavioralData, ...
            initTrialConditions('repeatType','random','movementTime','late','specificRTs',[.8 3]));            
    trainTrials = allTrials;
    trainTrials(ismember(allTrials,testTrials)) = [];
    testTrials = intersect(allTrials,testTrials);
    [~,trainIdx] = intersect(allTrials,trainTrials);
    [~,testIdx] = intersect(allTrials,testTrials);

    %collect neural activity into the X matrix
    % X: neural predictor matrix (trials x neurons x timepoints)
    X.test = nan(length(testTrials),size(neuralData.eta.alignedResps{ETA},3),length(timeRange));
    X.train = nan(length(trainTrials),size(neuralData.eta.alignedResps{ETA},3),length(timeRange));
    for t = 1:length(timeRange)
        X.train(:,:,t) = squeeze(neuralData.eta.alignedResps{ETA}(trainTrials,timeRange(t),:));
        X.test(:,:,t) = squeeze(neuralData.eta.alignedResps{ETA}(testTrials,timeRange(t),:));
    end
    
    switch Ytype
        case 'side_contrast'
            %collect values for the Y matrices
            Y.true.train = trueSide(trainTrials);
            Y.true.test = trueSide(testTrials);
            Y.pseudo.train = shiftSide(:,trainIdx);
            Y.pseudo.test = shiftSide(:,testIdx);

            D.true.train = choiceOneHot(trainTrials,:);
            D.true.test = choiceOneHot(testTrials,:);
            D.pseudo.train = cohShifts(:,trainIdx);
            D.pseudo.test = cohShifts(:,testIdx);
            family = 'binomial';
        case  'choice_contrast'
            %collect values for the Y matrices
            Y.true.train = trueChoices(trainTrials);
            Y.true.test = trueChoices(testTrials);
            Y.pseudo.train = shiftChoices(:,trainIdx);
            Y.pseudo.test = shiftChoices(:,testIdx);

            D.true.train = stimOneHot(trainTrials,:);
            D.true.test = stimOneHot(testTrials,:);
            D.pseudo.train = sohShifts(:,trainIdx,:);
            D.pseudo.test = sohShifts(:,testIdx,:);
            family = 'binomial';
    end
    
    %% Predict the true data: full model vs naive model

    %For these two Ytypes, we compare a full model that uses X (neural
    %data) and D (one-hot encoder holding trial stimulus or choice
    %identity), against a naive model that ONLY uses D. 

    %Both models are trying to predict true Y data
    clear options
    options.alpha = 1;
    options.nlambda = 20;
    options.standardize = 'false';
    nf = 3;

    %X is always subject to a penalty (ones), but D is never
    %subject to a penalty (zeros)
    options.penalty_factor = [ones(1,size(X.train,2)) zeros(1,size(D.true.train,2))];

    %For the full model, we don't exclude any predictors. We want to
    %use both X and D
    options.exclude = [];
    prog = 0;
    fprintf(1,'fitting true (full model): %3d%%\n',prog);
    for t = 1:length(timeRange)
        preds.train = [X.train(:,:,t) D.true.train];
        preds.test = [X.test(:,:,t) D.true.test];
        fit.full{t} = cvglmnet(preds.train, Y.true.train', family, options, 'deviance',nf,[],false);
        Y.full(t,:) = cvglmnetPredict(fit.full{t},preds.test,'lambda_min','response')';
        prog = ( 100*(t/length(timeRange)));
        fprintf(1,'\b\b\b\b%3.0f%%',prog);
    end
    fprintf(1,'\n')
    
    %For the naive model, we exclude X and see how we do with D alone
    options.exclude = [1:size(X.test,2)]';
    options.lambda = linspace(0,1,options.nlambda); %generate some meaningless lambdas to avoid glmnet bug
    prog = 0;
    fprintf(1,'fitting true (naive model): %3d%%\n',prog);
    
    preds.train = [X.train(:,:,1) D.true.train];
    preds.test = [X.test(:,:,1) D.true.test];
    fit.naive{1} = cvglmnet(preds.train, Y.true.train', family, options, 'deviance',nf,[],false);
    yy = cvglmnetPredict(fit.naive{1},preds.test,'lambda_min','response')';
    Y.naive = repmat(yy, [length(timeRange) 1]);
    prog = ( 100*(t/length(timeRange)));
    fprintf(1,'\b\b\b\b%3.0f%%',prog);
    fprintf('\n');
    
    %compute LL for true vs full model
    for t = 1:length(timeRange)
        term1 = Y.true.test.*log(Y.full(t,:));
        term2 = (1-Y.true.test).*log(1-Y.full(t,:));
        ll_full(t) = sum(term1 + term2);
    end

    %we generated the naive model outcomes earlier
    %so compute LL for true vs naive model
    for t = 1:length(timeRange)
        term1 = Y.true.test.*log(Y.naive(t,:));
        term2 = (1-Y.true.test).*log(1-Y.naive(t,:));
        ll_naive(t) = sum(term1 + term2);
    end
    info_full = ll_full./size(Y.full,2)./log(2);
    info_naive = ll_naive./size(Y.naive,2)./log(2);

    %% Predict the pseudo data: full model vs naive model

    np = size(Y.pseudo.test,1);

    %For these two Ytypes, we compare a full model that uses X (neural
    %data) and D (one-hot encoder holding trial stimulus or choice
    %identity), against a naive model that ONLY uses D. 

    %Both models are trying to predict true Y data
    clear options
    options.alpha = 1;
    options.nlambda = 20;
    options.standardize = 'false';
    nf = 3;

    %X is always subject to a penalty (ones), but D is never
    %subject to a penalty (zeros)
    options.penalty_factor = [ones(1,size(X.train,2)) zeros(1,size(D.true.train,2))];

    %For the full model, we don't exclude any predictors. We want to
    %use both X and D
    options.exclude = [];
    
    prog = 0;
    fprintf(1,'fitting pseudo (full model): %3d%%\n',prog);
    
    for p = 1:np
        for t = 1:length(timeRange)
            switch Ytype
                case 'side_contrast'
                    preds.train = [X.train(:,:,t) D.pseudo.train(p,:)'];
                    preds.test = [X.test(:,:,t) D.pseudo.test(p,:)'];
                case 'choice_contrast'
                    preds.train = [X.train(:,:,t) squeeze(D.pseudo.train(p,:,:))];
                    preds.test = [X.test(:,:,t) squeeze(D.pseudo.test(p,:,:))];
            end
            fit.pseudo.full{p,t} = cvglmnet(preds.train,Y.pseudo.train(p,:),family,options,'deviance',nf,[],false);
            Y.pseudo_full(t,p,:) = cvglmnetPredict(fit.pseudo.full{p,t},preds.test,'lambda_min','response')';
        end
        prog = (100*(p/np));
        fprintf(1,'\b\b\b\b%3.0f%%',prog);
    end
    fprintf('\n');
    
    
    %For the naive model, we exclude X and see how we do with D alone
    options.exclude = [1:size(X.test,2)]';
    options.lambda = linspace(0,1,options.nlambda); %generate some meaningless lambdas to avoid glmnet bug
    prog = 0;
    fprintf(1,'fitting pseudo (naive model): %3d%%\n',prog);
    
    for p = 1:np
        switch Ytype
            case 'side_contrast'
                preds.train = [X.train(:,:,1) D.pseudo.train(p,:)'];
                preds.test = [X.test(:,:,1) D.pseudo.test(p,:)'];
            case 'choice_contrast'
                preds.train = [X.train(:,:,1) squeeze(D.pseudo.train(p,:,:))];
                preds.test = [X.test(:,:,1) squeeze(D.pseudo.test(p,:,:))];
        end
        fit.pseudo.naive{p} = cvglmnet(preds.train,Y.pseudo.train(p,:),family,options,'deviance',nf,[],false);
        yy = cvglmnetPredict(fit.pseudo.naive{p},preds.test,'lambda_min','response')';
        Y.pseudo_naive(:,p,:) = repmat(yy,[length(timeRange), 1]);
        
        prog = (100*(p/np));
        fprintf(1,'\b\b\b\b%3.0f%%',prog);
    end
    fprintf('\n');
    
    for p = 1:np
        %compute LL for pseudo vs full model
        for t = 1:length(timeRange)
            term1 = Y.pseudo.test(p,:).*log(squeeze(Y.pseudo_full(t,p,:))');
            term2 = (1-Y.pseudo.test(p,:)).*log(1-squeeze(Y.pseudo_full(t,p,:))');
            llf(p,t) = sum(term1 + term2);
            
            term1 = Y.pseudo.test(p,:).*log(squeeze(Y.pseudo_naive(t,p,:))');
            term2 = (1-Y.pseudo.test(p,:)).*log(1-squeeze(Y.pseudo_naive(t,p,:))');
            lln(p,t) = sum(term1 + term2);
        end
    end

    info_pfull = llf./size(Y.pseudo_full,3)./log(2);
    info_pnaive = lln./size(Y.pseudo_naive,3)./log(2);
    
    
    
    
    %% save
    
    allX.(matlab.lang.makeValidName(strcat('con',num2str(contrasts(c))))) = X;
    allY.(matlab.lang.makeValidName(strcat('con',num2str(contrasts(c))))) = Y;
    allD.(matlab.lang.makeValidName(strcat('con',num2str(contrasts(c))))) = D;
    
    fits.true.full.(matlab.lang.makeValidName(strcat('con',num2str(contrasts(c))))) = fit.full;
    fits.true.naive.(matlab.lang.makeValidName(strcat('con',num2str(contrasts(c))))) = fit.naive;
    fits.pseudo.full.(matlab.lang.makeValidName(strcat('con',num2str(contrasts(c))))) = fit.pseudo.full;
    fits.pseudo.naive.(matlab.lang.makeValidName(strcat('con',num2str(contrasts(c))))) = fit.pseudo.naive;
    
    %the mutual information informs how well the full model does over the naive model
    mutual_info.true.(matlab.lang.makeValidName(strcat('con',num2str(contrasts(c))))) = info_full - info_naive;
    mutual_info.pseudo.(matlab.lang.makeValidName(strcat('con',num2str(contrasts(c))))) = info_pfull - info_pnaive;
end%%







