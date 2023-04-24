function [fits, mutual_info, allX, allY, allD] = neuralDecoder(expInfo, behavioralData, neuralData, Ytype, ETA)

%% Predict a task variable from neural activity
% This takes a 2D vector of neural activity at a particular trial timepoint (X) and uses
% this in lasso logistic regression to try and predict a 1D Y outcome
% variable. A model is fitted for every timepoint (n = 16)

try
    parpool();
catch
end

%% Collect task data for Y vectors containing true data
nt = length(behavioralData.eventTimes(1).daqTime);

% extract stimulus, choice, feedback, value, and block values
trueStimuli = expInfo.block.events.contrastValues(1:nt);
trialCorrectChoice = expInfo.block.events.correctResponseValues(1:nt);
trueBlocks = expInfo.block.events.highRewardSideValues(1:nt);
trueChoices = expInfo.block.events.responseValues(1:nt);
trueFeedback = double(expInfo.block.events.feedbackValues(1:nt));

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

% fetch RT, impulsive trials, and maxVel
trueVel = abs(behavioralData.wheelMoves.epochs(5).peakVel(1:nt));
trueRT = behavioralData.wheelMoves.epochs(5).onsetTimes(1:nt) - behavioralData.eventTimes(1).daqTime(1:nt);
[~, impTrials] = selectCondition(expInfo, getUniqueContrasts(expInfo), behavioralData, ...
    initTrialConditions('responseType','all','movementTime','early','specificRTs',[0.1 .8]));     
trueImp = false(1,nt);
trueImp(impTrials) = true;

%% Generate Y vectors containing pseudo data
%shuffle method for independent vars like stim
%pseudo-generated from task code for block and value
%these go into Y.pseudo, below
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

%generate pseudoblocks
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

%generate pseudo high/low trials based on the pseudoblocks + true stim
for p = 1:np
    pseudoValue(p,pseudoBlocks(p,:).*sign(trueStimuli) == -1) = 0;
    pseudoValue(p,pseudoBlocks(p,:).*sign(trueStimuli) == 1) = 1;
end

%% Convert to contra/ipsi

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

%% Generate other Y vectors containing pseudo data
%use linear-shift for non-independent vars with unknown distribution
%these go into Y.pseudo, below
trimLength = 20;

switch Ytype
    case {'RT' 'velocity' 'impulsive'}
        [~, whichTrials] = selectCondition(expInfo, getUniqueContrasts(expInfo), behavioralData, ...
            initTrialConditions('repeatType','random','specificRTs',[.1 3]));
        for l = 1:trimLength*2+1
            ss = l;
            es = length(mostTrials) - (trimLength*2-l+1);
            shiftVel(l,:) = trueVel(mostTrials(ss:es));
            shiftRT(l,:) = trueRT(mostTrials(ss:es));
            shiftImp(l,:) = trueImp(mostTrials(ss:es));
        end
        shiftVel(trimLength+1,:) = [];
        shiftRT(trimLength+1,:) = [];
        shiftImp(trimLength+1,:) = [];

    otherwise
        [~, whichTrials] = selectCondition(expInfo, getUniqueContrasts(expInfo), behavioralData, ...
            initTrialConditions('repeatType','random','movementTime','late','specificRTs',[.8 3]));
        for l = 1:trimLength*2+1
            ss = l;
            es = length(whichTrials) - (trimLength*2-l+1);
            shiftChoices(l,:) = trueChoices(whichTrials(ss:es));
            shiftSide(l,:) = trueSide(whichTrials(ss:es));
            shiftFeedback(l,:) = trueFeedback(whichTrials(ss:es));
            sohShifts(l,:,:) = stimOneHot(whichTrials(ss:es),:);
            cohShifts(l,:,:) = choiceOneHot(whichTrials(ss:es),:);
        end
        shiftChoices(trimLength+1,:) = [];
        shiftSide(trimLength+1,:) = [];
        shiftFeedback(trimLength+1,:) = [];
        sohShifts(trimLength+1,:,:) = [];
        cohShifts(trimLength+1,:,:) = [];
end

truncIdx = trimLength+1:length(whichTrials)-trimLength;

%% Set up X, Y, and D arrays for fitting

%subsample time range to save time
st = -1.5;
et = 3;
[~,si] = min(abs(neuralData.eta.eventWindow - st));
[~,ei] = min(abs(neuralData.eta.eventWindow - et));
timeRange = si:2:ei;  
% timeRange = 6:2:36;

%collect neural activity into the X matrix
% X: neural predictor matrix (trials x neurons x timepoints)
switch Ytype 
    case {'side' 'stimulus' 'block' 'value'}
        X = nan(length(whichTrials),size(neuralData.eta.alignedResps{ETA},3),length(timeRange));
        for t = 1:length(timeRange)
            X(:,:,t) = squeeze(neuralData.eta.alignedResps{ETA}(whichTrials,timeRange(t),:));
        end
    otherwise
        X = nan(length(truncIdx),size(neuralData.eta.alignedResps{ETA},3),length(timeRange));
        for t = 1:length(timeRange)
            X(:,:,t) = squeeze(neuralData.eta.alignedResps{ETA}(whichTrials(truncIdx),timeRange(t),:));
        end
end

%collect values for the Y matrices (and D arrays, if using)
% D: design matrix with task info, sometimes appended to X
% Y.true: outcome variable containing actual trial values
% Y.naive: predicted outcomes for a naive model. This contains a single
%     constant that equals the proportion of binary trials across the session
%     (e.g., prop right choices, mean RT, prop. correct trials, etc.)
% Y.pseudo: pseudo-generated outcomes (methods vary depending on the stats
%     of the variable)
% Y.pseudo_naive: predicted outcomes for a naive model using pseudo data
%     (same method as for Y.naive)

switch Ytype
    case 'side'
        Y.true = trueSide(whichTrials); 
        Y.naive = ones(1,length(whichTrials))*mean(trueSide(whichTrials));
        Y.pseudo = pseudoSide(:,whichTrials);
        Y.pseudo_naive = ones(1,length(whichTrials)).*mean(pseudoSide(:,whichTrials),2);
        family = 'binomial';

    case 'stimulus'
        Y.true = trueStimuli(whichTrials);
        Y.naive = ones(1,length(whichTrials))*mean(trueStimuli(whichTrials));
        Y.pseudo = pseudoStimuli(:,whichTrials);
        Y.pseudo_naive = ones(1,length(whichTrials)).*mean(pseudoStimuli(:,whichTrials),2);
        family = 'gaussian';

    case 'choice'
        Y.true = trueChoices(whichTrials(truncIdx));
        Y.naive = ones(1,length(truncIdx))*mean(trueChoices(whichTrials(truncIdx)));
        Y.pseudo = shiftChoices;
        Y.pseudo_naive = ones(1,length(truncIdx)).*mean(shiftChoices,2);
        family = 'binomial';

    case 'choiceGivenStim'
        Y.true = trueChoices(whichTrials(truncIdx));
        Y.pseudo = shiftChoices;
        D.true = stimOneHot(whichTrials(truncIdx),:);
        D.pseudo = sohShifts;
        family = 'binomial';

    case 'sideGivenChoice'
        Y.true = trueSide(whichTrials(truncIdx));
        Y.pseudo = shiftSide;
        D.true = choiceOneHot(whichTrials(truncIdx),:);
        D.pseudo = cohShifts;
        family = 'binomial';

    case 'feedback'
        Y.true = trueFeedback(whichTrials(truncIdx));
        Y.naive = ones(1,length(truncIdx))*mean(trueFeedback(whichTrials(truncIdx)));
        Y.pseudo = shiftFeedback;
        Y.pseudo_naive = ones(1,length(truncIdx)).*mean(shiftFeedback,2);
        family = 'binomial';

    case 'velocity'
        Y.true = trueVel(whichTrials(truncIdx));
        Y.naive = ones(1,length(truncIdx))*mean(trueVel(whichTrials(truncIdx)));
        Y.pseudo = shiftVel;
        Y.pseudo_naive = ones(1,length(truncIdx)).*mean(shiftVel,2);
        family = 'gaussian';

    case 'RT'
        Y.true = trueRT(whichTrials(truncIdx));
        Y.naive = ones(1,length(truncIdx))*mean(trueRT(whichTrials(truncIdx)));
        Y.pseudo = shiftRT;
        Y.pseudo_naive = ones(1,length(truncIdx)).*mean(shiftRT,2);
        family = 'gaussian';

    case 'impulsive'
        Y.true = trueImp(whichTrials(truncIdx));
        Y.naive = ones(1,length(truncIdx))*mean(trueImp(whichTrials(truncIdx)));
        Y.pseudo = shiftImp;
        Y.pseudo_naive = ones(1,length(truncIdx)).*mean(shiftImp,2);
        family = 'binomial';

    case 'block'
        Y.true = trueBlocks(whichTrials);
        Y.naive = ones(1,length(whichTrials))*mean(trueBlocks(whichTrials));
        Y.pseudo = pseudoBlocks(:,whichTrials);
        Y.pseudo_naive = ones(1,length(whichTrials)).*mean(pseudoBlocks(:,whichTrials),2);
        family = 'binomial';

    case 'value'
        Y.true = trueValue(whichTrials);
        Y.naive = ones(1,length(whichTrials))*mean(trueValue(whichTrials));
        Y.pseudo = pseudoValue(:,whichTrials);
        Y.pseudo_naive = ones(1,length(whichTrials)).*mean(pseudoValue(:,whichTrials),2);
        family = 'binomial';
end

%% Predict the true data: full model vs naive model

switch Ytype
    case {'choiceGivenStim' 'sideGivenChoice'}
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
        options.penalty_factor = [ones(1,size(X,2)) zeros(1,size(D.true,2))];

        %For the full model, we don't exclude any predictors. We want to
        %use both X and D
        options.exclude = [];
        prog = 0;
        fprintf(1,'fitting true (full model): %3d%%\n',prog);
        for t = 1:length(timeRange)
            preds = [X(:,:,t) D.true];
            fit.full{t} = cvglmnet(preds, Y.true', family, options, 'deviance',nf,[],false);
            info_full(t) = fit.full{t}.cvm(fit.full{t}.lambda == fit.full{t}.lambda_min)  / (-2 * log(2));
            prog = ( 100*(t/length(timeRange)));
            fprintf(1,'\b\b\b\b%3.0f%%',prog);
        end
        fprintf(1,'\n')
        
        %For the naive model, we exclude X and see how we do with D alone
        options.exclude = [1:size(X,2)]';
        options.lambda = ones(1,options.nlambda); %generate some meaningless lambdas to avoid glmnet bug
        prog = 0;
        fprintf(1,'fitting true (naive model): %3d%%\n',prog);
        preds = [X(:,:,1) D.true];
        fit.naive{1} = cvglmnet(preds, Y.true', family, options, 'deviance',nf,[],false);
        in = fit.naive{1}.cvm(1) / (-2 * log(2));
        info_naive = repmat(in, [1 length(timeRange)]);
        prog = ( 100*(t/length(timeRange)));
        fprintf(1,'\b\b\b\b%3.0f%%',prog);
        fprintf('\n');   

    otherwise
    %for all other Ytypes, we compare a full model that uses X to
    %predict true Y values, versus a naive model that uses a single
    %constant for all values of Y.
        clear options
        options.alpha = 1;
        options.nlambda = 20;
        options.standardize = 'false';
        nf = 3;
        
        %the full model uses X data and true Y data
        prog = 0;
        fprintf(1,'fitting true: %3d%%\n',prog);
        for t = 1:length(timeRange)
            fit.full{t} = cvglmnet(X(:,:,t),Y.true',family,options,'deviance',nf,[],true);
            Y.full(t,:) = cvglmnetPredict(fit.full{t},X(:,:,t),'lambda_min','response')';
            prog = ( 100*(t/length(timeRange)));
            fprintf(1,'\b\b\b\b%3.0f%%',prog);
        end
        fprintf('\n');      
        
        %compute LL for true vs full model
        for t = 1:length(timeRange)
            term1 = Y.true.*log(Y.full(t,:));
            term2 = (1-Y.true).*log(1-Y.full(t,:));
            ll_full(t) = sum(term1 + term2);
        end
        
        %we generated the naive model outcomes earlier
        %so compute LL for true vs naive model
        for t = 1:length(timeRange)
            term1 = Y.true.*log(Y.naive);
            term2 = (1-Y.true).*log(1-Y.naive);
            ll_naive(t) = sum(term1 + term2);
        end
        info_full = ll_full./size(Y.full,2)./log(2);
        info_naive = ll_naive./size(Y.naive,2)./log(2);
end

%% Predict the pseudo data: full model vs naive model

np = size(Y.pseudo,1);

switch Ytype
    case {'choiceGivenStim' 'sideGivenChoice'}
    %For these two Ytypes, we compare a full model that uses X (neural
    %data) and D (one-hot encoder holding trial stimulus or choice
    %identity), against a naive model that ONLY uses D.

    %But now both models are trying to predict pseudo Y data
        clear options
        options.alpha = 1;
        options.nlambda = 20;
        options.standardize = 'false';
        nf = 3;

        %X is always subject to a penalty (ones), but D is never
        %subject to a penalty (zeros)
        options.penalty_factor = [ones(1,size(X,2)) zeros(1,size(D.pseudo,3))];
        prog = 0;
        fprintf(1,'fitting pseudos: %3d%%\n',prog);
        for p = 1:np

            %For the full model, we don't exclude any predictors. We want to
            %use both X and D
            options.exclude = [];
            for t = 1:length(timeRange)
                try
                    preds = [X(:,:,t) squeeze(D.pseudo(p,:,:))];
                catch
                    preds = [X(:,:,t) squeeze(D.pseudo(p,:,:))'];
                end
                fit.pseudo{p,t} = cvglmnet(preds,Y.pseudo(p,:),family,options,'deviance',nf,[],false);
                info_pfull(p,t) = fit.pseudo{p,t}.cvm(fit.pseudo{p,t}.lambda == fit.pseudo{p,t}.lambda_min)  / (-2 * log(2));    
            end
            prog = (100*(p/np));
            fprintf(1,'\b\b\b\b%3.0f%%',prog);
        end
        fprintf('\n');
        
        %For the naive model, we exclude X and see how we do with D alone.
        %This is identical to the naive model we used on true data since D and Y arrays
        %don't shift against each other during pseudo-data generation
        info_pnaive = repmat(info_naive,[np 1]);
        
    otherwise
    %for all other Ytypes, we compare a full model that uses X to
    %predict pseudo Y values, versus a naive model that uses a single
    %constant for all values of Y.
        clear options
        options.alpha = 1;
        options.nlambda = 20;
        options.standardize = 'false';
        nf = 3;
        
        %the full model uses X data and pseudo Y data
        prog = 0;
        fprintf(1,'fitting pseudo: %3d%%\n',prog);
        for p = 1:np
            for t = 1:length(timeRange)
                fit.pseudo{p,t} = cvglmnet(X(:,:,t),Y.pseudo(p,:),family,options,'deviance',nf,[],true);
                Y.pseudo_full(t,p,:) = cvglmnetPredict(fit.pseudo{p,t},X(:,:,t),'lambda_min','response')';
            end
            prog = (100*(p/np));
            fprintf(1,'\b\b\b\b%3.0f%%',prog);
        end
        fprintf('\n');                
        
        for p = 1:np
            %compute LL for pseudo vs full model
            for t = 1:length(timeRange)
                term1 = Y.pseudo(p,:).*log(squeeze(Y.pseudo_full(t,p,:))');
                term2 = (1-Y.pseudo(p,:)).*log(1-squeeze(Y.pseudo_full(t,p,:))');
                llf(p,t) = sum(term1 + term2);
            end
            
            %we generated the naive model outcomes earlier
            %so compute LL for pseudo vs naive model
            for t = 1:length(timeRange)
                term1 = Y.pseudo(p,:).*log(Y.pseudo_naive(p,:));
                term2 = (1-Y.pseudo(p,:)).*log(1-Y.pseudo_naive(p,:));
                lln(p,t) = sum(term1 + term2);
            end
        end

        info_pfull = llf./size(Y.pseudo_full,3)./log(2);
        info_pnaive = lln./size(Y(1).pseudo_naive,2)./log(2);     
        
end

%% save

allX = X;
allY = Y;
try
    allD = D;
catch
    allD = [];
end

fits.true.full = fit.full;
try
    fits.true.naive = fit.naive;
catch
    fits.true.naive = 'constant model';
end
mutual_info.true = info_full - info_naive;

fits.pseudo.full = fit.pseudo;
try
    fits.pseudo.naive = fit.naive;
catch
    fits.pseudo.naive = 'constant model';
end
mutual_info.pseudo = info_pfull - info_pnaive;

%%
% 
% figure;
% ew = neuralData.eta.eventWindow(timeRange); 
% plot(ew,mutual_info_pseudo','Color',[.6 .6 .6])
% hold on 
% plot(ew,mutual_info_true,'k','LineWidth',2)
% line([0 0],[-.1 1.5],'Color',[.5 .5 .5],'LineStyle','--')
% ylim([-.05 1])
% xlim([ew(1) ew(end)])
% prettyPlot(gca)
% xlabel('Time from go cue (s)')
% ylabel('Mutual information')
% title(Ytype)
% %%
% 






