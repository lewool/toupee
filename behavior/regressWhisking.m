%% 

for ex = 1:length(expInfo)
    numTrials(ex) = length(expInfo(ex).block.events.endTrialValues);
end

trimLength = min(numTrials);
trimLength = 909;
%% LINEAR REGRESSION

% assign X. this is a 1 x n vector of values from your predictor feature 
% (i.e., pre-trial whisking). you should have 1 value per trial 
% (so n = number of trials in that session)
X =  mean(eyeData(1).eta.alignedFace{1}(:,91:101,2),2);
% X =  (mean(eyeData(1).eta.alignedFace{1}(1:trimLength,91:101,2),2));
% X(:,2) =  expInfo.block.events.contrastValues(1:909);

% assign Y. this is a 1 x n vector of outcomes. for the linear regression,
% use the RT (firstMov - stimOn) on each trial
Y = behavioralData(1).wheelMoves.epochs(5).onsetTimes' - behavioralData(1).eventTimes(1).daqTime';
% Y = (behavioralData(1).wheelMoves.epochs(5).onsetTimes(1:trimLength)' - behavioralData(1).eventTimes(1).daqTime(1:trimLength)');

% set up the partitions for 10-fold cross-validation
rng('default')
C = cvpartition(size(X,1),'KFold', 10);

B_fold = [];
dev_fold = [];
stats_fold = [];

for iCV = 1:C.NumTestSets
    [B_fold(:,iCV), dev_fold(iCV), stats_fold{iCV}] = glmfit(X(training(C,iCV),:),Y(training(C,iCV)));
    Y_pred{iCV} = glmval(B_fold(:,iCV),X(test(C,iCV),:),'identity');
    SSE_fold(iCV) = sum((Y_pred{iCV}-Y(test(C,iCV))).^2);
end

SSE = nanmean(SSE_fold);
dev = nanmean(dev_fold);
weights = nanmean(B_fold,2);


for iPerm = 2:length(expInfo)
    X =  diff(mean(eyeData(iPerm).eta.alignedFace{1}(1:trimLength,91:101,2),2));
    rng('default')
    C = cvpartition(size(X,1),'KFold', 10);

    B_fold = [];
    dev_fold = [];
    stats_fold = [];

    for iCV = 1:C.NumTestSets
        [B_fold(:,iCV), dev_fold(iCV), stats_fold{iCV}] = glmfit(X(training(C,iCV),:),Y(training(C,iCV)));
        Y_pred{iCV} = glmval(B_fold(:,iCV),X(test(C,iCV),:),'identity');
    end

    perm_dev(iPerm) = nanmean(dev_fold);
    perm_weights(:,iPerm) = nanmean(B_fold,2);
end

figure;
histogram(perm_weights(2,:));
hold on
line([weights(2,:) weights(2,:)],[0 5],'LineWidth',2,'Color',[1 .25 0])
    
%% LOGISTIC REGRESSION

% assign X. this is a 1 x n vector of values from your predictor feature 
% (i.e., pre-trial whisking). you should have 1 value per trial 
% (so n = number of trials in that session)
X =  mean(eyeData.eta.alignedFace{1}(:,91:101,2),2);
X = shuffleVec(X);

% assign Y. this is a 1 x n vector of outcomes. for the logistic regression,
% use early (1) vs late (0)
Y = behavioralData.wheelMoves.epochs(2).isMoving';
% Y = shuffleVec(Y);
% set up the partitions for 10-fold cross-validation
    rng('default')
    C = cvpartition(size(X,1),'KFold', 10);
    
    B_fold = [];
    dev_fold = [];
    for iCV = 1:C.NumTestSets
        [B_fold(:,iCV), dev_fold(iCV), stats] = glmfit(X(training(C,iCV)),Y(training(C,iCV)),'binomial');
        Y_pred{iCV} = glmval(B_fold(:,iCV),X(test(C,iCV)),'logit');
        predAcc(iCV) = sum((Y_pred{iCV} > .5) == Y(test(C,iCV))) / length(Y(test(C,iCV)));
    end
    
    accuracy = nanmean(predAcc);
    dev = nanmean(dev_fold);
    weights = nanmean(B_fold,2);
    
    
    