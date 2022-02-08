function neuralChoiceDecoder(neuralData,behavioralData)

[baselineResps, stimResps, pmovResps, movResps, rewResps, preCueResps] = getEpochResps(neuralData(1).eta);

%%
%select neural data. observations x predictors, aka trials x cells (pick a
%single time point to do the analysis)
% Xtest = rewResps(2:2:end,:);
% Xtrain = rewResps(1:2:end,:);

%select binomial data. i use the direction of the first movement on each
%trial. 1 = right, 0 = left
Y = double((behavioralData.wheelMoves.epochs(5).moveDir == 1)');
[~, lateTrials] = selectCondition(expInfo, getUniqueContrasts(expInfo), behavioralData, initTrialConditions('movementTime','late'));

Ylate = Y(lateTrials);
Ytest = Ylate(2:2:end,:);
Ytrain = Ylate(1:2:end,:);

%%

%     Y = (expInfo.block.events.highRewardSideValues(1:401) == 1)';

lambdas = [0, 1e-10, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-3, 1e-2, 1e-1, 1];
% Set up glmnet options. Pick our lambdas manually to be sure of a broad range
options = glmnetSet(); options.lambda = lambdas;
timerange = 1:41;
for timebin = 5:length(timerange)
    Xlate = squeeze(neuralData.eta.alignedResps{1}(lateTrials,timerange(timebin),:));
    Xtest = Xlate(2:2:end,:);
    Xtrain = Xlate(1:2:end,:);

    fit{timebin} = glmnet(Xtrain,Ytrain,'binomial',options);

    %Pick the output with the best deviance. This means the weakest
    % regularization that still converged
    dev_by_lambda = (1-fit{timebin}.dev) * fit{timebin}.nulldev;
    [best_dev, best_ind] = min(dev_by_lambda);
            
 % Record weights and devs
    weights = fit{timebin}.beta(:,best_ind);
    devs = best_dev;
    nulldevs = fit{timebin}.nulldev;

    pred{timebin} = glmnetPredict(fit{timebin},Xtest,[],'response');
    
    Ypred = pred{timebin}(:,best_ind) <= 0.5;
    acc(timebin) = sum(Ypred == Ytest)/length(Ytest);
    pause
end

figure;
plot(neuralData.eta.eventWindow,smooth(1-acc,1))

%%

timerange = 21;
for timebin = 1:length(timerange)
    Xlate = squeeze(neuralData.eta.alignedResps{1}(lateTrials,timerange(timebin),:));
    Xtest = Xlate(2:2:end,:);
    Xtrain = Xlate(1:2:end,:);

    cvfit{timebin} = cvglmnet(Xlate,Ylate,'binomial');

%     %Pick the output with the best deviance. This means the weakest
%     % regularization that still converged
%     dev_by_lambda = (1-cvfit{timebin}.glmnet_fit.dev) * cvfit{timebin}.glmnet_fit.nulldev;
%     [best_dev, best_ind] = min(dev_by_lambda);
%             
%  % Record weights and devs
%     weights = cvfit{timebin}glmnet_fit.beta(:,best_ind);
%     devs = best_dev;
%     nulldevs = cvfit{timebin}.glmnet_fit.nulldev;

    pred{timebin} = cvglmnetPredict(cvfit{timebin},Xlate,[],'response');
    
    Ypred = pred{timebin}(:,timebin) <= 0.5;
    acc(timebin) = sum(Ypred == Ytest)/length(Ytest);
end



%%


for m = 1:length(mouseList)
    
    mouseName = char(mouseList{m});
    expDate = char(expList{m}{1});
    expNum = expList{m}{2};
    [expRef, ~] = data.constructExpRef(mouseName,expDate,expNum);
    fprintf('fitting %s...\n',expRef)
    [behavioralData, expInfo, neuralData] = data.loadDataset(mouseName, expDate, expNum);
    [baselineResps, stimResps, pmovResps, movResps, rewResps, preCueResps] = getEpochResps(neuralData.eta);
    
    %set up binomial Y outcome vector (chose L vs R)
%     Y_all = double((behavioralData.wheelMoves.epochs(5).moveDir == 1)');
    
    nt = length(behavioralData.eventTimes(1).daqTime);
    if hemisphere > 0
        Y_all = double(expInfo.block.events.highRewardSideValues(1:nt) == -1);
    else
        Y_all = double(expInfo.block.events.highRewardSideValues(1:nt) == 1);
    end
    
    %filter out the early trials when animal broke quiescence
%     [~, lateTrials] = selectCondition(expInfo, getUniqueContrasts(expInfo), behavioralData, initTrialConditions('movementTime','late'));
%     [~, earlyTrials] = selectCondition(expInfo, getUniqueContrasts(expInfo), behavioralData, initTrialConditions('movementTime','early'));
% 
%     Y = Y_all(lateTrials);
    Y=Y_all;

    %split into test/train for simple CV
    Ytest = Y(2:2:end,:);
    Ytrain = Y(1:2:end,:);

     %set up some lambdas
    lambdas = [0, 1e-10, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-3, 1e-2, 1e-1, 1];

    %range of timebins in the trial to use (idx=21 is event onset)
    timerange = 1:41;
    for t = 1:length(timerange)
        %retrieve the activity of each cell on every trial, at timebin t
%         X = squeeze(neuralData.eta.alignedResps{2}(lateTrials,timerange(t),:));
        X = 
        
        %split into test/train for simple CV
        Xtest = X(2:2:end,:);
        Xtrain = X(1:2:end,:);

        [B,FitInfo] = lassoglm(Xtrain,Ytrain,'binomial','Lambda',lambdas,'Alpha',1,'Standardize',false);
    %     [B,FitInfo] = glmfit(Xtrain,Ytrain,'binomial');
    %     bestIdx = max(find(sum(B>0)>0));

        Y_hat_all = glmval([FitInfo.Intercept;B],Xtest,'logit');
        acc_all = sum((Y_hat_all>0.5) == Ytest)/length(Ytest);
        [~, bestIdx] = max(acc_all);
        accuracy(m, t) = acc_all(bestIdx);
        fprintf('timebin %d completed\n',t)

    end

    clearvars -except mouseList expList accuracy
end

% figure;
% plot(neuralData.eta.eventWindow, accuracy)











 %%
for tIdx = 1:length(neuralData.eta.eventWindow)
    %select neural data. observations x predictors, aka trials x cells (pick a
    %single time point to do the analysis)
    Xtest = squeeze(neuralData.eta.alignedPCs{2}(:,tIdx,1:20));

    %select binomial data. i use the direction of the first movement on each
    %trial. 1 = right, 0 = left
    Y = (behavioralData.wheelMoves.epochs(5).moveDir == 1)';
%     Y = (expInfo.block.events.highRewardSideValues(1:401) == 1)';

    % set up the partitions for 10-fold cross-validation
    rng('default')
    C = cvpartition(size(Xtest,1),'KFold', 10);

    for iCV = 1:C.NumTestSets
        [B,FitInfo] = lassoglm(Xtest(training(C,iCV),:),Y(training(C,iCV)),'binomial');
        indx = 1;
        B0 = B(:,indx);
        cnst = FitInfo.Intercept(indx);
        B1(tIdx,iCV,:) = [cnst;B0];
        Y_pred{iCV} = glmval(squeeze(B1(tIdx,iCV,:)),Xtest(test(C,iCV),:),'logit');
        predAcc(iCV) = sum((Y_pred{iCV} > .5) == Y(test(C,iCV))) / length(Y(test(C,iCV)));
    end
    
    accuracy(tIdx) = mean(predAcc);
end
weights = squeeze(mean(B1,2));

figure;
plot(neuralData.eta.eventWindow, smooth(accuracy))

%% 
for tIdx = 1:length(neuralData.eta.eventWindow)

    Xtest = squeeze(neuralData.eta.alignedResps{2}(2,tIdx,:));
    Xtest = squeeze(neuralData.eta.alignedPCs{2}(:,tIdx,1:20));
    DV(tIdx,:) = weights(tIdx,1) + weights(tIdx, 2:end)*Xtest';
%     meanX(tIdx) = mean(X); 
end
%%
figure;
plot(neuralData.eta.eventWindow,nanmean(DV(:,behavioralData.wheelMoves.epochs(5).moveDir == -1),2))
hold on;
