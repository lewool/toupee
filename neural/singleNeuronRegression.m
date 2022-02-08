function singleNeuronRegression(X,Y)
% X is an array of behavioral predictors, such as contrasts or choices (one per trial).
% Y is a vector of cell/PC values (one per trial).

for tIdx = 1:size(Y,2)

    % set up the partitions for 10-fold cross-validation
    rng('default')
    C = cvpartition(size(X,1),'KFold', 10);
    
    B_fold = [];
    dev_fold = [];
    for iCV = 1:C.NumTestSets
        [B_fold(:,iCV), dev_fold(iCV), stats] = glmfit(X(training(C,iCV),:),Y(training(C,iCV),tIdx));
        Y_pred{iCV} = glmval(B_fold(:,iCV),X(test(C,iCV),:),'identity');
        SSE_fold(iCV) = sum((Y_pred{iCV}-Y(test(C,iCV),tIdx)).^2);
    end
    
    SSE(tIdx) = nanmean(SSE_fold);
    dev(tIdx) = nanmean(dev_fold);
    weights(:,tIdx) = nanmean(B_fold,2);
end
%%
figure;
plot(neuralData.eta.eventWindow, (weights(2,:)))
hold on
plot(neuralData.eta.eventWindow, (weights(3,:)))
plot(neuralData.eta.eventWindow, (weights(4,:)))