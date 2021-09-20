function neuralChoiceDecoder(neuralData,behavioralData)

for tIdx = 1:length(neuralData.eta.eventWindow)
    %select neural data. observations x predictors, aka trials x cells (pick a
    %single time point to do the analysis)
    X = squeeze(neuralData.eta.alignedPCs{2}(:,tIdx,1:20));

    %select binomial data. i use the direction of the first movement on each
    %trial. 1 = right, 0 = left
    Y = (behavioralData.wheelMoves.epochs(5).moveDir == 1)';
%     Y = (expInfo.block.events.highRewardSideValues(1:401) == 1)';

    % set up the partitions for 10-fold cross-validation
    rng('default')
    C = cvpartition(size(X,1),'KFold', 10);

    for iCV = 1:C.NumTestSets
        [B,FitInfo] = lassoglm(X(training(C,iCV),:),Y(training(C,iCV)),'binomial');
        indx = 1;
        B0 = B(:,indx);
        cnst = FitInfo.Intercept(indx);
        B1(tIdx,iCV,:) = [cnst;B0];
        Y_pred{iCV} = glmval(squeeze(B1(tIdx,iCV,:)),X(test(C,iCV),:),'logit');
        predAcc(iCV) = sum((Y_pred{iCV} > .5) == Y(test(C,iCV))) / length(Y(test(C,iCV)));
    end
    
    accuracy(tIdx) = mean(predAcc);
end
weights = squeeze(mean(B1,2));

figure;
plot(neuralData.eta.eventWindow, smooth(accuracy))

%% 
for tIdx = 1:length(neuralData.eta.eventWindow)

    X = squeeze(neuralData.eta.alignedResps{2}(2,tIdx,:));
    X = squeeze(neuralData.eta.alignedPCs{2}(:,tIdx,1:20));
    DV(tIdx,:) = weights(tIdx,1) + weights(tIdx, 2:end)*X';
%     meanX(tIdx) = mean(X); 
end
%%
figure;
plot(neuralData.eta.eventWindow,nanmean(DV(:,behavioralData.wheelMoves.epochs(5).moveDir == -1),2))
hold on;
