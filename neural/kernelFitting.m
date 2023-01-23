function [kernelFunctions, predictedActivity, optR, maxEV, cellFeatureStrength] = kernelFitting(behavioralData, expInfo, neuralData, features)

%% setup  
% design matrix
[predictors, windows] = getPredictors(...
    expInfo, ...
    behavioralData, ...
    features, ...
    0.1);

% set up regression matrices
X = makeToeplitz(neuralData.respTimes, predictors, windows);
Y = neuralData.cellResps;

%% RRR prep
% factorize to find basis functions for RRR
[B,~,~] = kFactorize(X, Y);

% turn on parpool
try 
    parpool();
catch
end

%% fit kernels for each neuron
%for each neuron, determine the optimal rank of XB to use to find its
%kernels, then compute how much variance is explained by this model

nFold = 5;
nComp = size(X,2);

% initialize variables
maxEV = zeros(size(Y,2),1);
optR = zeros(size(Y,2),1);
kernelFunctions = zeros(nComp,size(Y,2));

parfor c = 1:size(Y,2)
    
    %full-rank fit (nComp = [])
    [~, weightsK, explVarAll] = rrFit(X,Y(:,c),B,[],nFold,0);

    %determine the rank with the highest EV
    [maxEV(c,:), optR(c,:)] = max(explVarAll);

    %generate the kernels based on this rank
    kernelFunctions(:,c) = B(:,1:optR(c,:))*weightsK(:,1:optR(c,:))';
end

predictedActivity = X*kernelFunctions;


%% determine feature selectivity for each neuron

cellFeatureStrength = zeros(size(Y,2),length(features));

for f = 1:length(features)

    %determine the 'feature of interest' to test
    foi = features(f);

    %remove the foi from the feature set
    otherFeatures = features;
    otherFeatures(f) = [];

    %generate a predictor matrix without the foi
    [predictors_without, windows_without] = getPredictors(...
        expInfo, ...
        behavioralData, ...
        otherFeatures, ...
        0.1);
    
    %now generate a predictor matrix with ONLY the foi
    [predictors_foi, windows_foi] = getPredictors(...
        expInfo, ...
        behavioralData, ...
        foi, ...
        0.1);

    %predictor matrices
    X_without = makeToeplitz(neuralData.respTimes, predictors_without, windows_without);
    X_foi = makeToeplitz(neuralData.respTimes, predictors_foi, windows_foi);

    %factorize
    [B_without,~,~] = kFactorize(X_without, Y);
    [B_foi,~,~] = kFactorize(X_foi, Y);
    
    parfor c = 1:size(Y,2) 
                
        %fit without
        [Yhat_without, ~, ~] = rrFit(X_without, Y(:,c), B_without, optR(c), nFold, 0);

        %collect the residuals between true data and fitted
        residualRate = Y(:,c) - Yhat_without;

        %fit again using the residuals as Y
        [~, ~, EV_foi] = rrFit(X_foi, residualRate, B_foi, optR(c), nFold, 0);

        %collect the explained variance for this foi
        cellFeatureStrength(c,f) = EV_foi;

    end
end


%% WORKBENCH

%     for i = 1:nFold
% 
%         %split timepoints for cv
%         trainIdx = ismember(T, T(cvp.training(i)));
%         testIdx = ismember(T, T(cvp.test(i)));
% 
%         % outcome variable is a single neuron activity vector
%         Y_train = Y(trainIdx,c);
% 
%         % fit X*B (RR basis functions) to Y using the training set
%         fit = cvglmnet(X(trainIdx,:)*B,Y_train,'gaussian',options,'deviance',5,[],true);
% 
%         % evaluate fit weights at the desired lambda value
%         coefs = cvglmnetCoef(fit, 'lambda_min');
% 
%         %get rid of the intercept
%         weightsK = coefs(2:end)';
% 
%         %for a number of 'ranks' (from 1 to nComp), evaluate X*B for the test set
%         for nc = 1:nComp
%             cvPB(testIdx,nc) = X(testIdx,:)*(B(:,1:nc)*weightsK(:,1:nc)');
%         end
% 
%     end
% 
%     %compute the EV for each rank
%     explVarAll = zeros(nComp,1); corrAll = zeros(nComp,1);
%     for nc = 1:nComp
%         explVarAll(nc) = ...
%             1-(var(Y(:,c)-cvPB(:,nc)))./...
%             var(Y(:,c));
%             corrAll(nc) = corr(Y(:,c), cvPB(:,nc));
%     end







