function [Y_hat, weightsK, EV] = rrFit(X,Y,B,nComp,nFold,useCV)

% [B,~,~] = kFactorize(X, Y);

% try 
%     parpool();
% catch
% end

options.alpha = 0.5;
options.nlambda = 20;
options.standardize = true;

% partition for CV
cvp = cvpartition(size(X,1),'Kfold',nFold);
T = 1:size(X,1);
if ~isempty(nComp)
    if nComp > size(X,2)
        nComp = size(X,2);
    end
end

for i = 1:nFold
    
    %split timepoints for cv
    trainIdx = ismember(T, T(cvp.training(i)));
    testIdx = ismember(T, T(cvp.test(i)));
    
    % outcome variable is a single neuron activity vector
    Y_train = Y(trainIdx,:);

    if useCV
        % fit X*B (RR basis functions) to Y using the training set
        fit = cvglmnet(X(trainIdx,:)*B,Y_train,'gaussian',options,'deviance',3,[],true);

        % evaluate fit weights at the desired lambda value
        coefs = cvglmnetCoef(fit, 'lambda_min');
        
        %get rid of the intercept
        weightsK = coefs(2:end)';
    else 
        fit = glmnet(X(trainIdx,:)*B,Y_train,'gaussian',options);
        coefs = glmnetCoef(fit,[]); 
        weightsK = coefs(2:end,end)';
    end
    
    
    if isempty(nComp)
        %for a number of 'ranks' (from 1 to nc), evaluate X*B for the test set
        for nc = 1:size(B,2)
            Y_hat(testIdx,nc) = X(testIdx,:)*(B(:,1:nc)*weightsK(:,1:nc)');
        end
    else
        Y_hat(testIdx,:) = X(testIdx,:)*(B(:,1:nComp)*weightsK(:,1:nComp)');
    end
end

if isempty(nComp)
    EV = zeros(size(B,2),1); corrAll = zeros(size(B,2),1);
    for nc = 1:size(B,2)
        EV(nc) = ...
            1-(var(Y-Y_hat(:,nc)))./...
            var(Y);
        corrAll(nc) = corr(Y, Y_hat(:,nc));
    end
else
    EV = ...
            1-(var(Y-Y_hat))./...
            var(Y);
    corrAll = corr(Y, Y_hat);
end

% for nc = 1:nComp
%     EV(nc) = ...
%         1-(var(Y(:,c)-Y_hat(:,nc)))./...
%         var(Y(:,c));
%         corrAll(nc) = corr(Y(:,c), Y_hat(:,nc));
% end
