function [meanResp, semResp, rewColors] = computeRewardCurves(cellArray,expInfo,behavioralData,trialCondition)

% Cell array is a 2D trials x window array of responses for a single cell

contrasts = getUniqueContrasts(expInfo);
rews = {'correct' 'incorrect'};

for d = 1:length(rews)
    trialCondition.responseType = rews{d};
    % select trials based on high-reward side x contrast
    [~, condIdx] = selectCondition(expInfo, contrasts, behavioralData, trialCondition);

    %compute mean+sem responses
    meanResp(d,:) = nanmean(cellArray(condIdx,:),1);
    semResp(d,:) = nanstd(cellArray(condIdx,:))/sqrt(length(condIdx));
end

rewColors = [.1 .7 .1; .75 0 0];