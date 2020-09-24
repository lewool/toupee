function [meanResp, semResp, rewColors] = computeRewardCurves(cellArray,expInfo,behavioralData,trialCondition,whichContrasts)

% Cell array is a 2D trials x window array of responses for a single cell

contrasts = getUniqueContrasts(expInfo);
rews = {'correct' 'incorrect'};

if strcmp(whichContrasts,'balanced')
    for d = 1:length(rews)
        for c = 1:length(contrasts)
            [~, condIdx{c,d}] = selectCondition(expInfo, contrasts(c), behavioralData, initTrialConditions('movementTime','late','responseType',rews{d}));
        end
    end
    for c = 1:length(condIdx)
        minc = min([length(condIdx{c,1}) length(condIdx{c,2})]);
        condIdx{c,1} = randsample(condIdx{c,1},minc);
        condIdx{c,2} = randsample(condIdx{c,2},minc);
    end
    if strcmp(trialCondition.movementDir,'cw')
        lIdx = cat(2, condIdx{contrasts<0,1});
        zIdx = cat(2, condIdx{contrasts==0,1});
        rIdx = cat(2, condIdx{contrasts>0,1});
        meanResp(1,:) = nanmean(cellArray(lIdx,:),1);
        semResp(1,:) = nanstd(cellArray(lIdx,:))/sqrt(length(lIdx));
        meanResp(2,:) = nanmean(cellArray(zIdx,:),1);
        semResp(2,:) = nanstd(cellArray(zIdx))/sqrt(length(zIdx));
        meanResp(3,:) = nanmean(cellArray(rIdx,:),1);
        semResp(3,:) = nanstd(cellArray(rIdx,:))/sqrt(length(rIdx));
    elseif strcmp(trialCondition.movementDir,'ccw')
        lIdx = cat(2, condIdx{contrasts<0,2});
        zIdx = cat(2, condIdx{contrasts==0,2});
        rIdx = cat(2, condIdx{contrasts>0,2});
        meanResp(1,:) = nanmean(cellArray(lIdx,:),1);
        semResp(1,:) = nanstd(cellArray(lIdx,:))/sqrt(length(lIdx));
        meanResp(2,:) = nanmean(cellArray(zIdx,:),1);
        semResp(2,:) = nanstd(cellArray(zIdx))/sqrt(length(zIdx));
        meanResp(3,:) = nanmean(cellArray(rIdx,:),1);
        semResp(3,:) = nanstd(cellArray(rIdx,:))/sqrt(length(rIdx));
    end
    
    contrastColors = [0 .4 1; .75 .75 .75; 1 0 0];
elseif strcmp(whichContrasts,'each')
for d = 1:length(rews)
    trialCondition.responseType = rews{d};
    % select trials based on high-reward side x contrast
    [~, condIdx] = selectCondition(expInfo, contrasts, behavioralData, trialCondition);

    %compute mean+sem responses
    meanResp(d,:) = nanmean(cellArray(condIdx,:),1);
    semResp(d,:) = nanstd(cellArray(condIdx,:))/sqrt(length(condIdx));
end
end

rewColors = [.1 .7 .1; .75 0 0];