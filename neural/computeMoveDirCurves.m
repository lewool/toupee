function [meanResp, semResp, dirColors] = computeMoveDirCurves(cellArray,expInfo,behavioralData,trialCondition,whichContrasts)

% Cell array is a 2D trials x window array of responses for a single cell

contrasts = getUniqueContrasts(expInfo);
dirs = {'cw' 'ccw'};

if nargin < 5
    for d = 1:length(dirs)
        trialCondition.movementDir = dirs{d};
        % select trials based on high-reward side x contrast
        [~, condIdx] = selectCondition(expInfo, contrasts, behavioralData, trialCondition);

        %compute mean+sem responses
        meanResp(d,:) = nanmean(cellArray(condIdx,:),1);
        semResp(d,:) = nanstd(cellArray(condIdx,:))/sqrt(length(condIdx));
    end
else
    for d = 1:length(dirs)
        trialCondition.movementDir = dirs{d};
        if strcmp(whichContrasts,'left')
        cc = contrasts(contrasts<0);
        elseif strcmp(whichContrasts,'right')
            cc = contrasts(contrasts<0);
        end
        for c = 1:length(cc)
            [~, condIdx{c,d}] = selectCondition(expInfo, cc(c), behavioralData, trialCondition);
        end
    end
    for c = 1:length(condIdx)
        minc = min([length(condIdx{c,1}) length(condIdx{c,2})]);
        condIdx{c,1} = randsample(condIdx{c,1},minc);
        condIdx{c,2} = randsample(condIdx{c,2},minc);
    end
    for d = 1:length(dirs)
        meanResp(d,:) = nanmean(cellArray(cat(2, condIdx{:,d}),:),1);
        semResp(d,:) = nanstd(cellArray(cat(2, condIdx{:,d}),:))/sqrt(length(cat(2, condIdx{:,d})));
    end
end

dirColors = [0 .4 1; 1 0 0];