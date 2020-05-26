function [meanResp, semResp, contrastColors] = computeContrastCurves(cellArray,expInfo,behavioralData,trialCondition,whichContrasts)

% Cell array is a 2D trials x window array of responses for a single cell
dirs = {'cw' 'ccw'};
contrasts = getUniqueContrasts(expInfo);
if strcmp(whichContrasts,'each')
    for c = 1:length(contrasts)
        % select trials based on contrast
        [~, condIdx] = selectCondition(expInfo, contrasts(c), behavioralData, trialCondition);

        %compute mean+sem responses
        meanResp(c,:) = nanmean(cellArray(condIdx,:),1);
        semResp(c,:) = nanstd(cellArray(condIdx,:))/sqrt(length(condIdx));
    end

    allColors = [0 0 .5; 0 0 1; 0 .4 1; .6 .8 1; .75 .75 .75; .8 .45 .45; 1 0 0; .5 0 0; .25 0 0];

    zeroIdx = find(contrasts == 0);
    walkup = length(contrasts) - zeroIdx;
    walkback = zeroIdx - 1;
    zeroGray = find(allColors(:,1) == .75);
    contrastColors = allColors(zeroGray-walkback:zeroGray + walkup,:);
    

    
elseif strcmp(whichContrasts,'side')
    cIdx = [contrasts < 0; contrasts == 0; contrasts > 0];
    for c = 1:size(cIdx,1)
        % select trials based on side
        [~, condIdx] = selectCondition(expInfo, contrasts(cIdx(c,:)), behavioralData, trialCondition);
        
        %compute mean+sem responses
        meanResp(c,:) = nanmean(cellArray(condIdx,:),1);
        semResp(c,:) = nanstd(cellArray(condIdx,:))/sqrt(length(condIdx));
    end
    
    contrastColors = [0 .4 1; .75 .75 .75; 1 0 0];
    
elseif strcmp(whichContrasts,'balanced')
    for d = 1:length(dirs)
        for c = 1:length(contrasts)
            [~, condIdx{c,d}] = selectCondition(expInfo, contrasts(c), behavioralData, trialCondition);
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
end
