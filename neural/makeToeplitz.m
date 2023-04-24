function toeplitzMatrix = makeToeplitz(respTimes, predictors, windows)
%Extract the trial-by-trial activity for ROIs in each imaging plane
%and align to a particular trial event (stim on, movement on, etc.)

%% make matrix
featureList = fieldnames(predictors);
    
for f = 1:length(featureList)
    if contains(featureList{f},'stimulus')
        wd = windows.stimulus;
    elseif contains(featureList{f},'action')
        wd = windows.action;
    elseif contains(featureList{f},'choice')
        wd = windows.choice;
    elseif contains(featureList{f},'velocity')
        wd = windows.velocity;
    elseif contains(featureList{f},'block')
        wd = windows.block;
    elseif contains(featureList{f},'outcome')
        wd = windows.outcome;
    elseif contains(featureList{f},'reward')
        wd = windows.outcome;
    elseif contains(featureList{f},'value')
        wd = windows.value;
    elseif contains(featureList{f},'RT')
        wd = windows.RT;
    end

    if contains(featureList{f},'velocity')
        vels = interp1(predictors.(featureList{f}).times, predictors.(featureList{f}).values, respTimes);
        vels(isnan(vels)) = 0;
        tplz{1,f} = zeros(size(respTimes,2),length(wd));
        for w = 1:length(wd)
            tplz{1,f}(:,w) = circshift(vels,wd(w));
            tplz{2,f} = wd;
        end
    else
        %interpolate predictor times to planeFrameTimes
        predTimes = interp1(respTimes, respTimes, predictors.(featureList{f}).times, 'nearest');
        %index which elements of planeFrameTimes correspond to
        %predictor times
        [ptimes , planeIdx] = intersect(respTimes,predTimes);
        % in theory this should be the same as above but sometimes there
        % are NaNs or duplicate values in predTimes which are
        % dropped...this tracks which ones we kept 
        [~, predIdx] = intersect(predTimes,ptimes);
        tplz{1,f} = zeros(size(respTimes,2),length(wd));
        for w = 1:length(wd)
            tplz{1,f}(planeIdx+wd(w),w) = predictors.(featureList{f}).values(predIdx);
            tplz{2,f} = wd;
        end
    end
end

toeplitzMatrix = [cat(2,tplz{1,1:length(featureList)})];

    
