function [featureList, fitKernels] = unpackKernels(kernelFunctions, predictors, windows)

featureList = fieldnames(predictors);
for c = 1:size(kernelFunctions,2)
    ww = 1;
    for f = 1:length(featureList)
        if contains(featureList{f},'stimulus')
            wd = length(windows.stimulus);
        elseif contains(featureList{f},'action')
            wd = length(windows.action);
        elseif contains(featureList{f},'choice')
            wd = length(windows.choice);
        elseif contains(featureList{f},'velocity')
            wd = length(windows.velocity);
        elseif contains(featureList{f},'block')
            wd = length(windows.block);
        elseif contains(featureList{f},'outcome')
            wd = length(windows.outcome);
        elseif contains(featureList{f},'reward')
            wd = length(windows.outcome);
        elseif contains(featureList{f},'interaction')
            wd = length(windows.outcome);
        elseif contains(featureList{f},'value')
            wd = length(windows.value);
        elseif contains(featureList{f},'RT')
            wd = length(windows.RT);
        end
        fitKernels{c,f} = kernelFunctions(ww:ww+wd-1,c);
        ww = ww + wd;
    end
end