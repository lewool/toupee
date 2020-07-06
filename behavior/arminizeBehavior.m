function bigArminMatrix = arminizeBehavior(varargin)

allContrastValues = [];
allRewardSideValues = [];
allRepeatValues = [];
allHitValues = [];
allChoiceValues = [];
rs = [];
allBlockTrialNum = [];

    for s = 1:length(varargin{1})
        expInfo = varargin{1}(s);
        expInfo = data.loadExpData(expInfo);
        block = expInfo.block;
        
        trialLimit = length(block.events.endTrialValues);
        contrastValues = block.events.contrastValues(1:trialLimit);
        
        allContrastValues = [allContrastValues, contrastValues];
        
        try
            allRewardSideValues = [allRewardSideValues, block.events.likelyRewardSideValues(1:trialLimit)];
        catch
            try
                allRewardSideValues = [allRewardSideValues, block.events.highRewardSideValues(1:trialLimit)];
            catch
                try
                    rs(mod(block.events.contingencyPeriodValues,2) == 0) = -1;
                    rs(mod(block.events.contingencyPeriodValues,2) == 1) = 1;
                    allRewardSideValues = [allRewardSideValues, rs(1:trialLimit)];
                catch
                end   
            end
        end
        
        blockEnds = [find(diff(expInfo(1).block.events.blockSwitchesValues)) trialLimit];
        blockStarts = [1 find(diff(expInfo(1).block.events.blockSwitchesValues))+1];

        blockLengths = [blockEnds(1) diff(blockEnds)];
        blockNum = [];
        for b = 1:length(blockLengths)
            allBlockTrialNum = [allBlockTrialNum 1:blockLengths(b)];
        end

        allRepeatValues = [allRepeatValues, block.events.repeatNumValues(1:trialLimit)];
        allHitValues = [allHitValues, block.events.feedbackValues(1:trialLimit)];
        allChoiceValues = [allChoiceValues, block.events.responseValues(1:trialLimit)];
    end

allRewardSideValues(allRewardSideValues == 1) = 2;
allRewardSideValues(allRewardSideValues == -1) = 1;
    
bigArminMatrix = nan(length(allContrastValues),10);
bigArminMatrix(:,1) = allBlockTrialNum;
bigArminMatrix(:,2) = allContrastValues;
bigArminMatrix(:,3) = allChoiceValues;
bigArminMatrix(:,4) = allHitValues;
bigArminMatrix(:,8) = allRewardSideValues;
bigArminMatrix(:,10) = allHitValues;