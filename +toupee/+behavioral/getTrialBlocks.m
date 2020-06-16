function trialBlocks = getTrialBlocks(expInfo)
tbl = 0;
for ex = 1:length(expInfo)
    
    nt = length(expInfo(ex).block.events.endTrialValues);
    allBlocks = unique(expInfo(ex).block.events.blockSwitchesValues);

    for b = 1:length(allBlocks)
        trialBlocks{tbl + b,1} =  find(expInfo(ex).block.events.blockSwitchesValues(1:nt) == allBlocks(b));
        trialBlocks{tbl + b,2} = expInfo(ex).block.events.highRewardSideValues(trialBlocks{b,1}(1));
    end

    tbl = size(trialBlocks,1);

end
