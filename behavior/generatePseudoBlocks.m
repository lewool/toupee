function pseudoBlocks = generatePseudoBlocks(expInfo, np, blockStart)

nt = size(expInfo.block.events.endTrialValues,2);
pseudoBlocks = nan(np,nt);
for p = 1:np
    if strcmp(blockStart,'fixed')
        firstSide = expInfo.block.paramsValues(1).firstHighSide;
    elseif strcmp(blockStart,'rand')
        firstSide = randsample([-1, 1],1,true);
    end
    b=nan(1,nt);
    switches = cumsum(125+randi(100,1,20));
    for s = 1:length(switches)
        if s == 1
            b(1:switches(s)-1) = firstSide;
        elseif mod(s,2) == 1
            b(switches(s-1):switches(s)-1) = firstSide;
        elseif mod(s,2) == 0
            b(switches(s-1):switches(s)-1) = -firstSide;
        end
    end
    pseudoBlocks(p,:) = b(1:nt);
end

