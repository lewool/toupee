function [cp, cp_shuffDist, ifSig] = getCP(iCell, resps, trialArray)

%this function retrives the choice probability of a cell's responses at
%a particular epoch or timestamp, comparing trials in trialArray using a
%Mann-Whitney U statistic (getMannWhitU.m)

%iCell: a cellID from your experiment (this is an integer)

%resps: a vector of responses from either a specific epoch (e.g., from
%getEpochResps.m) or a particular ETA timestamp (e.g., 500ms post-stimON).
%This is a 1 x N vector, where N = numTrials in your experiment

%trialArray: a cell array that has split trials by stimCondition and
%choiceDirection. You can get this from getTrialTypes.m 
%(example: trialTypes.intVar.all.contrast_direction). You can use other
%conditions besides stimulus but it's not been well tested for this purpose
%so YMMV.

numShuffles = 2000;
alphaValue = 0.05;
CIIdx = [alphaValue/2 * numShuffles numShuffles - (alphaValue/2 * numShuffles)];

for stimCond = 1:size(trialArray, 1)
    %pick the trial IDs corresponding to the two groups you want to compare
    trials1 = trialArray{stimCond,1};
    trials2 = trialArray{stimCond,2};
    trialsAll = cat(2, trials1,trials2);
    
    %find the cell responses for those trial IDs
    group1 = resps(trials1,iCell);
    group2 = resps(trials2,iCell);
    
    %compute the MWU value for the two groups
    [u1(stimCond), u2(stimCond), n(stimCond)] = getMannWhitU(group1,group2);
    
    %compute the MWU value for 2000 shuffles
    for iShuff = 1:numShuffles
        %permute
        rp = trialsAll(randperm(length(trialsAll),length(trialsAll)));
        
        %split into shuffled groups and find cell responses
        shuff1 = resps(rp(1:length(trials1)),iCell);
        shuff2 = resps(rp(length(trials1)+1:end),iCell);
        
        %compute the MWU value for the two groups
        [u1_shuff(stimCond,iShuff), u2_shuff(stimCond,iShuff), n_shuff(stimCond,iShuff)] = getMannWhitU(shuff1,shuff2);
    end
    
end

%sum all MWU stats (u1) and all possible comparisons (n) across conditions,
%then report the ratio
cp = sum(u1)/sum(n);

%do the same across conditions for the shuffles
cp_shuffDist = sort(sum(u1_shuff)./sum(n_shuff));

%test significance (outside the 95% CI)
ifSig = cp < cp_shuffDist(CIIdx(1)) | cp > cp_shuffDist(CIIdx(2));
