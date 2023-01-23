function [cp, cp_shuffDist, cpSig, bp, bp_shuffDist, bpSig] = getCP(plotCells, resps, timeRange, behavioralData, expInfo)

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

%set some params
stimRange = 1:9;
numShuffles = 1000;

trialTypes = getTrialTypes(expInfo, behavioralData, 'late');

%initialize vars
cpu1 = nan(length(plotCells),length(stimRange), length(timeRange));
% cpu2 = nan(length(plotCells),length(stimRange), length(timeRange));
cpn = nan(length(plotCells),length(stimRange), length(timeRange));

cpu1_shuff = nan(length(plotCells),length(stimRange), length(timeRange), numShuffles);
% cpu2_shuff = nan(length(plotCells),length(stimRange), length(timeRange), numShuffles);
cpn_shuff = nan(length(plotCells),length(stimRange), length(timeRange), numShuffles);

bpu1 = nan(length(plotCells),length(stimRange), length(timeRange));
% bpu2 = nan(length(plotCells),length(stimRange), length(timeRange));
bpn = nan(length(plotCells),length(stimRange), length(timeRange));

bpu1_shuff = nan(length(plotCells),length(stimRange), length(timeRange), numShuffles);
% bpu2_shuff = nan(length(plotCells),length(stimRange), length(timeRange), numShuffles);
bpn_shuff = nan(length(plotCells),length(stimRange), length(timeRange), numShuffles);


for s = 1:length(stimRange)
    stimCond = stimRange(s);
    
    stimTrials = trialTypes.singleVar.contrast{stimCond};
        
    %pick the trial IDs corresponding to the two groups you want to compare
    cpTrials1 = stimTrials(behavioralData.wheelMoves.epochs(5).moveDir(stimTrials) == -1);
    cpTrials2 = stimTrials(behavioralData.wheelMoves.epochs(5).moveDir(stimTrials) == 1);
    cpTrialsAll = cat(2, cpTrials1,cpTrials2);
    
    bpTrials1 = stimTrials(expInfo.block.events.highRewardSideValues(stimTrials) == -1);
    bpTrials2 = stimTrials(expInfo.block.events.highRewardSideValues(stimTrials) == 1);
    bpTrialsAll = cat(2, bpTrials1,bpTrials2);
    
    %shuffles
    clear cprp bprp
    for iShuff = 1:numShuffles
        cprp(iShuff,:) = cpTrialsAll(randperm(length(cpTrialsAll),length(cpTrialsAll)));
        bprp(iShuff,:) = bpTrialsAll(randperm(length(bpTrialsAll),length(bpTrialsAll)));
    end
    
    tic
    for t = 1:length(timeRange)
        
        parfor iCell = 1:length(plotCells)
            
            %find the cell responses for those trial IDs and
            %compute the MWU value for the two groups
            [cpu1(iCell,s,t), ~, cpn(iCell,s,t)] = getMannWhitU(...
                resps(cpTrials1,timeRange(t),plotCells(iCell)),... %group1
                resps(cpTrials2,timeRange(t),plotCells(iCell))); %group2
            
               [bpu1(iCell,s,t), ~, bpn(iCell,s,t)] = getMannWhitU(...
                resps(bpTrials1,timeRange(t),plotCells(iCell)),... %group1
                resps(bpTrials2,timeRange(t),plotCells(iCell))); %group2

            %compute the MWU value for 2000 shuffles
            
            for iShuff = 1:numShuffles

                %split into shuffled groups and find cell responses and
                %compute the MWU value for the two groups

                [cpu1_shuff(iCell,s,t,iShuff), ...
                    ~, ...
                    cpn_shuff(iCell,s,t,iShuff)] = getMannWhitU(...
                    resps(cprp(iShuff, 1:length(cpTrials1)),timeRange(t),plotCells(iCell)),...
                    resps(cprp(iShuff, length(cpTrials1)+1:end),timeRange(t),plotCells(iCell)));
                
                [bpu1_shuff(iCell,s,t,iShuff), ...
                    ~, ...
                    bpn_shuff(iCell,s,t,iShuff)] = getMannWhitU(...
                    resps(bprp(iShuff, 1:length(bpTrials1)),timeRange(t),plotCells(iCell)),...
                    resps(bprp(iShuff, length(bpTrials1)+1:end),timeRange(t),plotCells(iCell)));
            end
            
        end
       
    end 
    toc
end


%sum all MWU stats (u1) and all possible comparisons (n) across conditions,
%then report the ratio
cp = squeeze(sum(cpu1,2)./sum(cpn,2));
bp = squeeze(sum(bpu1,2)./sum(bpn,2));

%do the same across conditions for the shuffles
cp_shuffDist = squeeze(sum(cpu1_shuff,2))./squeeze(sum(cpn_shuff,2));
bp_shuffDist = squeeze(sum(bpu1_shuff,2))./squeeze(sum(bpn_shuff,2));

%test significance (outside the 95% CI)
for iCell = 1:length(plotCells)
    cpSig(iCell,:) = cp(iCell,:) < prctile(cp_shuffDist(iCell,:),2.5) | cp(iCell) > prctile(cp_shuffDist(iCell,:),97.5);
    bpSig(iCell,:) = bp(iCell,:) < prctile(bp_shuffDist(iCell,:),2.5) | bp(iCell) > prctile(bp_shuffDist(iCell,:),97.5);
end

%%