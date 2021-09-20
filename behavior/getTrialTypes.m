function trialTypes = getTrialTypes(expInfo, behavioralData, movementTime)
%This function divides and counterbalances all trial 
%types based on certain task variables (stimulus, mvmt direction, outcome,
%and reward block). You can then call one of these cell arrays anywhere
%you need a list of trial numbers 

% single Var = trials sorted by single variables
% intVar = trials sorted by interacting variables
%   all = no counterbalancing, all trials that pass criteria
%   cb3D = 3-way counterbalancing based on stim x direction x block
%   cb2D = 2-way counterbalancing based on stim x direction

outcomes = {'correct' 'incorrect'};
dirs = {'cw' 'ccw'};
blocks = {'left' 'right'};
contrasts = getUniqueContrasts(expInfo);

%%%%%%% SINGLE VARIABLES %%%%%%%%%%%%%%%%
%collects all trials (that pass the inputVar criteria) by the following
%task dimensions:
if strcmp(movementTime,'>0.5')
    %by stimulus contrast
    for i = 1:length(contrasts)
        [~, co{i,1}] = selectCondition(expInfo, contrasts(i), behavioralData, initTrialConditions('movementTime','all','specificRTs',[.5 2]));
    end
    %by stimulus side
    [~, si{1,1}] = selectCondition(expInfo, contrasts(contrasts<0), behavioralData, initTrialConditions('movementTime','all','specificRTs',[.5 2]));
    [~, si{2,1}] = selectCondition(expInfo, contrasts(contrasts==0), behavioralData, initTrialConditions('movementTime','all','specificRTs',[.5 2]));
    [~, si{3,1}] = selectCondition(expInfo, contrasts(contrasts>0), behavioralData, initTrialConditions('movementTime','all','specificRTs',[.5 2]));
    %by movement direction
    for i = 1:length(dirs)
        [~, di{i,1}] = selectCondition(expInfo, contrasts, behavioralData, initTrialConditions('movementTime','all','specificRTs',[.5 2],'movementDir',dirs{i}));
    end
    %by outcome
    for i = 1:length(outcomes)
        [~, ou{i,1}] = selectCondition(expInfo, contrasts, behavioralData, initTrialConditions('movementTime','all','specificRTs',[.5 2],'responseType',outcomes{i}));
    end
    %by block
    for i = 1:length(blocks)
        [~, bl{i,1}] = selectCondition(expInfo, contrasts, behavioralData, initTrialConditions('movementTime','all','specificRTs',[.5 2],'highRewardSide',blocks{i}));
    end
else
    %by stimulus contrast
    for i = 1:length(contrasts)
        [~, co{i,1}] = selectCondition(expInfo, contrasts(i), behavioralData, initTrialConditions('movementTime',movementTime));
    end
    %by stimulus side
    [~, si{1,1}] = selectCondition(expInfo, contrasts(contrasts<0), behavioralData, initTrialConditions('movementTime',movementTime));
    [~, si{2,1}] = selectCondition(expInfo, contrasts(contrasts==0), behavioralData, initTrialConditions('movementTime',movementTime));
    [~, si{3,1}] = selectCondition(expInfo, contrasts(contrasts>0), behavioralData, initTrialConditions('movementTime',movementTime));
    %by movement direction
    for i = 1:length(dirs)
        [~, di{i,1}] = selectCondition(expInfo, contrasts, behavioralData, initTrialConditions('movementTime',movementTime,'movementDir',dirs{i}));
    end
    %by outcome
    for i = 1:length(outcomes)
        [~, ou{i,1}] = selectCondition(expInfo, contrasts, behavioralData, initTrialConditions('movementTime',movementTime,'responseType',outcomes{i}));
    end
    %by block
    for i = 1:length(blocks)
        [~, bl{i,1}] = selectCondition(expInfo, contrasts, behavioralData, initTrialConditions('movementTime',movementTime,'highRewardSide',blocks{i}));
    end
end
%log
trialTypes.singleVar = struct(...
    'contrast',{co},...
    'side',{si},...
    'direction',{di},...
    'outcome',{ou},...
    'block',{bl}...
    );

%%%%%%% INTERACTIONS (UNBALANCED) %%%%%%%%%%%%%%%%%%%

%split trials as CONTRAST x DIRECTION (nContrasts x nDirs)       
for d = 1:length(dirs)
    for c = 1:length(contrasts)
        if strcmp(movementTime,'>0.5')
            [~, co_di{c,d}] = selectCondition(expInfo, contrasts(c), behavioralData, initTrialConditions('movementTime','all','specificRTs',[.5 Inf],'movementDir',dirs{d}));
        else
            [~, co_di{c,d}] = selectCondition(expInfo, contrasts(c), behavioralData, initTrialConditions('movementTime',movementTime,'movementDir',dirs{d}));
        end
    end
end
for d = 1:size(co_di,2)
    si_di{1,d} = cat(2, co_di{contrasts<0,d});
    si_di{2,d} = cat(2, co_di{contrasts==0,d});
    si_di{3,d} = cat(2, co_di{contrasts>0,d});
end

%split trials as CONTRAST x OUTCOME (nContrasts x nOutcomes)       
for o = 1:length(outcomes)
    for c = 1:length(contrasts)
        if strcmp(movementTime,'>0.5')
            [~, co_ou{c,o}] = selectCondition(expInfo, contrasts(c), behavioralData, initTrialConditions('movementTime','all','specificRTs',[.5 Inf],'responseType',outcomes{o}));
        else
            [~, co_ou{c,o}] = selectCondition(expInfo, contrasts(c), behavioralData, initTrialConditions('movementTime',movementTime,'responseType',outcomes{o}));
        end
    end
end
for d = 1:size(co_di,2)
    si_ou{1,d} = cat(2, co_ou{contrasts<0,d});
    si_ou{2,d} = cat(2, co_ou{contrasts==0,d});
    si_ou{3,d} = cat(2, co_ou{contrasts>0,d});
end

%split trials as CONTRAST x BLOCK (nContrasts x nBlocks)
for d = 1:length(blocks)
    for c = 1:length(contrasts)
        if strcmp(movementTime,'>0.5')
            [~, co_bl{c,d}] = selectCondition(expInfo, contrasts(c), behavioralData, initTrialConditions('responseType','correct','movementTime','all','specificRTs',[.5 Inf],'highRewardSide',blocks{d}));
        else
            [~, co_bl{c,d}] = selectCondition(expInfo, contrasts(c), behavioralData, initTrialConditions('responseType','correct','movementTime',movementTime,'highRewardSide',blocks{d}));
        end
    end
end
for d = 1:size(co_di,2)
    si_bl{1,d} = cat(2, co_bl{contrasts<0,d});
    si_bl{2,d} = cat(2, co_bl{contrasts==0,d});
    si_bl{3,d} = cat(2, co_bl{contrasts>0,d});
end

%split as OUTCOME x DIRECTION (counterbalanced)
for d = 1:size(di,1)
    ou_di{1,d} = intersect(di{d,:},find(~isnan(behavioralData.eventTimes(4).daqTime)));
    ou_di{2,d} = intersect(di{d,:},find(isnan(behavioralData.eventTimes(4).daqTime)));
end

%split trials as CONTRAST x DIRECTION x BLOCK (nContrasts x nDirs x nBlocks)    
for d = 1:length(dirs)
    for h = 1:length(blocks)
        for c = 1:length(contrasts)
            if strcmp(movementTime,'>0.5')
                [~, co_di_bl{c,d,h}] = selectCondition(expInfo, contrasts(c), behavioralData, initTrialConditions('movementTime','all','specificRTs',[.5 Inf],'movementDir',dirs{d},'highRewardSide',blocks{h}));
            else
                [~, co_di_bl{c,d,h}] = selectCondition(expInfo, contrasts(c), behavioralData, initTrialConditions('movementTime',movementTime,'movementDir',dirs{d},'highRewardSide',blocks{h}));
            end
        end
    end
end
for d = 1:size(co_di_bl,2)
    for h = 1:size(co_di_bl,3)
        si_di_bl{1,d,h} = cat(2, co_di_bl{contrasts<0,d,h});
        si_di_bl{2,d,h} = cat(2, co_di_bl{contrasts==0,d,h});
        si_di_bl{3,d,h} = cat(2, co_di_bl{contrasts>0,d,h});
    end
end

trialTypes.intVar.all = struct(...
    'contrast_direction',{co_di},...
    'contrast_block',{co_bl},...
    'contrast_outcome',{co_ou},...
    'contrast_direction_block',{co_di_bl},...
    'side_direction',{si_di},...
    'side_outcome',{si_ou},...
    'side_block',{si_bl},...
    'side_direction_block',{si_di_bl},...
    'outcome_direction',{ou_di}...
    );

%%%%%%% INTERACTIONS (3-WAY COUNTERBALANCED) %%%%%%%%%%%%%%%%%%%

%counterbalance CONTRAST x DIRECTION x BLOCK
sizec = cellfun('length', co_di_bl);
mincc = min(squeeze(min(sizec,[],2)),[],2);
for c = 1:size(co_di_bl,1)
    for d = 1:size(co_di_bl,2)
        for h = 1:size(co_di_bl,3)
            co_di_bl{c,d,h} = randsample(co_di_bl{c,d,h},mincc(c));
            if size(co_di_bl{c,d,h},1) == 0
                co_di_bl{c,d,h} = co_di_bl{c,d,h}';
            end
        end
    end
end

%split as SIDE x DIRECTION x BLOCK (counterbalanced)
for d = 1:size(co_di_bl,2)
    for h = 1:size(co_di_bl,3)
        si_di_bl{1,d,h} = cat(2, co_di_bl{contrasts<0,d,h});
        si_di_bl{2,d,h} = cat(2, co_di_bl{contrasts==0,d,h});
        si_di_bl{3,d,h} = cat(2, co_di_bl{contrasts>0,d,h});
    end
end

%project to CONTRAST x BLOCK (counterbalanced)
for c = 1:size(co_di_bl,1)
    co_bl_3{c,1} = cat(2,co_di_bl{c,:,1});
    co_bl_3{c,2} = cat(2,co_di_bl{c,:,2});
end

%split as CONTRAST ONLY (counterbalanced)
for c = 1:size(co_bl_3,1)
    co_3{c,1} = cat(2,co_bl_3{c,:});
end

%split as SIDE ONLY (counterbalanced)
si_3{1,1} = cat(2, co_3{contrasts<0});
si_3{2,1} = cat(2, co_3{contrasts==0});
si_3{3,1} = cat(2, co_3{contrasts>0});

%project to SIDE x BLOCK (counterbalanced)
for d = 1:size(co_bl_3,2)
    si_bl_3{1,d} = cat(2, co_bl_3{contrasts<0,d});
    si_bl_3{2,d} = cat(2, co_bl_3{contrasts==0,d});
    si_bl_3{3,d} = cat(2, co_bl_3{contrasts>0,d});
end

%project to CONTRAST x DIRECTION (counterbalanced)
for c = 1:length(contrasts)
    co_di_3{c,1} = cat(2,co_di_bl{c,1,:});
    co_di_3{c,2} = cat(2,co_di_bl{c,2,:});
end

%project to SIDE x DIRECTION (counterbalanced)
for d = 1:size(co_di_3,2)
    si_di_3{1,d} = cat(2, co_di_3{contrasts<0,d});
    si_di_3{2,d} = cat(2, co_di_3{contrasts==0,d});
    si_di_3{3,d} = cat(2, co_di_3{contrasts>0,d});
end

%split as DIRECTION ONLY (counterbalanced)
di_3{:,1} = cat(2,si_di_3{:,1});
di_3{:,2} = cat(2,si_di_3{:,2});

%split as OUTCOME x DIRECTION (counterbalanced)
for d = 1:size(di_3,2)
    ou_di_3{1,d} = intersect(di_3{:,d},find(~isnan(behavioralData.eventTimes(4).daqTime)));
    ou_di_3{2,d} = intersect(di_3{:,d},find(isnan(behavioralData.eventTimes(4).daqTime)));
end

%split as OUTCOME ONLY (counterbalanced)
for d = 1:size(ou_di_3,2)
    ou_3{d,1} = cat(2,ou_di_3{d,:});
end

%split as BLOCK ONLY (counterbalanced)
bl_3{:,1} = cat(2,si_bl_3{:,1});
bl_3{:,2} = cat(2,si_bl_3{:,2});

%split as OUTCOME x BLOCK (counterbalanced)
for d = 1:size(bl_3,2)
    ou_bl_3{1,d} = intersect(bl_3{:,d},find(~isnan(behavioralData.eventTimes(4).daqTime)));
    ou_bl_3{2,d} = intersect(bl_3{:,d},find(isnan(behavioralData.eventTimes(4).daqTime)));
end

%split trials as DIRECTION x BLOCK (2 x 2)
di_bl_3{1,1} = cat(2,co_di_bl{:,1,1});
di_bl_3{1,2} = cat(2,co_di_bl{:,1,2});
di_bl_3{2,1} = cat(2,co_di_bl{:,2,1});
di_bl_3{2,2} = cat(2,co_di_bl{:,2,2});

trialTypes.intVar.cb3D = struct(...
    'contrast_direction',{co_di_3},...
    'side_direction',{si_di_3},...
    'outcome_direction',{ou_di_3},...
    'direction_block',{di_bl_3},...
    'contrast_block',{co_bl_3},...
    'side_block',{si_bl_3},...
    'outcome_block',{ou_bl_3},...
    'contrast',{co_3},...
    'side',{si_3},...
    'direction',{di_3'},...
    'block',{bl_3'},...
    'outcome',{ou_3},...
    'contrast_direction_block',{co_di_bl},...
    'side_direction_block',{si_di_bl}...
    );

%%%%%%% INTERACTIONS (2-WAY COUNTERBALANCED) %%%%%%%%%%%%%%%%%%%

%counterbalance CONTRAST x DIRECTION
sizec = cellfun('length', co_di);
mincc = min(squeeze(min(sizec,[],2)),[],2);
for c = 1:size(co_di,1)
    for d = 1:size(co_di,2)
        for h = 1:size(co_di,3)
            co_di_2{c,d,h} = randsample(co_di{c,d,h},mincc(c));
        end
    end
    if size(co_di_2{c,1})-size(co_di_2{c,2}) ~= [0 0]
        co_di_2{c,2} = co_di_2{c,2}';
    end
        
end

%split as CONTRAST ONLY (counterbalanced)
for c = 1:size(co_di_2,1)
    co_2{c,1} = cat(2,co_di_2{c,:});
end

%split as SIDE ONLY (counterbalanced)
si_2{1,1} = cat(2, co_2{contrasts<0});
si_2{2,1} = cat(2, co_2{contrasts==0});
si_2{3,1} = cat(2, co_2{contrasts>0});

%split as SIDE x DIRECTION (counterbalanced)
for d = 1:size(co_di_2,2)
    si_di_2{1,d} = cat(2, co_di_2{contrasts<0,d});
    si_di_2{2,d} = cat(2, co_di_2{contrasts==0,d});
    si_di_2{3,d} = cat(2, co_di_2{contrasts>0,d});
end

%split as DIRECTION ONLY (counterbalanced)
di_2{:,1} = cat(2,si_di_2{:,1});
di_2{:,2} = cat(2,si_di_2{:,2});

%split as OUTCOME x DIRECTION (counterbalanced)
for d = 1:size(di_2,2)
    ou_di_2{1,d} = intersect(di_2{:,d},find(~isnan(behavioralData.eventTimes(4).daqTime)));
    ou_di_2{2,d} = intersect(di_2{:,d},find(isnan(behavioralData.eventTimes(4).daqTime)));
end

%split as OUTCOME ONLY (counterbalanced)
for d = 1:size(ou_di_2,2)
    ou_2{d,1} = cat(2,ou_di_2{d,:});
end

trialTypes.intVar.cb2D = struct(...
    'contrast_direction',{co_di_2},...
    'side_direction',{si_di_2},...
    'outcome_direction',{ou_di_2},...
    'contrast',{co_2},...
    'side',{si_2},...
    'direction',{di_2'},...
    'outcome',{ou_2}...
    );