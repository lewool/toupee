function plotBlockShifts(varargin)
% 15 Nov 2019: LEW adapted to take inputs from either 2AFC or B2AFC experiments

% This function takes a list of mice/exps and generates a psychometric
% curve based on all the concatenated trials. You can feed one mouse and
% multiple exps, or multiple mice and their associated exps (lengths must
% match)

%INPUT:
% 1. plotPsychometric({{'LEW025'}}, {{'2019-11-15',1}});
% 2. plotPsychometric(mouseList, expList), where
%       mouseList = {{'Mouse1'}} or {{'Mouse1'},{'Mouse2'},{MouseN'}} 
%       expList = {{'2018-06-10',2,[2 3]},{'2019-03-27',1,[1]}}
% 3. plotPsychometric(expInfo), where expInfo is a struct of one/many
%       experiments (see initExpInfo.m)


%%%% get data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

allContrastValues = [];
allRewardSideValues = [];
allRepeatValues = [];
allHitValues = [];
allChoiceValues = [];
rs = [];

if nargin == 1
    for s = 1:length(varargin{1})
        expInfo = varargin{1}(s);
        expInfo = data.loadExpData(expInfo);
        block = expInfo.block;
        
        contrastValues = block.events.contrastValues;
        trialLimit = length(block.events.endTrialValues);
        allContrastValues = [allContrastValues, contrastValues(1:trialLimit)];
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

        allRepeatValues = [allRepeatValues, block.events.repeatNumValues(1:trialLimit)];
        allHitValues = [allHitValues, block.events.feedbackValues(1:trialLimit)];
        allChoiceValues = [allChoiceValues, block.events.responseValues(1:trialLimit)];
    end

elseif nargin > 1
    for d = 1:length(varargin{2})
        expList = varargin{2};
        mouseList = varargin{1};
        if length(mouseList) == 1
            mouseList = repelem(mouseList,length(expList));
        end

        expInfo.mouseName = mouseList{d}{:};
        expInfo.expDate = expList{d}{1};
        expInfo.expNum = expList{d}{2};

        expInfo = data.loadExpData(expInfo);
        block = expInfo.block;

        contrastValues = block.events.contrastValues;
        trialLimit = length(block.events.endTrialValues);
        allContrastValues = [allContrastValues, contrastValues(1:trialLimit)];
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

        allRepeatValues = [allRepeatValues, block.events.repeatNumValues(1:trialLimit)];
        allHitValues = [allHitValues, block.events.feedbackValues(1:trialLimit)];
        allChoiceValues = [allChoiceValues, block.events.responseValues(1:trialLimit)];

    end
end

difficultTrials = find((abs(allContrastValues) < 1.1));
difficultChoices = allChoiceValues(difficultTrials) == 1;

colors = [0.1 0.7 0.1; 1 .6 0; 0 .4 1; 1 0 0];

figure;
hold on;

blockTrials = expInfo.block.events.highRewardSideValues;
blockBounds = [1 find(diff(expInfo.block.events.blockSwitchesValues) ~= 0) max(difficultTrials)];
for b = 2:length(blockBounds)
    if mean(blockTrials(blockBounds(b-1):blockBounds(b))) > 0
        color = colors(2,:);
    elseif mean(blockTrials(blockBounds(b-1):blockBounds(b))) < 0
        color = colors(1,:);
    end
    plotBlock = fill(...
                    [0 0 1 1],...
                    [blockBounds(b-1) blockBounds(b) blockBounds(b) blockBounds(b-1)],...
                    color,'LineStyle','none');
    alpha(0.3)
end

plot(smooth(difficultChoices,20), difficultTrials, 'LineWidth',1,'Color', [.25 .25 .25])

ax = gca;
ax.TickDir = 'out';
ax.YDir = 'reverse';
% ax.Position = [1240 560 170 420];
axis([0 1 1 max(difficultTrials)]);

ylabel('Trial number')
xlabel('P(right choice)')
