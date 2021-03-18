function plotTimingOverTrials(varargin)

allContrastValues = [];
allRewardSideValues = [];
allRepeatValues = [];
allHitValues = [];
allChoiceValues = [];
allPreviousChoiceValues = [];
allPreviousContrastValues = [];
rs = [];
allRTs = [];
allMoveTypes = [];
allTrialNumValues = [];

if nargin > 2
    range = varargin{3}.expIdx;
else
    range = 1:length(varargin{1});
end
        
for s = range
    expInfo = varargin{1}(s);
    try %to copy existing block data already in expInfo, 
        block = expInfo.block;
    catch %otherwise load it
        expInfo = data.loadExpData(expInfo);
        block = expInfo.block;
    end
    %load behavior (needed for options)
    behavioralData = varargin{2}(s);
    et = behavioralData.eventTimes;
    wm = behavioralData.wheelMoves;

    trialLimit = length(block.events.endTrialValues);
    contrastValues = abs(block.events.contrastValues(1:trialLimit));
    previousContrastValues = [NaN block.events.contrastValues(1:trialLimit-1)];
    RTs = et(1).daqTime - wm.epochs(5).onsetTimes < 5;
    trialNumValues = block.events.trialNumValues(1:trialLimit);
    
    % 1 if late, 0 if early
	moveTypes = ~isnan(wm.epochs(3).onsetTimes) & isnan(wm.epochs(2).onsetTimes);

    previousChoiceValues = [0 block.events.responseValues(1:trialLimit-1)];


    allContrastValues = [allContrastValues, contrastValues];
    allPreviousContrastValues = [allPreviousContrastValues, previousContrastValues];
    allRTs = [allRTs, RTs]; 
    allMoveTypes = [allMoveTypes, moveTypes];
    allPreviousChoiceValues = [allPreviousChoiceValues, previousChoiceValues];
    allTrialNumValues = [allTrialNumValues, trialNumValues];

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
    allChoiceValues = [allChoiceValues, wm.epochs(5).moveDir];

end

lineColor = 'k';

ss = 0:50:1200;
pp = nan(size(ss));
nn = nan(size(ss));

for seg = 2:length(ss)
    segIdxs = [];
    
    %determine which trials fell within the given segment
    segIdxs = allTrialNumValues <= ss(seg) & allTrialNumValues > ss(seg-1) & allRepeatValues == 1 & allRTs == 1;
    
    %number of trials in that session segment
    nn(seg) = sum(segIdxs); 
    try
        pp(seg) = sum(allMoveTypes(segIdxs) == 0)/sum(segIdxs);
    catch
        if nn(seg) == 0
             pp(seg) = 0;
        end
    end
end

ss(isnan(pp)) = [];
nn(isnan(pp)) = [];
pp(isnan(pp)) = [];

figure;

xlabel('Session time (trials)')
ylabel('P(impulsive choice)')
hold on;

% compute CIs
alpha = 0.1; % 90% CIs - optional argument? 
[~, pci] = binofit(round(pp.*nn), nn, alpha); % get CIs of the binomial distribution
errorbar(ss, pp, pp-pci(:,1)', pp-pci(:,2)', ...
    'k', 'Color', lineColor,'MarkerFaceColor',lineColor,'MarkerEdgeColor','w', 'Marker','o','MarkerSize', 8, ...
    'LineWidth', 1,'capsize',0)
hold on;
box off
ax3 = gca;
    ax3.TickDir = 'out';