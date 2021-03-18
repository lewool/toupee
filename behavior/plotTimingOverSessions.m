function plotTimingOverSessions(varargin)

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


tt = nan(size(range));
pp = nan(size(tt));
nn = nan(size(tt));

for s = 1:length(range)
    expInfo = varargin{1}(range(s));
    try %to copy existing block data already in expInfo, 
        block = expInfo.block;
    catch %otherwise load it
        expInfo = data.loadExpData(expInfo);
        block = expInfo.block;
    end
    %load behavior (needed for options)
    behavioralData = varargin{2}(range(s));
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
    choiceValues = wm.epochs(5).moveDir;

    tIdx = block.events.repeatNumValues(1:trialLimit) == 1 & RTs == 1;
    tt(s) = s;
    nn(s) = sum(tIdx);
    try
        pp(s) = sum(moveTypes(tIdx) == 0)/sum(tIdx);
    catch
        if nn(s) == 0
             pp(s) = 0;
        end
    end
end

lineColor = 'k';

figure;

xlabel('Session number (days)')
ylabel('P(impulsive choice)')
hold on;

% compute CIs
alpha = 0.1; % 90% CIs - optional argument? 
[~, pci] = binofit(round(pp.*nn), nn, alpha); % get CIs of the binomial distribution
errorbar(tt, pp, pp-pci(:,1)', pp-pci(:,2)', ...
    'k', 'Color', lineColor,'MarkerFaceColor',lineColor,'MarkerEdgeColor','w', 'Marker','o','MarkerSize', 8, ...
    'LineWidth', 1,'capsize',0)
xlim([.5 tt(end)+.5])
hold on;
box off
ax3 = gca;
ax3.TickDir = 'out';