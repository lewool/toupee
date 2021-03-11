function plotTimingByContrast(varargin)

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
    
    % 1 if late, 0 if early
	moveTypes = ~isnan(wm.epochs(3).onsetTimes) & isnan(wm.epochs(2).onsetTimes);

    previousChoiceValues = [0 block.events.responseValues(1:trialLimit-1)];


    allContrastValues = [allContrastValues, contrastValues];
    allPreviousContrastValues = [allPreviousContrastValues, previousContrastValues];
    allRTs = [allRTs, RTs]; 
    allMoveTypes = [allMoveTypes, moveTypes];
    allPreviousChoiceValues = [allPreviousChoiceValues, previousChoiceValues];


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

cc = unique(allContrastValues);
pp = nan(size(cc));
nn = nan(size(cc));
for cont = 1:length(cc)
    contIdxs = [];

    %determine which trials had contrast = cc(cont)
    contIdxs = allContrastValues == cc(cont) & allRepeatValues == 1 & allRTs == 1;

    %number of trials with contrast = cc(cont)
    nn(cont) = sum(contIdxs); 
    try
        pp(cont) = sum(allMoveTypes(contIdxs) == 0)/sum(contIdxs);
    catch
        if nn(cont) == 0
             pp(cont) = 0;
        end
    end 
end

%this should only happen for CW but just in case
if all(cc <= 1)
        cc = cc*100;
end

figure;

if size(cc,2) >= 2

    % fit and plot a psychometric curve (asymmetric lapse rate)
    parstart = [mean(cc), 3, 0.05, 0.05 ]; %threshold, slope, gamma1, gamma2
    parmin = [min(cc) 0 0 0];
    parmax = [max(cc) 30 0.40 0.40];
    c = -max(cc):max(cc);
    xlabel('Contrast (%)')
    ylabel('P(impulsive choice)')
    text(75,0,strcat(num2str(length(allChoiceValues)),{' trials'}));
    hold on;

    % compute CIs
    alpha = 0.1; % 90% CIs - optional argument? 
    [~, pci] = binofit(round(pp.*nn), nn, alpha); % get CIs of the binomial distribution
    errorbar(cc, pp, pp-pci(:,1)', pp-pci(:,2)', ...
        'k', 'Color', lineColor,'MarkerFaceColor',lineColor,'MarkerEdgeColor','w', 'Marker','o','MarkerSize', 8, ...
        'LineWidth', 1,'capsize',0)
    hold on;
    box off
else
    disp('Error: Not enough contrasts to plot performance')

end
    xlim([-2 102])
    set(gca, 'XTick', [0, 12, 50, 100])
    set(gca, 'XTickLabels', {'0', '12', '50', '100'})
%     set(gca, 'YTick', [0, .5])
%     set(gca, 'YTickLabels', {'0', '0.5'})
    ax3 = gca;
    ax3.TickDir = 'out'; 
