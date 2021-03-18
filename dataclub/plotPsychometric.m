function plotPsychometric(varargin)
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

% options struct:
% expIdx = feed a range of sessions you want to plot (e.g., [3:10])
% choiceHistory = plot psycho based on choice made @t-1
% diffHistory = plot psycho based on difficulty (contrast hi/lo) @t-1
% movType = 'early', 'late', or 'all'
%%%% get data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

if isstruct(varargin{1}) %if expInfo
    if nargin == 1 %if no options struct
        for s = 1:length(varargin{1}) %basic load from block
            expInfo = varargin{1}(s);
            try %to copy existing block data already in expInfo, 
                block = expInfo.block;
            catch %or else load it
                expInfo = data.loadExpData(expInfo);
                block = expInfo.block;
            end

            trialLimit = length(block.events.endTrialValues);
            contrastValues = block.events.contrastValues(1:trialLimit);
            RTs = (block.events.feedbackTimes(1:trialLimit) - block.events.stimulusOnTimes(1:trialLimit)) < 3.2;
            moveTypes = true(1,trialLimit);

            allContrastValues = [allContrastValues, contrastValues];
            allRTs = [allRTs, RTs]; 
            allMoveTypes = [allMoveTypes, moveTypes];        

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
        
    elseif nargin == 3 %options struct
        options = varargin{3};
        moveType = options.movType;
        try 
            exps = options.expIdx;
        catch
            exps = 1:length(varargin{1});
        end
        for s = exps
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
            contrastValues = block.events.contrastValues(1:trialLimit);
            previousContrastValues = [NaN block.events.contrastValues(1:trialLimit-1)];
            RTs = et(1).daqTime - wm.epochs(5).onsetTimes < 3;
            
            if strcmp(moveType,'early')
                moveTypes = ~isnan(wm.epochs(2).onsetTimes);
            elseif strcmp(moveType,'late')
                moveTypes = ~isnan(wm.epochs(3).onsetTimes) & isnan(wm.epochs(2).onsetTimes);
            elseif strcmp(moveType,'all')
                moveTypes = true(1,trialLimit);
            else
                moveTypes = true(1,trialLimit);
            end
            
            if varargin{3}.choiceHistory == 1
                previousChoiceValues = [0 block.events.responseValues(1:trialLimit-1)];
                allPreviousChoiceValues = [allPreviousChoiceValues, previousChoiceValues];
            end
                

            allContrastValues = [allContrastValues, contrastValues];
            allPreviousContrastValues = [allPreviousContrastValues, previousContrastValues];
            allRTs = [allRTs, RTs]; 
            allMoveTypes = [allMoveTypes, moveTypes];

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
    end
            
elseif ~isstruct(varargin{1}) %if mouseList/expList
    if nargin == 2 %if no options struct
        for d = 1:length(varargin{2}) %basic load from block
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

            trialLimit = length(block.events.endTrialValues);
            contrastValues = block.events.contrastValues(1:trialLimit);
            RTs = (block.events.feedbackTimes(1:trialLimit) - block.events.stimulusOnTimes(1:trialLimit)) < 3.2;
            moveTypes = true(1,trialLimit);

            allContrastValues = [allContrastValues, contrastValues];
            allRTs = [allRTs, RTs]; 
            allMoveTypes = [allMoveTypes, moveTypes];
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
    elseif nargin == 4 %if options struct
        options = varargin{3};
        moveType = options.movType;
        try 
            exps = options.expIdx;
        catch
            exps = 1:length(varargin{1});
        end
        for d = exps %basic load from block
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
            
            %load behavior (needed for options)
            behavioralData = varargin{2}(s);
            et = behavioralData.eventTimes;
            wm = behavioralData.wheelMoves;

            trialLimit = length(block.events.endTrialValues);
            contrastValues = block.events.contrastValues(1:trialLimit);
            previousContrastValues = [NaN block.events.contrastValues(1:trialLimit-1)];
            RTs = (block.events.feedbackTimes(1:trialLimit) - block.events.stimulusOnTimes(1:trialLimit)) < 3.2;
            
            if strcmp(moveType,'early')
                moveTypes = ~isnan(wm.epochs(2).onsetTimes);
            elseif strcmp(moveType,'late')
                moveTypes = ~isnan(wm.epochs(3).onsetTimes) & isnan(wm.epochs(2).onsetTimes);
            elseif strcmp(moveType,'all')
                moveTypes = true(1,trialLimit);
            else
                moveTypes = true(1,trialLimit);
            end
            
            if varargin{3}.choiceHistory == 1
                previousChoiceValues = [0 block.events.responseValues(1:trialLimit-1)];
                allPreviousChoiceValues = [allPreviousChoiceValues, previousChoiceValues];
            end

            allContrastValues = [allContrastValues, contrastValues];
            allPreviousContrastValues = [allPreviousContrastValues, previousContrastValues];
            allRTs = [allRTs, RTs]; 
            allMoveTypes = [allMoveTypes, moveTypes];
        end
    end
end

%%%% plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

colors = [0.1 0.7 0.1; 1 .6 0; 0 .4 1; 1 0 0];
figure;
hold on;
line([0 0],[0 1],'LineStyle',':','Color',[.5 .5 .5])
line([-100 100],[.5 .5],'LineStyle',':','Color',[.5 .5 .5])

% 2AFC EXPERIMENT OR ONE-SIDED B2AFC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if length(unique(allRewardSideValues)) < 2
    if mean(allRewardSideValues) > 0
        lineColor = [1 0 0];
    elseif mean(allRewardSideValues) < 0
        lineColor = [0 .4 1];
    else
        lineColor = [0 0 0];
    end
    cc = unique(allContrastValues);
    pp = nan(size(cc));
    nn = nan(size(cc));
    for cont = 1:length(cc)
        contIdxs = [];

        %determine which trials had contrast = cc(cont)
        contIdxs = allContrastValues == cc(cont) & allRepeatValues == 1 & allRTs == 1 & allMoveTypes == 1;

        %number of trials with contrast = cc(cont)
        nn(cont) = sum(contIdxs); 
        try
            pp(cont) = sum(allChoiceValues(contIdxs) == 1)/sum(contIdxs);
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

    % plot
    if size(cc,2) >= 2

        % fit and plot a psychometric curve (asymmetric lapse rate)
        parstart = [mean(cc), 3, 0.05, 0.05 ]; %threshold, slope, gamma1, gamma2
        parmin = [min(cc) 0 0 0];
        parmax = [max(cc) 30 0.40 0.40];
        [pars, L] = mle_fit_psycho([cc; nn; pp],'erf_psycho_2gammas', parstart, parmin, parmax);
        c = -max(cc):max(cc);
%         plot(c, erf_psycho_2gammas(pars, c), 'Color', lineColor, 'LineWidth', 1,'LineStyle','-')
        xlabel('Left Contrast (%)       Right Contrast (%)')
        ylabel('P(rightward choice)')
        text(75,0,strcat(num2str(length(allChoiceValues)),{' trials'}));
        hold on;

        % compute CIs
        alpha = 0.1; % 90% CIs - optional argument? 
        [~, pci] = binofit(round(pp.*nn), nn, alpha); % get CIs of the binomial distribution
        errorbar(cc, pp, pp-pci(:,1)', pp-pci(:,2)', ...
            'k', 'Color', lineColor,'MarkerFaceColor',lineColor,'MarkerEdgeColor','w', 'MarkerSize', 8, ...
            'LineWidth', 1,'capsize',0)
        hold on;
        box off
        else
            disp('Error: Not enough contrasts to plot performance')

    end
    axis([cc(1)*1.05 cc(end)*1.05 -0.05 1.05])
    set(gca, 'XTick', [-100, -50, -25, 0, 25, 50, 100])
    set(gca, 'XTickLabels', {'100', '50', '25', '0', '25', '50', '100'})
    set(gca, 'YTick', [0, 0.5, 1])
    set(gca, 'YTickLabels', {'0', '0.5', '1.0'})
    ax3 = gca;
    ax3.TickDir = 'out'; 

% B2AFC EXPERIMENT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif length(unique(allRewardSideValues)) == 2 && ~options.choiceHistory
    % plot trials with high reward on the left
    lineColor = colors(1,:);
    cc = unique(allContrastValues);
    pp = nan(size(cc));
    nn = nan(size(cc));
    for cont = 1:length(cc)
        contIdxs = [];

        %determine which trials had contrast = cc(cont)
        contIdxs = (allRewardSideValues < 0) & allContrastValues == cc(cont) & allRepeatValues == 1 & allRTs == 1 & allMoveTypes == 1;

        %number of trials with contrast = cc(cont)
        nn(cont) = sum(contIdxs); 
        try
            pp(cont) = sum(allChoiceValues(contIdxs) == 1)/sum(contIdxs);
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

    % plot
    if size(cc,2) >= 2

        % fit and plot a psychometric curve (asymmetric lapse rate)
        parstart = [mean(cc), 3, 0.05, 0.05 ]; %threshold, slope, gamma1, gamma2
        parmin = [min(cc) 0 0 0];
        parmax = [max(cc) 30 0.40 0.40];
        [pars, L] = mle_fit_psycho([cc; nn; pp],'erf_psycho_2gammas', parstart, parmin, parmax);
        c = -max(cc):max(cc);
%         plot(c, erf_psycho_2gammas(pars, c), 'Color', lineColor, 'LineWidth', 1,'LineStyle','-')
        hold on;

        % compute CIs
        alpha = 0.1; % 90% CIs - optional argument? 
        [~, pci] = binofit(round(pp.*nn), nn, alpha); % get CIs of the binomial distribution
        errorbar(cc, pp, pp-pci(:,1)', pp-pci(:,2)', ...
            'Marker','o', 'Color', lineColor,'MarkerFaceColor',lineColor,'MarkerEdgeColor','w', 'MarkerSize', 8, ...
            'LineWidth', 1,'capsize',0)
        
        hold on;

        box off
        else
            disp('Error: Not enough contrasts to plot performance')

    end

    % plot trials with high reward on the right
    lineColor = colors(2,:);
    cc = unique(allContrastValues);
    pp = nan(size(cc));
    nn = nan(size(cc));
    for cont = 1:length(cc)
        contIdxs = [];

        %determine which trials had contrast = cc(cont)
        contIdxs = (allRewardSideValues > 0) & allContrastValues == cc(cont) & allRepeatValues == 1 & allRTs == 1 & allMoveTypes == 1;

        %number of trials with contrast = cc(cont)
        nn(cont) = sum(contIdxs); 
        try
            pp(cont) = sum(allChoiceValues(contIdxs) == 1)/sum(contIdxs);
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

    % plot
    if size(cc,2) >= 2

        % fit and plot a psychometric curve (asymmetric lapse rate)
        parstart = [mean(cc), 3, 0.05, 0.05 ]; %threshold, slope, gamma1, gamma2
        parmin = [min(cc) 0 0 0];
        parmax = [max(cc) 30 0.40 0.40];
        [pars, L] = mle_fit_psycho([cc; nn; pp],'erf_psycho_2gammas', parstart, parmin, parmax);
        c = -max(cc):max(cc);
%         plot(c, erf_psycho_2gammas(pars, c), 'Color', lineColor, 'LineWidth', 1,'LineStyle','-')

        % compute CIs
        alpha = 0.1; % 90% CIs - optional argument? 
        [~, pci] = binofit(round(pp.*nn), nn, alpha); % get CIs of the binomial distribution
        errorbar(cc, pp, pp-pci(:,1)', pp-pci(:,2)', ...
            'Marker','o', 'Color', lineColor,'MarkerFaceColor',lineColor,'MarkerEdgeColor','w', 'MarkerSize', 8, ...
            'LineWidth', 1,'capsize',0)
        box off
        
    else
        disp('Error: Not enough contrasts to plot performance')
    end
    
    axis([cc(1)*1.05 cc(end)*1.05 -0.05 1.05])
    set(gca, 'XTick', [-100, -50, -25, 0, 25, 50, 100])
    set(gca, 'XTickLabels', {'100', '50', '25', '0', '25', '50', '100'})
    set(gca, 'YTick', [0, 0.5, 1])
    set(gca, 'YTickLabels', {'0', '0.5', '1.0'})        
    xlabel('Left contrast (%)             Right contrast (%)')
    ylabel('P(right choice)')
    text(75,0,strcat(num2str(min([sum(allRTs) sum(allMoveTypes)])),{' trials'}));

    ax3 = gca;
    ax3.TickDir = 'out';

% CHOICE HISTORY ONLY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif options.choiceHistory && ~options.diffHistory
    
    % plot trials where previous trial was a left choice
    lineColor = colors(3,:);
    cc = unique(allContrastValues);
    pp = nan(size(cc));
    nn = nan(size(cc));
    for cont = 1:length(cc)
        contIdxs = [];

        %determine which trials had contrast = cc(cont)
        contIdxs = allContrastValues == cc(cont) & allRepeatValues == 1 & allRTs == 1 & allMoveTypes == 1 & allPreviousChoiceValues < 0;

        %number of trials with contrast = cc(cont)
        nn(cont) = sum(contIdxs); 
        try
            pp(cont) = sum(allChoiceValues(contIdxs) == 1)/sum(contIdxs);
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

    % plot
    if size(cc,2) >= 2

        % fit and plot a psychometric curve (asymmetric lapse rate)
        parstart = [mean(cc), 3, 0.05, 0.05 ]; %threshold, slope, gamma1, gamma2
        parmin = [min(cc) 0 0 0];
        parmax = [max(cc) 30 0.40 0.40];
        [pars, L] = mle_fit_psycho([cc; nn; pp],'erf_psycho_2gammas', parstart, parmin, parmax);
        c = -max(cc):max(cc);
%         plot(c, erf_psycho_2gammas(pars, c), 'Color', lineColor, 'LineWidth', 1,'LineStyle','-')
        hold on;

        % compute CIs
        alpha = 0.1; % 90% CIs - optional argument? 
        [~, pci] = binofit(round(pp.*nn), nn, alpha); % get CIs of the binomial distribution
        errorbar(cc, pp, pp-pci(:,1)', pp-pci(:,2)', ...
            'Marker','o', 'Color', lineColor,'MarkerFaceColor',lineColor,'MarkerEdgeColor','w', 'MarkerSize', 8, ...
            'LineWidth', 1,'capsize',0)
        
        hold on;

        box off
        else
            disp('Error: Not enough contrasts to plot performance')

    end

    % plot trials where previous trial was a right choice
    lineColor = colors(4,:);
    cc = unique(allContrastValues);
    pp = nan(size(cc));
    nn = nan(size(cc));
    for cont = 1:length(cc)
        contIdxs = [];

        %determine which trials had contrast = cc(cont)
        contIdxs = allContrastValues == cc(cont) & allRepeatValues == 1 & allRTs == 1 & allMoveTypes == 1 & allPreviousChoiceValues > 0;

        %number of trials with contrast = cc(cont)
        nn(cont) = sum(contIdxs); 
        try
            pp(cont) = sum(allChoiceValues(contIdxs) == 1)/sum(contIdxs);
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

    % plot
    if size(cc,2) >= 2

        % fit and plot a psychometric curve (asymmetric lapse rate)
        parstart = [mean(cc), 3, 0.05, 0.05 ]; %threshold, slope, gamma1, gamma2
        parmin = [min(cc) 0 0 0];
        parmax = [max(cc) 30 0.40 0.40];
        [pars, L] = mle_fit_psycho([cc; nn; pp],'erf_psycho_2gammas', parstart, parmin, parmax);
        c = -max(cc):max(cc);
%         plot(c, erf_psycho_2gammas(pars, c), 'Color', lineColor, 'LineWidth', 1,'LineStyle','-')

        % compute CIs
        alpha = 0.1; % 90% CIs - optional argument? 
        [~, pci] = binofit(round(pp.*nn), nn, alpha); % get CIs of the binomial distribution
        errorbar(cc, pp, pp-pci(:,1)', pp-pci(:,2)', ...
            'Marker','o', 'Color', lineColor,'MarkerFaceColor',lineColor,'MarkerEdgeColor','w', 'MarkerSize', 8, ...
            'LineWidth', 1,'capsize',0)
        box off
        
    else
        disp('Error: Not enough contrasts to plot performance')
    end
    
    axis([cc(1)*1.05 cc(end)*1.05 -0.05 1.05])
    set(gca, 'XTick', [-100, -50, -25, 0, 25, 50, 100])
    set(gca, 'XTickLabels', {'100', '50', '25', '0', '25', '50', '100'})
    set(gca, 'YTick', [0, 0.5, 1])
    set(gca, 'YTickLabels', {'0', '0.5', '1.0'})        
    xlabel('Left contrast (%)             Right contrast (%)')
    ylabel('P(right choice)')
    text(75,0,strcat(num2str(min([sum(allRTs) sum(allMoveTypes)])),{' trials'}));

    ax3 = gca;
    ax3.TickDir = 'out';
    
% CHOICE HISTORY & DIFFICULTY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif options.choiceHistory && options.diffHistory
    
    % plot trials where previous trial was a left choice & difficult
    lineColor = [.6 .8 1];
    cc = unique(allContrastValues);
    pp = nan(size(cc));
    nn = nan(size(cc));
    for cont = 1:length(cc)
        contIdxs = [];

        %determine which trials had contrast = cc(cont)
        contIdxs = ...
            allContrastValues == cc(cont) & ...
            allRepeatValues == 1 & ...
            allRTs == 1 & ...
            allMoveTypes == 1 & ...
            allPreviousChoiceValues < 0 & ...
            abs(allPreviousContrastValues) < .25;

        %number of trials with contrast = cc(cont)
        nn(cont) = sum(contIdxs); 
        try
            pp(cont) = sum(allChoiceValues(contIdxs) == 1)/sum(contIdxs);
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

    % plot
    if size(cc,2) >= 2

        % fit and plot a psychometric curve (asymmetric lapse rate)
        parstart = [mean(cc), 3, 0.05, 0.05 ]; %threshold, slope, gamma1, gamma2
        parmin = [min(cc) 0 0 0];
        parmax = [max(cc) 30 0.40 0.40];
        [pars, L] = mle_fit_psycho([cc; nn; pp],'erf_psycho_2gammas', parstart, parmin, parmax);
        c = -max(cc):max(cc);
%         plot(c, erf_psycho_2gammas(pars, c), 'Color', lineColor, 'LineWidth', 1,'LineStyle','-')
        hold on;

        % compute CIs
        alpha = 0.1; % 90% CIs - optional argument? 
        [~, pci] = binofit(round(pp.*nn), nn, alpha); % get CIs of the binomial distribution
        errorbar(cc, pp, pp-pci(:,1)', pp-pci(:,2)', ...
            'Marker','o', 'Color', lineColor,'MarkerFaceColor',lineColor,'MarkerEdgeColor','w', 'MarkerSize', 8, ...
            'LineWidth', 1,'capsize',0)
        
        hold on;

        box off
        else
            disp('Error: Not enough contrasts to plot performance')

    end

    % plot trials where previous trial was a left choice & easy
    lineColor = [0 .4 1];
    cc = unique(allContrastValues);
    pp = nan(size(cc));
    nn = nan(size(cc));
    for cont = 1:length(cc)
        contIdxs = [];

        %determine which trials had contrast = cc(cont)
        contIdxs = ...
            allContrastValues == cc(cont) & ...
            allRepeatValues == 1 & ...
            allRTs == 1 & ...
            allMoveTypes == 1 & ...
            allPreviousChoiceValues < 0 & ...
            abs(allPreviousContrastValues) > .25;

        %number of trials with contrast = cc(cont)
        nn(cont) = sum(contIdxs); 
        try
            pp(cont) = sum(allChoiceValues(contIdxs) == 1)/sum(contIdxs);
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

    % plot
    if size(cc,2) >= 2

        % fit and plot a psychometric curve (asymmetric lapse rate)
        parstart = [mean(cc), 3, 0.05, 0.05 ]; %threshold, slope, gamma1, gamma2
        parmin = [min(cc) 0 0 0];
        parmax = [max(cc) 30 0.40 0.40];
        [pars, L] = mle_fit_psycho([cc; nn; pp],'erf_psycho_2gammas', parstart, parmin, parmax);
        c = -max(cc):max(cc);
%         plot(c, erf_psycho_2gammas(pars, c), 'Color', lineColor, 'LineWidth', 1,'LineStyle','-')

        % compute CIs
        alpha = 0.1; % 90% CIs - optional argument? 
        [~, pci] = binofit(round(pp.*nn), nn, alpha); % get CIs of the binomial distribution
        errorbar(cc, pp, pp-pci(:,1)', pp-pci(:,2)', ...
            'Marker','o', 'Color', lineColor,'MarkerFaceColor',lineColor,'MarkerEdgeColor','w', 'MarkerSize', 8, ...
            'LineWidth', 1,'capsize',0)
        box off
        
    else
        disp('Error: Not enough contrasts to plot performance')
    end
    
    % plot trials where previous trial was a right choice & easy
    lineColor = [1 0 0];
    cc = unique(allContrastValues);
    pp = nan(size(cc));
    nn = nan(size(cc));
    for cont = 1:length(cc)
        contIdxs = [];

        %determine which trials had contrast = cc(cont)
        contIdxs = ...
            allContrastValues == cc(cont) & ...
            allRepeatValues == 1 & ...
            allRTs == 1 & ...
            allMoveTypes == 1 & ...
            allPreviousChoiceValues > 0 & ...
            abs(allPreviousContrastValues) > .25;

        %number of trials with contrast = cc(cont)
        nn(cont) = sum(contIdxs); 
        try
            pp(cont) = sum(allChoiceValues(contIdxs) == 1)/sum(contIdxs);
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

    % plot
    if size(cc,2) >= 2

        % fit and plot a psychometric curve (asymmetric lapse rate)
        parstart = [mean(cc), 3, 0.05, 0.05 ]; %threshold, slope, gamma1, gamma2
        parmin = [min(cc) 0 0 0];
        parmax = [max(cc) 30 0.40 0.40];
        [pars, L] = mle_fit_psycho([cc; nn; pp],'erf_psycho_2gammas', parstart, parmin, parmax);
        c = -max(cc):max(cc);
%         plot(c, erf_psycho_2gammas(pars, c), 'Color', lineColor, 'LineWidth', 1,'LineStyle','-')

        % compute CIs
        alpha = 0.1; % 90% CIs - optional argument? 
        [~, pci] = binofit(round(pp.*nn), nn, alpha); % get CIs of the binomial distribution
        errorbar(cc, pp, pp-pci(:,1)', pp-pci(:,2)', ...
            'Marker','o', 'Color', lineColor,'MarkerFaceColor',lineColor,'MarkerEdgeColor','w', 'MarkerSize', 8, ...
            'LineWidth', 1,'capsize',0)
        box off
        
    else
        disp('Error: Not enough contrasts to plot performance')
    end
    
    % plot trials where previous trial was a right choice & difficult
    lineColor = [.8 .45 .45];
    cc = unique(allContrastValues);
    pp = nan(size(cc));
    nn = nan(size(cc));
    for cont = 1:length(cc)
        contIdxs = [];

        %determine which trials had contrast = cc(cont)
        contIdxs = ...
            allContrastValues == cc(cont) & ...
            allRepeatValues == 1 & ...
            allRTs == 1 & ...
            allMoveTypes == 1 & ...
            allPreviousChoiceValues > 0 & ...
            abs(allPreviousContrastValues) < .25;

        %number of trials with contrast = cc(cont)
        nn(cont) = sum(contIdxs); 
        try
            pp(cont) = sum(allChoiceValues(contIdxs) == 1)/sum(contIdxs);
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

    % plot
    if size(cc,2) >= 2

        % fit and plot a psychometric curve (asymmetric lapse rate)
        parstart = [mean(cc), 3, 0.05, 0.05 ]; %threshold, slope, gamma1, gamma2
        parmin = [min(cc) 0 0 0];
        parmax = [max(cc) 30 0.40 0.40];
        [pars, L] = mle_fit_psycho([cc; nn; pp],'erf_psycho_2gammas', parstart, parmin, parmax);
        c = -max(cc):max(cc);
%         plot(c, erf_psycho_2gammas(pars, c), 'Color', lineColor, 'LineWidth', 1,'LineStyle','-')

        % compute CIs
        alpha = 0.1; % 90% CIs - optional argument? 
        [~, pci] = binofit(round(pp.*nn), nn, alpha); % get CIs of the binomial distribution
        errorbar(cc, pp, pp-pci(:,1)', pp-pci(:,2)', ...
            'Marker','o', 'Color', lineColor,'MarkerFaceColor',lineColor,'MarkerEdgeColor','w', 'MarkerSize', 8, ...
            'LineWidth', 1,'capsize',0)
        box off
        
    else
        disp('Error: Not enough contrasts to plot performance')
    end
    
    axis([cc(1)*1.05 cc(end)*1.05 -0.05 1.05])
    set(gca, 'XTick', [-100, -50, -25, 0, 25, 50, 100])
    set(gca, 'XTickLabels', {'100', '50', '25', '0', '25', '50', '100'})
    set(gca, 'YTick', [0, 0.5, 1])
    set(gca, 'YTickLabels', {'0', '0.5', '1.0'})        
    xlabel('Left contrast (%)             Right contrast (%)')
    ylabel('P(right choice)')
    text(75,0,strcat(num2str(min([sum(allRTs) sum(allMoveTypes)])),{' trials'}));

    ax3 = gca;
    ax3.TickDir = 'out';
end

end


