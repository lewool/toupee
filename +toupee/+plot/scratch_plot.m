eventNames = ...
    {{'interactiveOn', 'response'}};
eventWindows = [0, 0];
fs = 1000;
gradFn = @(x) gradient(movmean(x, 10));
wheelSpecs.radius = 0.031; wheelSpecs.res = 400; wheelSpecs.gain = 5;
[~, wheelMoves2] = ...
    toupee.behavioral.getWheelMoves(expInfo, 'eventNames', eventNames, ...
    'eventWindows', eventWindows, ...
    'fs', fs, 'gradFn', gradFn, ...
    'wheelSpecs', wheelSpecs);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nE = height(expInfo);  % number of experiment sessions
tT = sum(cellfun(@(x) numel(x), expInfo.behavioralData.earlyMovesTrials));  % total number of trial across all sessions

%% Can we predict preGoCue moves from preStimOn moves?

% create bar of preStimMovement times from experiment start,
% color-coded by those that are followed by a preGoCue movement, 
% color-coded by those that have movement in same direction (give
% both of these percentages as well) 

psmTimes = cell(nE, 1); % pre stim move
psmpgcmTimes = cell(nE, 1);  % pre stim move + pre go cue move
psmpgcmsdTimes = cell(nE, 1);  % last pre stim move + last pre go cue move same direction

for iE = 1:nE
    bd = expInfo{iE, 'behavioralData'};
    % for each experiment, find all trials with preStimOn moves.
    psmMask = bd.preStimOnMovesTrials{1};
    psmTimes{iE} = ...
        bd.daqDelay{1} ...
        + bd.eventTimes{1}{'stimulusOn', 'rigTimes'}{1}(psmMask);
    % for these trials, find those that also have a preGoCue move.
    % (get this %)
    psmpgcmMask = bd.preStimOnMovesTrials{1} & bd.preGoCueMovesTrials{1};
    psmpgcmTimes{iE} = ...
        bd.daqDelay{1} ...
        + bd.eventTimes{1}{'stimulusOn', 'rigTimes'}{1}(psmpgcmMask);
    % for these trials, find moves where dir(preStimOn) == dir(preGoCue).
    % (get this %)
    wm = bd.wheelMoves{1};
    psmDir = wm{'newTrial,stimulusOn: [0, 0]', 'moveDirection'}{1};
    % fill empty cells with '--'
    psmDir(cellfun(@(x) isempty(x), psmDir)) = {{'--'}};
    psmLastDir = ...
        cellfun(@(x) ...
            cellfun(@(y) y((find(y == '-', 1, 'last') + 1):end), ...
                    x(end), 'uni', 0), psmDir, 'uni', 0);
    pgcDir = wm{'stimulusOn,interactiveOn: [0, 0]', 'moveDirection'}{1};
    % fill empty cells with '--'
    pgcDir(cellfun(@(x) isempty(x), pgcDir)) = {{'--h'}};
    pgcLastDir = ...
        cellfun(@(x) ...
            cellfun(@(y) y((find(y == '-', 1, 'last') + 1):end), ...
                    x(end), 'uni', 0), pgcDir, 'uni', 0);
    psmpgcmsdMask = cellfun(@(x, y) strcmp(x,y), psmLastDir, pgcLastDir);
    psmpgcmsdTimes{iE} = ...
         bd.daqDelay{1} ...
        + bd.eventTimes{1}{'stimulusOn', 'rigTimes'}{1}(psmpgcmsdMask);
end

psmTimesAll = sort(cell2mat(psmTimes));
psmpgcmTimesAll = sort(cell2mat(psmpgcmTimes));
psmpgcmsdTimesAll = sort(cell2mat(psmpgcmsdTimes));

figure, zz = scatter(psmTimesAll, rand(numel(psmTimesAll), 1), ...
                     'SizeData', 24);
hold on, zz2 = scatter(psmpgcmTimesAll, rand(numel(psmpgcmTimesAll), 1), ...
                       'SizeData', 24);
hold on, zz3 = scatter(psmpgcmsdTimesAll, ...
                       rand(numel(psmpgcmsdTimesAll), 1), 'SizeData', 24);

legend(...
    sprintf('preStim (%i / %i)', numel(psmTimesAll), tT), ...
    sprintf('preStim + preGoCue (%i / %i)', ...
            numel(psmpgcmTimesAll), numel(psmTimesAll)), ...
    sprintf('preStim + preGoCue + sameDir (%i / %i)', ...
            numel(psmpgcmsdTimesAll), numel(psmpgcmTimesAll)));
ylim([-1, 2])
xlabel('Time (s)')
title('Times of early movements', 'FontSize', 14)


% for preStimOn move trials, create histogram of nMoves (0, 1, 2, >=3)
% within period per trial, and overlay with % of those followed by a 
% preGoCue move
npsm = zeros(nE, 4);  % number of pre stim moves
psmgcm = zeros(nE, 4);  % pre stim move after go cue move

for iE = 1:nE
    bd = expInfo{iE, 'behavioralData'};
    wm = bd.wheelMoves{1};
    zeroIdxs = cellfun(@(x) x == 0, wm{'newTrial,stimulusOn: [0, 0]', ...
                       'nMoves'}{1});
    oneIdxs = cellfun(@(x) x == 1, wm{'newTrial,stimulusOn: [0, 0]', ...
                                      'nMoves'}{1});
    twoIdxs = cellfun(@(x) x == 2, wm{'newTrial,stimulusOn: [0, 0]', ...
                                      'nMoves'}{1});
    threeIdxs = cellfun(@(x) x >= 3, wm{'newTrial,stimulusOn: [0, 0]', ...
                                        'nMoves'}{1});
    npsm(iE, 1) = numel(find(zeroIdxs));
    npsm(iE, 2) = numel(find(oneIdxs));
    npsm(iE, 3) = numel(find(twoIdxs));
    npsm(iE, 4) = numel(find(threeIdxs));
    psmgcm(iE, 1) = numel(find(cellfun(@(x) x > 0, ...
        wm{'stimulusOn,interactiveOn: [0, 0]', 'nMoves'}{1}(zeroIdxs))));
    psmgcm(iE, 2) = numel(find(cellfun(@(x) x > 0, ...
        wm{'stimulusOn,interactiveOn: [0, 0]', 'nMoves'}{1}(oneIdxs))));
    psmgcm(iE, 3) = numel(find(cellfun(@(x) x > 0, ...
        wm{'stimulusOn,interactiveOn: [0, 0]', 'nMoves'}{1}(twoIdxs))));
    psmgcm(iE, 4) = numel(find(cellfun(@(x) x > 0, ...
        wm{'stimulusOn,interactiveOn: [0, 0]', 'nMoves'}{1}(threeIdxs))));
    
end

figure, bar([0 1 2 3], ...
            [sum(npsm(:,1)), sum(npsm(:,2)), ...
             sum(npsm(:,3)), sum(npsm(:,4)); ...
             sum(psmgcm(:,1)), sum(psmgcm(:,2)), ...
             sum(psmgcm(:,3)), sum(psmgcm(:,4))]);

set(ax, 'YGrid', 'on')
set(ax, 'YTick', [0:250:8000]) 
xlabel('nMoves')
ylabel('nTrials')
legend('preStim', 'preGoCue post preStim')
title('number of moves during early periods', 'FontSize', 14)
%% Can we predict intOn move direction from early period move direction?

% (look at early, preStimOn, and preGoCue periods)
% (evts.responseValues: -1 = right, 1 = left)
% for each session, get trials with early move
% for these trials, find preferred direction (sign of sum of wheel displacement)
% for these trials, find % of interactiveOn moves whose sign of sum of wheel displacement matches

psmL = zeros(nE, 1);  % left pre stim
psmR = zeros(nE, 1);  % right pre stim
ioPsmR = zeros(nE, 1);  % int on right following pre stim right
ioPsmL = zeros(nE, 1);  % int on left following pre stim left
correctPsmL = zeros(nE, 1);  % correct trials following pre stim left
correctPsmR = zeros(nE, 1);  % correct trials following pre stim right
correctIoPsmL = zeros(nE, 1);
correctIoPsmR = zeros(nE, 1);

pgcL = zeros(nE, 1);  % left pre go cue
pgcR = zeros(nE, 1); 
ioPgcR = zeros(nE, 1); 
ioPgcL = zeros(nE, 1);
correctPgcL = zeros(nE, 1);
correctPgcR = zeros(nE, 1);
correctIoPgcL = zeros(nE, 1);
correctIoPgcR = zeros(nE, 1);

earlyL = zeros(nE, 1);  % left early
earlyR = zeros(nE, 1); 
ioEarlyR = zeros(nE, 1); 
ioEarlyL = zeros(nE, 1);
correctEarlyL = zeros(nE, 1);
correctEarlyR = zeros(nE, 1);
correctIoEarlyL = zeros(nE, 1);
correctIoEarlyR = zeros(nE, 1);

for iE = 1:nE
    bd = expInfo{iE, 'behavioralData'};
    wm = bd.wheelMoves{1};
    fb = expInfo{iE, 'BlockFile'}{1}{1, 'events'}{1, 'feedbackValues'}{1}';   % feedback
    
    psmWdAll = cellfun(@(x) sum(x), wm{'newTrial,stimulusOn: [0, 0]', 'moveDisplacement'}{1});
    psmWd = psmWdAll(bd.preStimOnMovesTrials{1});
    psmL(iE) = numel(find(psmWd < 0));
    psmR(iE) = numel(find(psmWd > 0));
    ioPsmWdAll = cellfun(@(x) sum(x), wm{'interactiveOn,response: [0, 0]', 'moveDisplacement'}{1});
    ioPsmL(iE) = numel(find(ioPsmWdAll < 0 & psmWdAll < 0));
    ioPsmR(iE) = numel(find(ioPsmWdAll > 0 & psmWdAll > 0));
    if numel(fb) > numel(psmWdAll), fb = fb(1:(end - 1)); end
    correctPsmL(iE) = numel(find(fb & (psmWdAll < 0)));
    correctPsmR(iE) = numel(find(fb & (psmWdAll > 0)));
    correctIoPsmL(iE) = numel(find(fb & (ioPsmWdAll < 0) & (psmWdAll < 0)));
    correctIoPsmR(iE) = numel(find(fb & (ioPsmWdAll > 0) & (psmWdAll > 0)));
    
    pgcWdAll = cellfun(@(x) sum(x), wm{'stimulusOn,interactiveOn: [0, 0]', 'moveDisplacement'}{1});
    pgcWd = pgcWdAll(bd.preGoCueMovesTrials{1});
    pgcL(iE) = numel(find(pgcWd < 0));
    pgcR(iE) = numel(find(pgcWd > 0));
    ioPgcWdAll = cellfun(@(x) sum(x), wm{'interactiveOn,response: [0, 0]', 'moveDisplacement'}{1});
    ioPgcL(iE) = numel(find(ioPgcWdAll < 0 & pgcWdAll < 0));
    ioPgcR(iE) = numel(find(ioPgcWdAll > 0 & pgcWdAll > 0));
    correctPgcL(iE) = numel(find(fb & (pgcWdAll < 0)));
    correctPgcR(iE) = numel(find(fb & (pgcWdAll > 0)));
    correctIoPgcL(iE) = numel(find(fb & (ioPgcWdAll < 0) & (pgcWdAll < 0)));
    correctIoPgcR(iE) = numel(find(fb & (ioPgcWdAll > 0) & (pgcWdAll > 0)));
    
    earlyWdAll = cellfun(@(x) sum(x), wm{'newTrial,interactiveOn: [0, 0]', 'moveDisplacement'}{1});
    earlyWd = earlyWdAll(bd.earlyMovesTrials{1});
    earlyL(iE) = numel(find(earlyWd < 0));
    earlyR(iE) = numel(find(earlyWd > 0));
    ioEarlyWdAll = cellfun(@(x) sum(x), wm{'interactiveOn,response: [0, 0]', 'moveDisplacement'}{1});
    ioEarlyL(iE) = numel(find(ioEarlyWdAll < 0 & (earlyWdAll < 0)));
    ioEarlyR(iE) = numel(find(ioEarlyWdAll > 0 & (earlyWdAll > 0)));
    correctEarlyL(iE) = numel(find(fb & (earlyWdAll < 0)));
    correctEarlyR(iE) = numel(find(fb & (earlyWdAll > 0)));
    correctIoEarlyL(iE) = numel(find(fb & (ioEarlyWdAll < 0) & (earlyWdAll < 0)));
    correctIoEarlyR(iE) = numel(find(fb & (ioEarlyWdAll > 0) & (earlyWdAll > 0)));
    
end

figure, bar([1, 2, 3, 4, 5, 6], ...
            [sum(psmL), sum(psmR), sum(pgcL), ...
             sum(pgcR), sum(earlyL), sum(earlyR); ...
             sum(correctPsmL), sum(correctPsmR), sum(correctPgcL), ...
             sum(correctPgcR), sum(correctEarlyL), sum(correctEarlyR); ...
             sum(ioPsmL), sum(ioPsmR), sum(ioPgcL), ...
             sum(ioPgcR), sum(ioEarlyL), sum(ioEarlyR); ...
             sum(correctIoPsmL), sum(correctIoPsmR), sum(correctIoPgcL), ...
             sum(correctIoPgcR), sum(correctIoEarlyL), sum(correctIoEarlyR); ...
             ]);
set(gca, 'xtick', [1:6], 'xticklabel', {'preStimL', 'preStimR', 'preGoCueL', 'preGoCueR', 'earlyL', 'earlyR'})
set(gca, 'YGrid', 'on')
set(gca, 'YTick', [0:250:3500])
ylabel('nTrials')
legend('total', 'correct', 'intOn same dir', 'intOn same dir correct')
title('Predicting intOn direction and response from early moves', 'FontSize', 14)

%% What do "wiggle" moves represent?

% get wiggle moves
% get contrast vals for wiggle moves
% get direction of first departure for wiggle moves, and compare to
% direction of correct choice ('correctResponse' == -1 for right turns)

% make 3 types of wig moves: 
% 1) >= 3 moves, at least one direction change (og), < 0.5s between moves
% 2) >= 2 moves, at least one direction change, < 0.5s between moves
% 3) >= 2 moves, < 0.5s between moves

wigMoveContrasts = cell(nE, 1);
nWigMoves = zeros(nE, 1);
nFwgScd = zeros(nE, 1);  % first wig dir, same correct dir

for iE = 1:nE
    bd = expInfo{iE, 'behavioralData'};
    wm = bd.wheelMoves{1};
    % mask for number of moves
    wigMaskN1 = cellfun(@(x) x >= 3, ...
                        wm{'interactiveOn,response: [0, 0]', 'nMoves'}{1});
    % mask for direction change
    wigMaskD1 = ...
        cellfun(@(x) numel(unique(sign(cell2mat(x(:,2))))) == 2, ...
                wm{'interactiveOn,response: [0, 0]', 'moveDirection'}{1});
    % mask for max time b/w moves
    wigMaskT1 = cellfun(@(x, y) all((x - y) < 0.5), ...
                        wm{'interactiveOn,response: [0, 0]', 'moveOff'}{1}, ...
                        wm{'interactiveOn,response: [0, 0]', 'moveOn'}{1});
    % wiggle mask as combination of above
    wigMask = wigMaskN1 & wigMaskD1 & wigMaskT1;
    contrasts = expInfo{iE, 'BlockFile'}{1}{1, 'events'}{1, 'contrastValues'}{1}';
    if numel(contrasts) > numel(wigMask), contrasts = contrasts(1:(end - 1)); end
    wigMoveContrasts{iE} = contrasts(wigMask);
    firstWigDir = cellfun(@(x) sign(x(1)), wm{'interactiveOn,response: [0, 0]', 'moveDisplacement'}{1}(wigMask));
    cr = expInfo{iE, 'BlockFile'}{1}{1, 'events'}{1, 'correctResponseValues'}{1}';  % correct responses
    if numel(cr) > numel(wigMask), cr = cr(1:(end - 1)); end
    cr = -cr;  % flip sign
    nWigMoves(iE) = numel(find(wigMask));
    nFwgScd(iE) = numel(find(cr(wigMask) == firstWigDir));
end

wigMoveContrastsAll = cell2mat(wigMoveContrasts);
xContrasts = unique(wigMoveContrastsAll);

figure, histogram(wigMoveContrastsAll, [-1.05:0.025:1.05])


%% How does early movement change as a function of task experience?

pctEarlyMoves = zeros(nE, 1);

for iE = 1:nE
    bd = expInfo{iE, 'behavioralData'};
    pctEarlyMoves(iE) = numel(find(bd.earlyMovesTrials{1})) / numel(bd.earlyMovesTrials{1});
end

figure, h = plot(pctEarlyMoves, 'marker', 'x');
h.MarkerSize = 8;
h.LineWidth = 1;
title('proportion of early move trials as a function of task xp', 'fontsize', 13)
xlabel('day')

%% Plot psychometrics

evts1 = expInfo{'2020-02-03_1_LEW031', 'BlockFile'}{1}.events;
nT1 = numel(evts1.endTrialValues{1});
x{1} = evts1.contrastValues{1}(1:nT1);
y{1} = evts1.responseValues{1}(1:nT1);
evts2 = expInfo{'2020-02-04_2_LEW031', 'BlockFile'}{1}.events;
nT2 = numel(evts2.endTrialValues{1});
x{2} = evts1.contrastValues{1}(1:nT2);
y{2} = evts2.responseValues{1}(1:nT2);


allContrastsMask = xC == xvals;
[zrow, zcol] = find(allContrastsMask);
cnts = histcounts(zrow, numel(xvals)); 
