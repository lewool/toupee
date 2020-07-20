function [moveOn, MoveOff, displacement, dir, class, peakVel, peakAcc] =...
    getMoves(x, t, varargin)
% Gets and classifies wheel moves
%
% Algorithm (pseudocode):
% >> for each sample
%     >> if displacement b/w (sample) : (`tThresh` samples) > `xThresh`
%           >> (sample) : (sample of max displacement) belong to a movement
% >> find all movement starts
% >> for all movements
%     >> merge consecutive, same-direction movements that are separated by 
%     less than `tMinGap` seconds
%     >> find more exact movement start by looking in future (from
%     predefined movement start) for first sample with diff greater than
%     `xOnThresh`
%     >> find more exact end of movement by looking in future (from
%     predefined movement end) for first sample with diff less than 
%     `xOffThresh` 
% >> get total displacement of each movement
% >> get direction of each movement
% >> get duration of each movement
% >> classify movement ("flinch" if movement duration < `tThresh`)
%
% 
% Inputs:
% -------
% x : double array
%   Contains the wheel position (in m) at each sample.
% 
% t : double array
%   Contains the timestamp (in s) of each sample.
% 
% 'fs' : int scalar (optional name-value pair)
%   The sampling frequency (in hz) of the wheel. (Default: 1000 hz)
%
% 'xThresh' : double scalar (optional name-value pair)
%   The minimum change in position (in m) of the wheel (within `tThresh`) 
%   for it to be classified as a movement. (Default: 0.001 m)
%
% 'tThresh' : double scalar (optional name-value pair)
%   The maximum duration (in s) in which the wheel position must change by
%   at least `xThresh` for it to be classified as a movement. 
%   (Default: 0.2 s)
%
% 'tMinGap' : double scalar (optional name-value pair)
%   The minimum duration (in s) between consecutive wheel movements for the
%   movements to be considered separate. (Default: 0.1 s)
% 
% 'xOnThresh' : double scalar (optional name-value pair)
%   The minimum change in position (in m) of the wheel for a more exact
%   wheel movement start time to be defined from an initial wheel movement
%   start time. (Default: 0.0005 m)
%
% 'xOffThresh' : double scalar (optional name-value pair)
%   The maximum allowable change in position (in m) of the wheel for a more
%   exact wheel movement end time to be defined from an initial wheel
%   movement end time. (Default: 0.0005 m)
%
% 'minDur' : double scalar (optional name-value pair)
%   The minimum duration (in s) a wheel movement must take place over for
%   it not to be discarded. (Default: 0 s)
%
% 'makePlots' : logical (optional name-value pair)
%   A flag to indicate whether to plot the detected wheel movements.
%
%
% Outputs:
% --------
% moveOn : double array
%   The times (in s) each movement began.
%
% moveOff : double array
%   The times (in s) each movement ended.
%
% displacement : double array 
%   The change in position (in m) of each movement.
% 
% dir : cell array
%   The direction ('left' or 'right') of each movement.
%
% class : cell array
%   The class ('flinch' or 'smooth') of each movement.
%
% peakVel : double array
%   The peak velocity (in m/s) of each movement.
%
% peakAcc : double array
%   The peak acceleration (in m/s^2) of each movement.
%
%
% Examples:
% ---------
%

%% Prerun checks.
% Validate inputs.
p = inputParser;
isValidX = @(y) isnumeric(y) && isvector(y) && numel(y) == numel(t);
isValidT = @(y) isnumeric(y) && isvector(y) && numel(y) == numel(x);
isValidParam = @(y) isnumeric(y) && isscalar(y) && (y > 0);

addRequired(p, 'x', isValidX);
addRequired(p, 't', isValidT);
addParameter(p, 'fs', 1000, isValidParam);
addParameter(p, 'xThresh', 0.001, isValidParam);
addParameter(p, 'tThresh', 0.2, isValidParam);
addParameter(p, 'tMinGap', 0.05, isValidParam);
addParameter(p, 'xOnThresh', 0.0005, isValidParam);
addParameter(p, 'xOffThresh', 0.0005, isValidParam);
addParameter(p, 'minDur', 0, isValidParam);
addParameter(p, 'makePlots', false, @(x) islogical(x) && isscalar(x));

parse(p, x, t, varargin{:});
p = p.Results;  % final parameters

%% Compute approximate movement onset and offset samples
% Convert the time threshold for detecting movements into a number of 
% samples threshold (given the sampling frequency)
sThresh = round(p.tThresh * p.fs);
nS = numel(t);  % total number of samples
% For each sample, see if it belongs to a movement.
moveMask = false(nS, 1);
for iS = 1:(nS - sThresh)
    % Find all current samples that pass thresh for movement.
    % displacement for current samples
    curDis = x(iS) - x(iS:(iS + sThresh));
    curMoveMask = abs(curDis) > p.xThresh;
    % If the samples pass thresh in both directions, just keep the move in
    % the first direction, the second will be caught later.
    % direction of displacement for current samples
    curDisDir = sign(diff(curDis(curMoveMask)));
    if numel(unique(curDisDir)) == 3
        iFirstMoveStop = find(curDisDir == -curDisDir(1), 1, 'first');
        iMoves = find(curMoveMask);
        curMoveMask(iMoves(iFirstMoveStop:end)) = 0;
    end
    % If a movement that passes thresh has position change in both
    % directions, don't count it; it will be caught later.
    iEndMove = find(curMoveMask, 1, 'last') - 1;
    % direction of position change for current move
    curMoveDir = sign(diff(x(iS:(iS + iEndMove))));
    if numel(find(curMoveDir == -1)) > 2 ...
       && numel(find(curMoveDir == 1)) > 2 
        continue
    % Else, keep the movement and mark these samples.
    else
        moveMask(iS:(iS + iEndMove)) = true;
    end
end
% Make sure final sample allows for a movement end.
moveMask(end) = false;


tic
warning('off', 'MATLAB:hankel:AntiDiagonalConflict')
displacement = zeros(nS, 1); % Initialize vector of displacement
iS = 0;  % index of current sample
while true
    iX = (1:p.batchSize) + iS;  % indices of `x` to process
    % Remove any overflowed indices for last batch.
    if iX(end) > length(t)
        iX(iX > length(t)) = []; 
    end
    w2e = hankel(x(iX), nan(1, sThresh));
    displacement(iX) = max(w2e, [], 2) - min(w2e, [], 2);
    iS = iS + p.batchSize - sThresh;
    if iX(end) == length(t), break, end
end

isMoving = displacement > p.xThresh;
isMoving(end) = false; % make sure we end on an offset
toc

% fill in small gaps - this is like a dilation/contraction
onsetSamps = find(~isMoving(1:end-1) & isMoving(2:end));
offsetSamps = find(isMoving(1:end-1) & ~isMoving(2:end));
tooShort = find((onsetSamps(2:end) - offsetSamps(1:end-1)) / p.fs < p.minGap);
for q = 1:length(tooShort)
    isMoving(offsetSamps(tooShort(q)):onsetSamps(tooShort(q)+1)) = true;
end

%% Compute precise movement onset and offset samples
% definition of move onset: 
% - starting from the end of the tThresh window, look back until you find
% one that's not different from the moveOnset, by some smaller threshold.
onsetSamps = find(~isMoving(1:end-1) & isMoving(2:end));

% Same as above, this is to replace the two lines below in a memory
% controlled way without looping on every sample
%  wheelToepOnsets = wheelToep(moveOnsetSamps,:);
%  wheelToepOnsetsDev = abs(bsxfun(@minus, wheelToepOnsets, wheelToepOnsets(:,1)));
wheelToepOnsetsDev = zeros(length(onsetSamps), sThresh);
iS = 0; cwt = 0;
while ~isempty(onsetSamps)
    iX = (1:p.batchSize) + iS;
    [icomm] = intersect( iX(1:end-sThresh-1), onsetSamps );
    [~, itpltz] = intersect( iX(1:end-sThresh-1), onsetSamps );
    iX(iX > length(t)) = [];
    if ~isempty(icomm)        
        w2e = hankel(x(iX), nan(1,sThresh));
        w2e = abs(bsxfun(@minus, w2e, w2e(:,1)));
        wheelToepOnsetsDev(cwt + (1:length(icomm)), :) =  w2e(itpltz, :);
        cwt = cwt + length(icomm);
    end
    iS = iS + p.batchSize - sThresh;
    if iX(end) >= onsetSamps(end), break, end
end
warning('on', 'MATLAB:hankel:AntiDiagonalConflict')

hasOnset = wheelToepOnsetsDev > p.xThreshOnset;
[a, b] = find(~fliplr(hasOnset));
onsetLags = sThresh-accumarray(a(:), b(:), [], @min);
onsetSamps = onsetSamps + onsetLags;
onsets = t(onsetSamps);

% we won't do the same thing for offsets, instead just take the actual end
% of the isMoving. This is because we're just not so concerned about being
% temporally precise with these. 
offsetSamps = find(isMoving(1:end-1) & ~isMoving(2:end));
offsets = t(offsetSamps);

moveDurs = offsets - onsets;
tooShort = moveDurs < p.minDur;
onsetSamps = onsetSamps(~tooShort);
onsets = onsets(~tooShort); 
offsetSamps = offsetSamps(~tooShort);
offsets = offsets(~tooShort); 

moveGaps = onsets(2:end) - offsets(1:end-1);
gapTooSmall = moveGaps < p.minGap;
% for these, drop the offending offset and onset, which effectively joins
% the two
if ~isempty(onsets)
    onsets = onsets([true ~gapTooSmall]); % always keep first onset
    onsetSamps = onsetSamps([true ~gapTooSmall]);
    offsets = offsets([~gapTooSmall true]); % always keep last offset
    offsetSamps = offsetSamps([~gapTooSmall true]); % always keep last offset
end
onsets = onsets(:); % return a column
offsets = offsets(:); % return a column
% Calculate displacement
displacement = x(offsetSamps) - x(onsetSamps);
displacement = displacement(:); % return a column
% Calculate peak velocity times and peak amplitudes
vel = conv(diff([0 x]), wheel.gausswin(10), 'same');
peakVelTimes = nan(size(onsets));
peakAmps = nan(size(onsets));
for m = 1:numel(onsets)
    thisV = abs(vel(onsetSamps(m):offsetSamps(m)));
    peakVelTimes(m) = onsets(m) + find(thisV == max(thisV), 1) / p.fs;
    % Get index of maximum absolute position relative to move onset
    [~,I] = max(abs(x(onsetSamps(m):offsetSamps(m)) - x(onsetSamps(m))));
    peakAmps(m) = x(onsetSamps(m) + I) - x(onsetSamps(m));
end

%% see how it looks
if p.makePlots
    figure('Name', 'Wheel movements'); 
    % Plot the wheel position
    ax1 = subplot(2,1,1);
    hold on; 
    on = plot(onsets, x(onsetSamps), 'go', 'DisplayName', 'onset');
    off = plot(offsets, x(offsetSamps), 'bo', 'DisplayName', 'offset');
    hold on; 
    inMove = logical(WithinRanges(t, [onsets offsets]));
    in = plot(t(inMove), x(inMove), 'r.', 'DisplayName', 'in movement');
    plot(t(~inMove), x(~inMove), 'k.');
    ylabel('position');
    legend([on off in], 'Location', 'SouthEast')
    
    % Plot the velocity trace
    ax2 = subplot(2,1,2);
    vel = wheel.computeVelAcc(x, t);
    hold on; 
    plot(onsets, vel(onsetSamps), 'go');
    plot(offsets, vel(offsetSamps), 'bo');
    plot(t(inMove), vel(inMove), 'r.');
    plot(t(~inMove), vel(~inMove), 'k.');
    ylabel('velocity');
    xlabel('time (sec)');
    
    linkaxes([ax1 ax2], 'x');
end

end