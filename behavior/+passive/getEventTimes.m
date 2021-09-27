function allEventTimes = getEventTimes(expInfo, signalsNames)
% for signalsNames, use the exact fieldnames from the block file or you'll get an error
% e.g. signalsNames = {'stimulusOnTimes' 'interactiveOnTimes' 'stimulusOffTimes'};

% 2018 Nov 21: Now collects feedback times (as well as a rough estimate of when,
% on incorrect trials, the valve would have fired)
% 2019 May 09: Collects feedback times for biased-likelihood exps when
% valve didn't fire
% 2019 Nov 18: expInfo struct input arg
% add multi expInfos
% 2021 July 20: PASSIVE

for ex = 1:length(expInfo)
%% PULL OUT EXPINFO VARIABLES

block = expInfo(ex).block;
Timeline = expInfo(ex).Timeline;

% check signals names
if ~isfield(block.events,signalsNames{1})
    allEventTimes{ex} = NaN;
else


%% DETECT PHD FLIPS
% since the photodiode doesn't always register a WHITE or BLACK signal
% during wheel movement, but also GRAY, we can use departures from GRAY to 
% BLACK or GRAY to WHITE to determine when threshold was reached and the
% feedback signal generated.

phdRaw = Timeline.rawDAQData(:,strcmp({Timeline.hw.inputs.name},'photoDiode'));

% apply a boxcar filter to get rid of some trace jaggedness
windowSize = 3;
b = (1/windowSize)*ones(1,windowSize);
a = 1;
phdFilt = filter(b, a, phdRaw);

%find the top three clusters of photodiode signal values 
[~, flipvals] = kmeans(phdFilt, 3); %[BLACK WHITE GRAY]

%interpolate the photodiode signal to one of these three cluster values
phdInterp = interp1(flipvals, flipvals, phdFilt, 'nearest', 'extrap');

%note when the phd flips occurred by capturing the diffs
whenFlips = abs(diff(phdInterp)) > 0;

% get rid of the timestamps when there was no flip
flipTimes = Timeline.rawDAQTimestamps(2:end); %good enough for govt work
flipTimes(whenFlips < 1) = [];

%% MATCH UP BLOCK TIMES TO STIM ONSET TIMES
% determine the closest value of stimWindowUpdateTimes to
% each relevant block event, and then map this to a nearby flip of the
% photodiode

numCompleteTrials = numel(block.events.endTrialValues);

% find the closest stimWindowUpdate timepoint to every block event
for iName = 1:length(signalsNames)
    signalsUpdates{iName} = interp1(block.stimWindowUpdateTimes, block.stimWindowUpdateTimes, block.events.(signalsNames{iName}), 'nearest', 'extrap');
end

% external TL now needs an offset to bring it close to signals times.
% do this by aligning the first 'real' flip to the first stim presentation
% plus a little buffer so TL is not brought ahead
tlOffset = flipTimes(2) - signalsUpdates{1}(1) - 0.1;

% set up the corrected-timestamp struct
eventTimes = struct('event',signalsNames);

for i = 1:length(signalsUpdates)
    updateTimes = signalsUpdates{i}(1:numCompleteTrials);
    updateTimes_offset = signalsUpdates{i}(1:numCompleteTrials) + tlOffset;
    
    for j = 1:numel(updateTimes_offset)

        % choose a timepoint in signals
        timepoint = updateTimes_offset(j);

        % find the earliest phd flip after the relevant stimWindowUpdate
        % timepoint
        firstFlip = first(flipTimes(flipTimes > timepoint));
        
        % collect these values for the minute
        try
            realTimes(j) = firstFlip;
        catch
            disp(j)
            error('oops')
        end

    end
    
    % save signals event times wrt to stimWindowUpdate times...
    eventTimes(i).signalsTime = getfield(block.events,signalsNames{i});
    eventTimes(i).signalsTime = eventTimes(i).signalsTime(1:numCompleteTrials);
    eventTimes(i).updateTime = updateTimes;
    
    % ... and also phd flip times
    eventTimes(i).daqTime = realTimes;
    
    % tidy up
    updateTimes_offset = [];
    realTimes = [];
end

allEventTimes{ex} = eventTimes;
end
end


