function combinedNeuralData = combineNeuralData(varargin)

%find data
expInfo = varargin{1};
behavioralData = varargin{2};
neuralData = varargin{3};

%find number of experiments/sessions
numExps = length(neuralData);

%find number of ETA events (typically 3, see alignResps.m)
numEvents = size(neuralData(1).eta.alignedResps,2);

%designate a holding cell for all exps x all ETA events
holdingCell = cell(numExps,numEvents);

%for a given ETA event...
for ev = 1:numEvents
    
    %...for each experiment...
    for ex = 1:numExps
    
        %find number of neurons in the exp
        numCells = size(neuralData(ex).eta.alignedResps{1},3);
        
        %initialize a holding variable big enough 
        holdingCell{ex,ev} = cell(numCells,1);
        
        %...for every neuron in that experiment
        for iCell = 1:numCells
            
            %...load the aligned trial-by-trial activity (trials x eventWindow double)
            holdingCell{ex,ev}{iCell} = neuralData(ex).eta.alignedResps{ev}(:,:,iCell);
            
            %...and record which exp it was from (needed to identify trial
            %types in, eg., selectCondition.m
            expID{ex,:}(iCell,:) = ex;
            
        end
        
        % collect the cell stats for that experiment (previously computed)
        holdingPs{ex,:} =  neuralData(ex).stats.pValues;
        holdingHs{ex,:} =  neuralData(ex).stats.bfcH;
        
    end
    
end

% for each ETA event, concatenate neurons from all experiments 
% output: cell array of size TOTAL neurons x ETA events (each cell 
% element holds a single neuron's activity, size = trials x eventWindow)
for ev = 1:numEvents
    bigCell(:,ev) = cat(1,holdingCell{:,ev});
end

% concatenate stats and expID (dim1 = total neurons)
bigPValues = cat(1,holdingPs{:});
bigH = cat(1,holdingHs{:});
bigID = cat(1,expID{:});

%collect everything into a new structure
combinedNeuralData.eta = struct(...
    'eachCellAligned',{bigCell},...
    'eventWindow',{neuralData(1).eta.eventWindow},...
    'events',{neuralData(1).eta.events}...
    );

combinedNeuralData.info = struct(...
    'expID',bigID...
    );

combinedNeuralData.stats = struct(...
    'PValues',{bigPValues},...
    'bfcH',bigH,...
    'labels',{neuralData(1).stats.labels},...
    'bfcAlpha',{neuralData(1).stats.bfcAlpha}...
    );

%additional steps if matched sessions
if expInfo(1).cellMatched == 1
    
    %the number of cells in a matched experiment will be the same 
    %across all experiments/events, choose one for reference
    numCells = size(neuralData(1).eta.alignedResps{1},3);
    ar = cell(1,numEvents);
    for ev = 1:numEvents
        catExps = [];
        for ex = 1:numExps
            catExps = cat(1,catExps,neuralData(ex).eta.alignedResps{ev});
            ar{ev} = catExps; % corresponds to neuralData(ex).eta.alignedResps
        end
    end
    
    % generate eachAlignedCell, but for combined trials
    for ev = 1:numEvents
        for iCell = 1:numCells
            
            hh{iCell,ev} = ar{ev}(:,:,iCell);
            
        end
    end

    %collect into a substruct 'matched'
    combinedNeuralData.matched.eta = struct(...
        'alignedResps',{ar},...
        'eachCellAligned',{hh},...
        'eventWindow',{neuralData(1).eta.eventWindow},...
        'events',{neuralData(1).eta.events}...
        );

    %recompute statistics using ALL trials across ALL exps
    combinedNeuralData = getSignificantActivity(expInfo, behavioralData, combinedNeuralData);

end

