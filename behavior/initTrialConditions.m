function trialConditions = initTrialConditions(varargin)
% This function sets up the types of trials you want in'selectCondition.m'. 
% Use with no inputs to generate a default struct with all trials included. 
%
% Name/value pairs available:
% repeatType: 'random','baited' 
% movementDir: 'ccw', 'cw'
% movementTime: 'early','late
% highRewardSide: 'left','right'
% responseType: 'correct','incorrect'
% rewardOutcome: 'rewarded','unrewarded'
% pastStimulus: 'right','left','zero'
% pastMovementDir: 'ccw','cw'
% pastResponseType: 'correct','incorrect'
% trialsBack: [any integer]

% Call specific name/value pairs to change them; e.g.,:
% initTrialConditions('highRewardSide','left','responseType','correct')
    
% count arguments
nArgs = length(varargin);
if round(nArgs/2)~=nArgs/2
   error('initTrialConditions needs name/value pairs')
end

% set up default struct
trialConditions = struct(...
    'repeatType',{'all'},...
    'movementDir',{'all'},...
    'movementTime',{'all'},...
    'highRewardSide',{'all'},...
    'preStimMovement',{'all'},...
    'responseType',{'all'},...
    'rewardOutcome',{'all'},...
    'pastStimulus',{'all'},...
    'pastMovementDir',{'all'},...
    'pastMovementTime',{'all'},...
    'pastResponseType',{'all'},...
    'trialsBack',1,...
    'switchBlocks',{'all'},...
    'whichTrials',{'all'},...
    'specificRTs',{'all'});

structNames = fieldnames(trialConditions);

% replace defaults with values
for pair = reshape(varargin,2,[]) %# pair is {propName;propValue}
   inpName = pair{1}; 
   if any(strcmp(inpName,structNames))
      trialConditions.(inpName) = pair{2};
   else
      error('%s is not a recognized parameter name',inpName)
   end
end