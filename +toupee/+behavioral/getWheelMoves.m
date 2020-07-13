function [expInfo, wheelMoves] = getWheelMoves(expInfo, sessions)
% Gets per-trial wheel information
%
%
% Inputs:
% -------
% expInfo : table
%   A table containing relevant information and data variables (columns) 
%   for particular experiment sessions (rows).
%
% sessions : int scalar OR int array OR char array OR cell array (optional)
%   Specific sessions for which to filter trials from, instead of from all
%   sessions.
%
%
% Outputs:
% --------
% expInfo : table
%   The updated `expInfo`.
%
% wheelMoves : table
%   Contains continuous position, continuous velocity, peak velocity, 
%   final movement direction, movement classifications per binned time, and 
%   movement relative to specified signals events for each trial in a 
%   session. The number of cells is equal to the number of specified
%   sessions.
%
%
% Examples:
% ---------
%
%
% See Also:
% ---------
% toupee.behavioral.getTrials
% toupee.behavioral.getTrialBlocks
% toupee.behavioral.getEventTimes
%

%% Prerun checks.


end