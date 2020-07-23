

function posOut = correctCounterDiscont(pos)
% Corrects for rotary encoder counter value discontinuities.
%
% The rotary encoder values are represented by a counter whose bit
% representation is uint32. This means that if we get to either the left or
% right edge of the counter (either 0 or (2^32 - 1)), then a subsequent
% left or right move, respectively, will push the counter to (2^32 -1) or
% 0, respectively, causing a discontinuity == (2^32 - 1). This code fixes
% these discontinuities by essentially transforming the uint32
% representation into int16.
%
%
% Inputs:
% -------
% pos : double array
%   The raw counter values of the rotary encoder.
%
%
% Outputs:
% --------
% posOut : double array
%   The corrected counter values.
%
%
% Examples:
% ---------
% x = block.inputs.wheelValues;
% xCorrected = wheel.correctCounterDiscont(x);
%

% Check for counter flow discontinuities in both directions (rightwards:
% (2^32 - 1) -> 0; leftwards: 0 -> (2^32 - 1))
posDiff = diff(pos(:)); 
posDiff(posDiff > 2^31) = posDiff(posDiff > 2^31) - 2^32;
posDiff(posDiff < -2^31) = posDiff(posDiff < -2^31) + 2^32;
% Shift back to initial position
posOut = cumsum([0; posDiff]) + pos(1);
