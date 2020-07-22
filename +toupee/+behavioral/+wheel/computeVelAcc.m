function [v, a, t] = computeVelAcc(x, t, varargin)
% Computes the velocity and acceleration of given position data
%
% Uses a user-specified custom function, `gradFn`, to compute the first and
% second derivatives of the position data. If `gradFn` is not supplied, the
% function used is the numerical gradient of a single-pass, moving-average,
% 10-sample window over the position.
%
% Returned units for velocity are in units of `x` per second.
% 
%
% Inputs:
% -------
% x : numeric array
%   Position data. Number of rows = number of datapoints, number of columns
%   = number of dimensions.
%
% t : numeric vector
%   Time data: the time associated with each sample.
%
% 'fs' : int scalar (optional name-value pair)
%   The sampling frequency of `x` in hertz.
%
% 'gradFn' : function handle (optional name-value pair)
%   A function handle containing the function to use on the position data
%   to compute the velocity. (Default: The numerical gradient of a 
%   single-pass, moving-average, 10-sample window over the position)
%
%
% Outputs:
% --------
% v : numeric array
%   The computed velocity at each timepoint in `t` for data in `x`.
%
% a : numeric array
%   The computed acceleration at each timepoint in `t` for data in `x`.
%
% t : numeric array
%   The linearly sampled (at `fs`) input arg, `t`.
%
% Examples:
% ---------
%
% 1) Compute velocity and acceleration using default gradient function.
%
% 2) Compute velocity and acceleration using custom gradient function.
%   gradFn = @(x) conv(diff(x), gausswin(10, 6), 'same');
%
% See Also:
% ---------
% toupee.behavioral.wheel.getWheelMoves
%

%% Prerun checks.
% Ensure input args are of proper type.
p = inputParser;
% must be numeric and the same size as `t`
isValidX = @(x) isnumeric(x) && all(size(x) == size(t));
% must be numeric and the same size as `x` and monotonically increasing
isValidT = @(t) isnumeric(t) && all(size(t) == size(x))...
                && all(diff(t) > 0);
% must be an int scalar
isValidFs = @(x) isscalar(x) && mod(x, 1) == 0;
% must be a function handle
isValidGradFn = @(x) isa(x, 'function_handle');

addRequired(p, 'x', isValidX);
addRequired(p, 't', isValidT);
addParameter(p, 'fs', 1000, isValidFs);
addParameter(p, 'gradFn', @(x) gradient(movmean(x, 10)), isValidGradFn);

parse(p, x, t, varargin{:});
x = p.Results.x;
t = p.Results.t;
fs = p.Results.fs;
gradFn = p.Results.gradFn;

%% Compute velocity and acceleration.
% Ensure `x` is sampled evenly in time.
if ~all(diff(t) == t(1) - t(2))
    tRaw = t;
    t = t(1):(1 / fs):t(end);
    x = interp1(tRaw, x, t(:), 'nearest');
end

% Compute derivatives.
v = gradFn(x) * fs;
a = gradFn(v) * fs;
    
end