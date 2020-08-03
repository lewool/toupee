function figList = psychometric(x, y, varargin)
% Creates a psychometric plot
%
%
% Inputs:
% -------
% x : array OR cell array
%   The 'x' values for the psychometric. If a cell array, each cell
%   represents the 'x' values for one curve to plot.
% y : array OR cell array
%   The response values for each element in `x` for the psychometric. 
%   Must be the same size as `x`.
% 'responseVal' : scalar numeric OR string (optional name-value pair)
%   The response value to use when creating the psychometric. (Default: 1)
% 'fitFn' : function handle (optional name-value pair)
%   The psychometric function to use as a fit to the data. (Default: none)
% 'figType' : char array (optional name-value pair)
%   The number of figures to return. Options:
%   'single' : a single figure containing one curve for each row in `x`
%   'multiple' : one figure each with a single curve for each row in `x`
%   (Default : 'single')
% 'cb' : double scalar OR char array
%   The confidence bound to plot around each unique value in `x`. Must be a
%   value (0,1) or 'none' to specify no bounds. (Default: 0.99)
% 'fig' : figure
%   Handle to figure to plot on top of. (Default: none)
% 'dataLabel' : logical
%   Whether to include a label for each datapoint in each psychometric
%   curve that gives the number of events that contribute to that point.
%   (Default : true)
%
%
% Ouputs:
% -------
% figList : figure array
%   Handle to created figure or array of figures.
%
%
% Examples:
% ---------
%
% 
% @todo add code and examples with psychometric fits
% @todo add optional plotting input args
%

%% Prerun checks.
% Imports.
import toupee.misc.*
% Validate inputs.
p = inputParser;
isValidX = @(z) all(size(z) == size(y));
isValidY = @(z) all(size(z) == size(x));
isValidResponseVal = @(z) isscalar(z);
isValidFitFn = @(z) isa(z, 'function_handle') || isempty(z);
isValidFigType = @(z) strcmpi(z, 'single') || strcmpi(z, 'multiple');
isValidCb = @(z) (isnumeric(z) && (z > 0 && z < 1)) || strcmpi(z, 'none');
isValidFig = @(z) isa(z, 'matlab.ui.Figure');

addRequired(p, 'x', isValidX);
addRequired(p, 'y', isValidY);
addParameter(p, 'responseVal', 1, isValidResponseVal);
addParameter(p, 'fitFn', [], isValidFitFn);
addParameter(p, 'figType', 'single', isValidFigType);
addParameter(p, 'cb', 0.99, isValidCb);
addParameter(p, 'fig', [], isValidFig);
addParameter(p, 'dataLabel', true, @islogical);

parse(p, x, y, varargin{:});
x = p.Results.x;
y = p.Results.y;
responseVal = p.Results.responseVal;
fitFn = p.Results.fitFn;
figType = p.Results.figType;
cb = p.Results.cb;
fig = p.Results.fig;
dataLabel = p.Results.dataLabel;

if strcmpi(cb, 'none'), cb = 0; end
if ~iscell(x), x = {x}; end; if ~iscell(y), y = {y}; end  % cellify

% - flags:
    	% 1) y-axis: p vs. z-score
    	% 2) fit or nofit
    	% 3) separate figures or one figure

%% Plot psychometrics
nC = numel(x);  % number of curves to plot
% Create figure if not supplied.
if isempty(fig)
    fig = figure;
else  % `figType` must be 'single' if `fig` was supplied
    if ~strcmpi(figType, 'single')
        errId = 'toupee:plot:psychometric:badFigType';
        errMsg = ['When providing a figure, `figType` must be set as ' ...
                  '''single'''];
        error(errId, errMsg);
    end
end
figList = iif(strcmpi(figType, 'single'), [fig], gobjects(nC, 1));
% Plot each curve.
for iC = 1:nC
    xC = x{iC};  % x vals for current curve
    yC = y{iC};  % y vals for current curve
    % Get unique x values.
    xvals = unique(xC)';
    xvals = xvals(:);
    nX = numel(xvals);  % number of unique xvals
    % For each unique x val, find the total number of associated y vals, 
    % the total number of hits for each associated y val, and the mean and 
    % ci bounds of the associated y vals.
    yCTot = zeros(nX, 1);
    yCHit = zeros(nX, 1);
    yCMean = zeros(nX, 1);
    yCBounds = zeros(nX, 2);
    for iX = 1:nX
        % y vals for current curve for current unique x val
        yCiX = yC(xC == xvals(iX));
        yCTot(iX) = numel(yCiX);
        yCHit(iX) = numel(find(yCiX == responseVal));
        [yCMean(iX), yCBounds(iX, :)] = ...
            binofit(yCHit(iX), yCTot(iX), (1 - cb));
    end
    % Plot curves.
    % Get figure to plot on (dependent on `figType`).
    try figure(figList(iC)), catch, figure(figList(1)), end; hold on
    % handle to psychometric graphic
    hp = plot(xvals, yCMean, 'Marker', 'o', 'LineWidth', 1.5);
    hp.MarkerFaceColor = hp.Color;
    % handle to bounds graphic
    if cb
        hb = fill([xvals; flipud(xvals)], [yCBounds(:,1); ...
                                           flipud(yCBounds(:,2))], ...
                  round(hp.Color, 1, 'significant'), ...
                  'edgecolor', 'none', 'facealpha', 0.2);
    end
    % Add labels to datapoints if specified.
    if dataLabel
        ht = text(xvals, yCMean, num2str(yCTot), 'FontSize', 9, ...
                  'VerticalAlignment', 'top', ...
                  'HorizontalAlignment', 'right');
    end
    % Assign `fig` into `figList` if returning multiple figures.
    if strcmpi(figType, 'multiple'), figList(iC) = fig; end
end


end