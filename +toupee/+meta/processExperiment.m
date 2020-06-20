function [expInfo] = processExperiment(details, specs)
% Gets experiment information and data for given session(s)
% 
%
% Inputs:
% -------
% details : cell array
%   The session(s) details. If getting info for multiple sessions, use a 
%   nested cell array for each session. The elements in the innermost cells
%   should be subject name (char array), experiment date (char array in 
%   datestr format), and experiment number (int scalar). Optionally, there
%   can be a 4th element: series (int array), which is the specs for the
%   filenames saved by suite2P.
%
% specs : cell array
%   Specifications for which files to load into the returned `expInfo`
%   struct. If getting info for multiple sessions, use a nested cell array
%   for each session. The elements in the innermost cells can be:
%   1) 'block': loads the block file
%   2) 'timeline': loads the timeline file
%   3) individual behavioral or neural data files. These files can
%   correspond to block data (e.g. 'wheel.position.npy'), timeline data 
%   (e.g. 'rewardvalve.raw.npy'), or suite2P data (e.g. 'SVD_plane1.mat')
% 
%
% Outputs:
% --------
% expInfo : struct
%   A struct array with each element containing fields with information for
%   a particular session. The fields for each struct element can include:
%   'subject', 'expDate', 'expNum', 'expRef', 'expSeries', 'block',
%   'timeline', 'numPlanes', 'numChannels', 'behavioralData', 'neuralData'.
%
%
% Examples:
% ---------
% 1) Return the bare experiment info for a single session from a single
% subject.
%   details = {'LEW031', '2020-02-03', 1};
%   expInfo = toupee.meta.processExperiment(details);
%
% 2) Return experiment info and the block + timeline files for all of 
% multiple sessions from multiple subjects.
%   details = {{'LEW031', '2020-02-03', 1},... 
%              {'LEW037', '2020-03-13', 1},...
%              {'LEW005', '2018-06-10', 2, [2 3]}};
%   specs = {'block', 'timeline'};
%   expInfo = toupee.meta.processExperiment(details, specs);
%
% 3) Return experiment info and the block + specific individual data files
% for particular sessions.
%   details = {{'LEW031', '2020-02-03', 1},... 
%              {'LEW037', '2020-03-13', 1},...
%              {'LEW005', '2018-06-10', 2, [2 3]}};
%   specs = {{'timeline'}, {'block', 'rewardvalve.raw.npy'}, {'block'}};
%   expInfo = toupee.meta.processExperiment(details, specs);
%
%
% See Also:
% ---------
% `toupee.meta.constructExpRef`
% `toupee.meta.loadDatafile`
%

% Import all other functions in this subpackage and `iif`.
import toupee.meta.*
import toupee.misc.iif

% Do some checks on input args.
if ~iscell(details)  % ensure `details` is cell.
    error('toupee:meta:processExperiment:badInput',...
          'The "details" input arg should be a cell array')
elseif ~iscell(details{1})  % convert to nested cell if not already
    details = {details};
end

% For each experiment session, set up `expInfo`.
for e = 1:numel(details)
    % Initialize expInfo struct.
    d = details{e};
    expInfo(e) = struct(...
        'subject', d{1},...
        'expDate', d{2},...
        'expNum', d{3},...
        'expRef', constructExpRef(d{1}, d{2}, d{3}),...
        'expSeries', [],...
        'block', [],...
        'timeline', [],...
        'numPlanes', [],...
        'numChannels', [],...
        'behavioralData', struct(),...
        'neuralData', struct());  %#ok<*AGROW>
    expInfo(e).expSeries = iif(numel(d) > 3, @() d{4}, []);
    % Throw warning if there are extraneous details.
    if numel(d) > 4  
        warning('toupee:meta:processExperiment:extraDetails',...
            ['There should only be a maximum of 4 types of details, '...
             'but this session, "%s", was provided with all of the '...
             'following:'], expInfo(e).expRef);
         disp(d);
    end
end

if nargin > 1  % try to load data files
    if ~iscell(specs)  % ensure `specs` is cell
        error('toupee:meta:processExperiment:badInput',...
              'The "specs" input arg should be a cell array')
    elseif ~iscell(specs{1})  % convert to nested cell if not already
        specs = {specs};
    end
    expInfo = loadDatafile(expInfo, specs);  % load data files
end

end

