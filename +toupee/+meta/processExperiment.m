function [expInfo, fdata] = processExperiment(details, files)
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
% files : cell array OR char array
%   File(s) to load into the `expInfo` struct. If loading files for
%   multiple sessions, use a nested cell array for each session. The
%   elements in the innermost cells can be:
%   1) 'block': loads the block file
%   2) 'timeline': loads the timeline file
%   3) Full names of individual behavioral or neural data files. These
%   files can correspond to block data (e.g. 'wheel.position.npy'),
%   timeline data (e.g. 'rewardvalve.raw.npy'), or suite2P data (e.g. 
%   'SVD_plane1.mat')
% 
%
% Outputs:
% --------
% expInfo : struct array
%   Each array element contains fields with information for a particular
%   session. The fields for each struct element can include:
%   'subject', 'expDate', 'expNum', 'expRef', 'expSeries', 'block',
%   'timeline', 'numPlanes', 'numChannels', 'behavioralData', 'neuralData'.
%
% fdata: struct array
%   Each array element contains the loaded datafiles specified in `files`
%   for the corresponding experiment session.
%
% Examples:
% ---------
% 1) For a single session: return the bare experiment info.
%   details = {'LEW031', '2020-02-03', 1};
%   expInfo = toupee.meta.processExperiment(details);
%
% 2) For multiple sessions: return each session's respective experiment
% info and block and timeline files.
%   details = {{'LEW031', '2020-02-03', 1},... 
%              {'LEW037', '2020-03-13', 1},...
%              {'LEW005', '2018-06-10', 2, [2 3]}};
%   files = {'block', 'timeline'};
%   expInfo = toupee.meta.processExperiment(details, files);
%
% 3) For multiple sessions: for the first session load just the timeline
% file, for the second session load the block file and the raw reward valve 
% data from timeline, and for the third session load just the block file.
%   details = {{'LEW031', '2020-02-03', 1},... 
%              {'LEW037', '2020-03-13', 1},...
%              {'LEW005', '2018-06-10', 2, [2 3]}};
%   files = {{'timeline'}, {'block', 'rewardvalve.raw.npy'}, {'block'}};
%   expInfo = toupee.meta.processExperiment(details, files);
%
%
% See Also:
% ---------
% toupee.meta.constructExpRef
% toupee.meta.loadDatafile
%

% Import all other functions in this subpackage and `iif`.
import toupee.meta.*
import toupee.misc.iif

% Do some checks on input args.
if ~iscell(details)  % ensure `details` is cell.
    error('toupee:meta:processExperiment:badInput',...
          'The "details" input arg must be a cell array')
elseif ~iscell(details{1})  % convert to nested cell if not already
    details = {details};
end

% For each experiment session, set up `expInfo`.
for e = 1:numel(details)
    % Initialize expInfo struct.
    d = details{e};
    if numel(d) < 3 || numel(d) > 4
        error('toupee:meta:processExperiment:badInput',...
              ['Each session specified in the "details" input arg cell '... 
               'array must contain exactly 3 or 4 elements'])
    end
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
end

if nargin > 1  % try to load data files
    if ~(iscell(files) || ischar(files))  % ensure `files` is cell or char
        error('toupee:meta:processExperiment:badInput',...
              'The "files" input arg must be a cell or char array')
    end
    [expInfo, fdata] = loadDatafile(expInfo, files);  % load data files
end

end

