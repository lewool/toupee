function [expInfo, fdata] = processExperiment(deats, files)
% Gets experiment information and data for given session(s)
% 
%
% Inputs:
% -------
% deats : char array OR cell array
%   The session(s) details. Details can either be provided for a single
%   session or multiple sessions. Details for each individual session
%   should be in the form of either: 1) a char array of the session's 
%   Rigbox expRef, or 2) a cell array of three or four cells that are, in
%   order: i) subject char array; ii) expDate char array; iii) expNum int
%   scalar; iv) (optionally) expSeries int array. Examples below:
%   Single session from expRef:
%   '<expRef>'
%   Single session from subject, expDate, expNum:
%   {'<subject>', '<expDate>', '<expNum>'}
%   Multiple sessions from expRef:
%   {'<expRef1>',...'<expRefN>'}
%   Multiple sessions from subject, expDate, expNum:
%   {{'<subject1>', '<expDate1>', '<expNum1>'},...
%    {'<subjectN>', '<expDateN>', '<expNum2N'}}
%
% files : char array OR cell array
%   File(s) to load into the `expInfo` table. If loading files for
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
% expInfo : table
%   Each row contains fields with information for a particular session. Row
%   names are expRefs. The column names can include: 'subject', 'expDate',
%   'expNum', 'expRef', 'expSeries', 'numPlanes', 'numChannels',
%   'behavioralData', 'neuralData', and additional columns whose names end
%   in 'File' and which contain data loaded from raw datafiles for the 
%   given session.
%
% fdata: table
%   Each row contains the loaded datafiles specified in `files` for the 
%   corresponding experiment session.
%
%
% Examples:
% ---------
% 1) For a single session: return the bare experiment info from a given
% expRef:
%   deats = {'2020-02-03_1_LEW031'};
%   expInfo = toupee.meta.processExperiment(deats);
%
% 2) For a single session: return the bare experiment info from a given
% subject, experiment date, and experiment number.
%   deats = {'LEW031', '2020-02-03', 1};
%   expInfo = toupee.meta.processExperiment(deats);
%
% 3) For multiple sessions: return each session's respective experiment
% info and block and timeline files from given expRefs:
%   deats = {'2020-02-03_1_LEW031',...
%            '2020-03-13_1_LEW037',...
%            '2018-06-10_2_LEW005'};
%   files = {'block', 'timeline'};
%   expInfo = toupee.meta.processExperiment(deats, files);
%
% 4) For multiple sessions: for the first session load just the timeline
% file, for the second session load the block file and the raw reward valve 
% data from timeline, and for the third session load just the block file,
% from given subjects, expDates, and expNums:
%   deats = {{'LEW031', '2020-02-03', 1},... 
%            {'LEW037', '2020-03-13', 1},...
%            {'LEW005', '2018-06-10', 2, [2 3]}};
%   files = {{'timeline'}, {'block', 'rewardvalve.raw.npy'}, {'block'}};
%   expInfo = toupee.meta.processExperiment(deats, files);
%
%
% See Also:
% ---------
% toupee.meta.isExpRef
% toupee.meta.constructExpRef
% toupee.meta.deconstructExpRef
% toupee.meta.loadDatafile
%

%% Prerun checks.
% Import all other functions in this subpackage and `iif`.
import toupee.meta.*
import toupee.misc.iif

% See if given 1) expRef(s), or 2) subject(s), expDate(s), and expNum(s).
% Check to make sure good expRefs can be made from all info in `deats`.
% Default assumptions are 1) provided with expRefs, 2) the provided 
% expRefs, or the details provided to construct expRefs, are good.
fromExpRefs = true;  
goodDeats = true;
% Check if given single expRef.
if ischar(deats) && isExpRef(deats)
    deats = {deats};  % convert to cell array
    nE = 1;  % number of experiment sessions
% check if given multiple expRefs    
elseif iscell(deats) && isExpRef(deats{1})  
    allExpRefs = cellfun(@(x) isExpRef(x), deats);
    if ~all(allExpRefs)
        goodDeats = false; 
    end
    nE = numel(deats);
% Check if given single subject, expDate, expNum.
elseif iscell(deats) && ~isExpRef(deats{1}) && ~iscell(deats{1})...
       && numel(deats) >= 3
    if ~isExpRef(constructExpRef(deats{1}, deats{2}, deats{3}))
        goodDeats = false;
    else
        fromExpRefs = false;
    end
    deats = {deats};  % convert to nested cell array
    nE = 1;
% Check if given multiple subject, expDate, expNum.
elseif iscell(deats) && ~isExpRef(deats{1}) && iscell(deats{1})...
       && numel(deats{1}) >= 3
   allExpRefs = cellfun(@(x) isExpRef(constructExpRef(x{1}, x{2}, x{3})),...
                        deats);
   if ~all(allExpRefs)
        goodDeats = false; 
   else
       fromExpRefs = false;
   end
   nE = numel(deats);
else  % any other aspect of `deats` is wrong
    goodDeats = false; 
end

% Throw error if `deats` not properly provided.
if ~goodDeats
    error('toupee:meta:processExperiment:badInput',...
          ['The "deats" input arg was not provided correctly. ',...
           'See the docstring for info and examples.']);
end

%% Set up `expInfo`.
% Make row names expRefs.
rowNames = iif(fromExpRefs, deats,...
               @() cellfun(@(d) constructExpRef(d{1}, d{2}, d{3}), deats,...
                           'uni', 0));
% Specify columns and create `expInfo` table.
columnNames = {'subject', 'expDate', 'expNum', 'expSeries',...
               'nPlanes', 'nChannels', 'behavioralData', 'neuralData'};
columnTypes = {'cell', 'cell', 'cell', 'cell',...
               'cell', 'cell', 'table', 'table'};
nC = numel(columnNames);  % number of columns
expInfo = table('Size', [nE, nC], 'VariableNames', columnNames,...
                'VariableTypes', columnTypes, 'RowNames', rowNames);
% Add info from `deats` to `expInfo`.
if fromExpRefs
    [subjects, expDates, expNums] =...
        cellfun(@(x) deconstructExpRef(x), deats, 'uni', 0);
    expInfo.subject = subjects';
    expInfo.expDate = expDates';
    expInfo.expNum = expNums';
else
    expInfo.subject = cellfun(@(d) d{1}, deats, 'uni', 0)';
    expInfo.expDate = cellfun(@(d) d{2}, deats, 'uni', 0)';
    expInfo.expNum = cellfun(@(d) d{3}, deats, 'uni', 0)';
    expInfo.expSeries =...
        cellfun(@(d) iif(numel(d) > 3, @() d{4}, []), deats, 'uni', 0)';
end

%% Load datafiles if specified.
if nargin > 1
    [expInfo, fdata] = loadDatafile(expInfo, files);  % load datafiles
end

end

