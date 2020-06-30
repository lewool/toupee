function [expInfo, fdata] = loadDatafile(expInfo, files)
% Loads data from files into the `expInfo` table
% 
%
% Inputs:
% -------
% expInfo : table
%   A table containing relevant information and data variables (columns) 
%   for particular experiment sessions (rows).
%
% files : char array OR cell array
%   File(s) to load into `expInfo`. If loading files for multiple sessions,
%   use a nested cell array for each session. The elements in the innermost
%   cells can be:
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
%   The updated `expInfo`.
%
% fdata: table
%   The loaded data from `files`. Each row corresponds to an experiment
%   session and contains loaded data from `files` in its columns.
%
%
% Examples:
% ---------
% 1) For a single session: load the block file.
%   details = {'LEW031', '2020-02-03', 1};
%   expInfo = toupee.meta.processExperiment(details);
%   files = 'block';
%   expInfo = toupee.meta.loadDatafile(expInfo, files);
%
% 2) For multiple sessions: load each session's block + timeline files.
%   details = {{'LEW031', '2020-02-03', 1},... 
%              {'LEW037', '2020-03-13', 1},...
%              {'LEW005', '2018-06-10', 2, [2 3]}};
%   expInfo = toupee.meta.processExperiment(details);
%   files = {'block', 'timeline'};
%   expInfo = toupee.meta.loadDatafile(expInfo, files);
%
% 3) For multiple sessions: for the first session load just the timeline
% file, for the second session load the block file and the raw reward valve 
% data from timeline, and for the third session load just the block file.
%   details = {{'LEW031', '2020-02-03', 1},... 
%              {'LEW037', '2020-03-13', 1},...
%              {'LEW005', '2018-06-10', 2, [2 3]}};
%   expInfo = toupee.meta.processExperiment(details);
%   files = {{'timeline'}, {'block', 'rewardvalve.raw.npy'}, {'block'}};
%   expInfo = toupee.meta.loadDatafile(expInfo, files);
%
%
% See Also:
% ---------
% toupee.meta.processExperiment
% toupee.meta.getPaths
%
%
% @todo make better specs for loading neural data
% @todo distinguish neural vs. behavioral datafiles
% @todo add support for binary (+ other?) file types
%

%% Prerun checks.
% Import all other functions in this subpackage and `+misc`.
import toupee.meta.*
import toupee.meta.npy.*
import toupee.misc.*
% Ensure `files` is cell or char.
if ~(iscell(files) || ischar(files))
    error('toupee:meta:loadDatafile:badInput',...
          'The "files" input arg must be a cell or char array')
% Convert to nested cell if not already.
elseif ischar(files)
    files = {{files}};
elseif ~iscell(files{1})
    files = {files};
end
% If there are multiple experiment sessions, repmat `files` if necessary.
nE = size(expInfo, 1);  % number of experiment sessions
if nE > 1 && ~(numel(files) > 1)
    files = repmat(files, [1, nE]);
end

%% Load datafiles for each experiment session.
% Initialize `fdata`.
fdata = table();
% Get paths to servers to search for datafiles.
allPaths = [getPaths().local, getPaths().server];
for iE = 1:nE  % for each experiment session
    subject = expInfo.('subject'){iE};
    expDate = expInfo.('expDate'){iE};
    expNum = expInfo.('expNum'){iE};
    expRef = expInfo.('Row'){iE};
    f = files{iE};  % files to be loaded for current experiment session
    % Get all possible directories of datafiles for this session.
    eDir = cellfun(@(p) fullfile(p, subject, expDate, num2str(expNum)),...
                   allPaths, 'uni', 0);
   % Get full paths to block and/or timeline files if specified. 
   if any(strcmpi(f, 'block'))
       f{strcmpi(f, 'block')} = strcat(expRef, '_Block.mat');
   end
   if any(strcmpi(f, 'timeline'))
       f{strcmpi(f, 'timeline')} = strcat(expRef, '_Timeline.mat');
   end
   % Load any specified misc individual data files.
   % Create full paths for files in `f`.
   fullPaths = cellfun(@(x) fullfile(eDir, x), f, 'uni', 0);
   % Find the first path within each of the possible paths for each file.
   pathIdxs =...
       cellfun(@(x) find(cell2mat(cellfun(@(y) isfile(y), x, 'uni', 0)),...
                         1), fullPaths, 'uni', 0);
   % Remove all paths for files not found.
   emptyPathIdxs = cellfun(@(x) isempty(x), pathIdxs);
   fullPaths(emptyPathIdxs) = [];
   pathIdxs(emptyPathIdxs) = [];
   finalPaths = cellfun(@(x, y) x{y}, fullPaths, pathIdxs, 'uni', 0);
   % Try to load data from files.
   % the loaded data from the datafiles for current expRef
   fdataE = cellfun(@(x) loadMiscFile(x), finalPaths, 'uni', 0);
   % If we loaded some data, then clean file names and assign to `expInfo`
   % and `fdata`
   if ~all(cellfun(@(x) isempty(x), fdataE))
       % Convert any structs in `fdataE` to tables, and cellify all data
       % vars within tables
       sIdxs = cellfun(@(x) isstruct(x), fdataE);
       if any(sIdxs)
           fdataE(sIdxs) = cellfun(@(x) struct2tableNested(x),...
                                   fdataE(sIdxs), 'uni', 0);
           fdataE(sIdxs) = cellfun(@(x) cellifyTableVarsNested(x),...
                                   fdataE(sIdxs), 'uni', 0);
       end
       % Ensure the filenames are table compatible, and add the files'
       % data to `expInfo`.
       [~, fnames, exts] =...
           cellfun(@(x) fileparts(x), finalPaths, 'uni', 0);
       colNames =...
           cellfun(@(x) cleanFilename(x, expRef), fnames, 'uni', 0);
       % Add 'File' suffix to each column name.
       colNames = cellfun(@(x) strcat(x, 'File'), colNames, 'uni', 0);
       % Add the data to the `expInfo` and `fdata` tables.
       expInfo(iE, colNames) = fdataE;
       fdata(iE, colNames) = fdataE;
       % Remove loaded files from `f`.
       rmIdxs = cellfun(@(x) any(strcmpi(x, strcat(fnames, exts))), f);
       f(rmIdxs) = [];
   end
   
   % Mention any files that weren't able to be found/loaded.
   if ~isempty(f)
       fprintf('\nThe following files for %s were unable to be found:\n',...
           expRef);
       disp(f);
   end
   
end

end


function x = loadMiscFile(filepath)
% Tries to load a single .npy or .mat datafile
%
%
% Inputs:
% -------
% filepath : char array
%   The path to a single .npy or .mat file
%
%
% Outputs:
% --------
% x : struct
%   Contains the loaded data from `filepath`
%
%
% Examples:
% ---------
% 1) Load an .npy file
%   x = loadMiscFile('path\to\numpy_file.npy');
%
%
% See Also:
% ---------
% load
% readNPY
%

x = [];  % initialize as empty
if isfile(filepath)  % ensure file exists
    [~, ~, ext] = fileparts(filepath);  % get file extension
    try
        if strcmp(ext, '.mat')  % use `load` if .mat
            loadFn = @load;
        elseif strcmp(ext, '.npy')  % use `readNPY` if .npy
            loadFn = @readNPY;
        end
        fprintf('\nLoading %s ...', filepath);
        x = loadFn(filepath);
        if isequal(loadFn, @load)  % then get out variable from struct
            fieldname = fieldnames(x);
            x = x.(fieldname{1});
        end
        fprintf('\nDone.\n');
    catch ex
        fprintf('\nCould not load %s. Full error message: %s\n',...
             filepath, ex.message);
        return
    end
end

end


function clean = cleanFilename(dirty, expRef)
% Cleans a filename to make it compatible as a struct field or table column
%
%
% Inputs:
% -------
% dirty : char array
%   The filename.
% 
% expRef : char array
%   The expRef for the given file.
%
%
% Outputs:
% --------
% clean : char array
%   The cleaned filename.
%
%
% Examples:
% ---------
% 1) Clean filename that begins with numbers and contains dashes:
%   clean = cleanFilename('2020-02-28_1_LEW032_eye',...
%                         '2020-02-28_1_LEW032');
%
% See Also:
% ---------
% genvarname
%

clean = dirty;  % set `clean` to `dirty`, then clean it.
% Remove expRef from filename.
clean = erase(clean, expRef);
% Replace `.` & '-' with `_`.
clean = strrep(strrep(clean, '.', '_'), '-', '_');
% Remove first letter if incompatible (underscore or number).
while strcmp(clean(1), '_') || ~isnan(str2double(clean(1)))
    clean(1) = [];
end

end
