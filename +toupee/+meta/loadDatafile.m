function expInfo = loadDatafile(expInfo, file)
% Loads data from a file into the `expInfo` struct
% 
%
% Inputs:
% -------
% expInfo : struct
%   A struct containing relevant information and data for particular
%   experiment sessions.
%
% file : cell array
%   File(s) to load into the `expInfo` struct. If loading files for
%   multiple sessions, use a nested cell array for each session. The
%   elements in the innermost cells can be:
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
%   A struct containing relevant information and data for particular
%   experiment sessions.
%
%
% Examples:
% ---------
% 1) Load the block file for a single session from a single subject.
%   details = {'LEW031', '2020-02-03', 1};
%   expInfo = toupee.meta.processExperiment(details);
%   file = {'block'};
%   expInfo = toupee.meta.loadDatafile(expInfo, file);
%
% 2) Load the block + timeline files for multiple sessions from multiple
% subjects.
%   details = {{'LEW031', '2020-02-03', 1},... 
%              {'LEW037', '2020-03-13', 1},...
%              {'LEW005', '2018-06-10', 2, [2 3]}};
%   expInfo = toupee.meta.processExperiment(details);
%   files = {'block', 'timeline'};
%   expInfo = toupee.meta.loadDatafile(expInfo, files);
%
% 3) Load the block + specific individual data files dependent on the
% particular session.
%   details = {{'LEW031', '2020-02-03', 1},... 
%              {'LEW037', '2020-03-13', 1},...
%              {'LEW005', '2018-06-10', 2, [2 3]}};
%   expInfo = toupee.meta.processExperiment(details);
%   files = {{'wheel.position.npy', 'wheel.timestamps.npy'},...
%             {'block', 'rewardvalve.raw.npy'}, {'block'}};
%   expInfo = toupee.meta.processExperiment(expInfo, files);
%
%
% See Also:
% ---------
% `toupee.meta.processExperiment`
% `toupee.meta.getPaths`
%
% @todo make better specs for loading neural data
% @todo distinguish neural vs. behavioral datafiles
% @todo add support for binary (+ other?) file types

% import all other functions in this subpackage and `iif`.
import toupee.meta.*
import toupee.meta.npy.*
import toupee.misc.iif

% Do some checks on input args.
% Ensure `file` is cell.
if ~iscell(file)
    error('toupee:meta:loadDatafile:badInput',...
          'The "file" input arg should be a cell array')
elseif ~iscell(file{1})  % convert to nested cell if not already
    file = {file};
end
% If there are multiple sessions, repmat `file` if necessary.
if numel(expInfo) > 1 && ~(numel(file) > 1)
    file = repmat(file, [1, numel(expInfo)]);
end

% For each experiment session, load datafiles.
for e = 1:length(expInfo)
    subject = expInfo(e).subject;
    expDate = expInfo(e).expDate;
    expNum = expInfo(e).expNum;
    expRef = expInfo(e).expRef;
    allPaths = [getPaths().server, getPaths().local];
    f = file{e};  % files to be loaded for current experiment session
    % Check each location for exp data.
    for i = 1:numel(allPaths)
        p = allPaths{i};
        % Get directory where datafiles would be for this location.
        eDir = fullfile(p, subject, expDate, num2str(expNum));
        % Load block file if specified
        if any(strcmpi(f, 'block'))  
            blockFilePath = fullfile(eDir, strcat(expRef, '_Block.mat'));
            if isfile(blockFilePath)  % load file and remove from `f`
                fprintf('\nLoading %s...', blockFilePath);
                block = load(blockFilePath);
                expInfo.block = block.block;
                fprintf('\nDone.');
                f(strcmpi(f, 'block')) = [];
            end
        end
        % Load timeline file if specified
        if any(strcmpi(s, 'timeline'))  
            timelineFilePath =... 
                fullfile(eDir, strcat(expRef, '_Timeline.mat'));
            if isfile(timelineFilePath)  % load file and remove from `f`
                fprintf('\nLoading %s...', timelineFilePath);
                timeline = load(timelineFilePath);
                expInfo.timeline = timeline.Timeline;
                fprintf('\nDone.');
                f(strcmpi(f, 'timeline')) = [];
            end
        end
        % Load any specified misc individual data files.
        % Create full paths for files in `f`.
        fullPaths = cellfun(@(x) fullfile(eDir, x), f, 'UniformOutput', 0);
        % Try to load data from files.
        fdata =...
            cellfun(@(x) loadMiscFile(x), fullPaths, 'UniformOutput', 0);
        if ~isempty(fdata)
            % Remove empty values for files data wasn't loaded from.
            nada = cellfun(@(x) isempty(x), fdata);
            fdata(nada) = [];
            loadedFiles = f(~nada);
            [~, fnames, ~] = cellfun(@(x) fileparts(x), loadedFiles,...
                                     'UniformOutput', 0);
            fnames = cellfun(@(x) strrep(x, '.', '_'), fnames,...
                             'UniformOutput', 0);  % replace `.` with `_`
            % Add data to `expInfo`.
            for i = 1:fnames
                expInfo(e).behavioralData.(fnames{i}) = fdata{i};
            end
            f(strcmpi(f, loadedFiles)) = [];  % rm loaded files from `f`.
        end
    end
    % Mention any files that weren't able to be found/loaded.
    if ~isempty(f)
        fprintf('The following files for %s were unable to be found:',...
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

x = [];  % initialize as empty
if isfile(filepath)  % ensure file exists
    [~, ~, ext] = fileparts(filepath);  % get file extension
    try
        if strcmp(ext, '.mat')  % use `load` if .mat
            loadFn = @load;
        elseif strcmp(ext, '.npy')  % use `readNPY` if .npy
            loadFn = @readNPY;
        end
        fprintf('\nLoading %s...', filepath);
        x = loadFn(filepath);
        if isequal(loadFn, @load)  % then get out variable from struct
            fieldname = fieldnames(x);
            x = x.(fieldname{1});
        end
        fprintf('\nDone.');
    catch ex
        fprintf('\nCould not load %s. Full error message: %s',...
             filepath, ex.message);
        return
    end
end

end