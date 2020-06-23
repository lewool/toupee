function [expInfo, fdata] = loadDatafile(expInfo, files)
% Loads data from a file into the `expInfo` struct
% 
%
% Inputs:
% -------
% expInfo : struct array
%   A struct containing relevant information and data for particular
%   experiment sessions.
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
%   A struct containing relevant information and data for particular
%   experiment sessions.
%
% fdata: struct array
%   Each array element contains the loaded datafiles specified in `files`
%   for the corresponding experiment session.
%
% Examples:
% ---------
% 1) For a single session: load the block file.
%   details = {'LEW031', '2020-02-03', 1};
%   expInfo = toupee.meta.processExperiment(details);
%   file = {'block'};
%   expInfo = toupee.meta.loadDatafile(expInfo, file);
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
% `toupee.meta.processExperiment`
% `toupee.meta.getPaths`
%
% @todo make better specs for loading neural data
% @todo distinguish neural vs. behavioral datafiles
% @todo add support for binary (+ other?) file types

% import all other functions in this subpackage.
import toupee.meta.*
import toupee.meta.npy.*
import toupee.misc.iif

% Do some checks on input args.
if ~(iscell(files) || ischar(files))  % ensure `files` is cell or char
    error('toupee:meta:loadDatafile:badInput',...
          'The "file" input arg must be a cell array')
% convert to nested cell if not already
elseif ischar(files)
    files = {{files}};
elseif ~iscell(files{1})
    files = files{1};
end
% If there are multiple sessions, repmat `file` if necessary.
if numel(expInfo) > 1 && ~(numel(files) > 1)
    files = repmat(files, [1, numel(expInfo)]);
end

% For each experiment session, load datafiles.
for e = 1:length(expInfo)
    subject = expInfo(e).subject;
    expDate = expInfo(e).expDate;
    expNum = expInfo(e).expNum;
    expRef = expInfo(e).expRef;
    allPaths = [getPaths().server, getPaths().local];
    f = files{e};  % files to be loaded for current experiment session
    % Check each location for exp data.
    for loc = 1:numel(allPaths)
        p = allPaths{loc};
        % Get directory where datafiles would be for this location.
        eDir = fullfile(p, subject, expDate, num2str(expNum));
        % Load block file if specified
        if any(strcmpi(f, 'block'))  
            blockFilePath = fullfile(eDir, strcat(expRef, '_Block.mat'));
            if isfile(blockFilePath)  % load file and remove from `f`
                fprintf('\nLoading %s...', blockFilePath);
                block = load(blockFilePath);
                expInfo(e).block = block.block;
                fdata(e).block = expInfo(e).block;
                fprintf('\nDone.\n');
                f(strcmpi(f, 'block')) = [];
            end
        end
        % Load timeline file if specified
        if any(strcmpi(f, 'timeline'))  
            timelineFilePath =... 
                fullfile(eDir, strcat(expRef, '_Timeline.mat'));
            if isfile(timelineFilePath)  % load file and remove from `f`
                fprintf('\nLoading %s...', timelineFilePath);
                timeline = load(timelineFilePath);
                expInfo(e).timeline = timeline.Timeline;
                fdata(e).timeline = expInfo(e).timeline;
                fprintf('\nDone.\n');
                f(strcmpi(f, 'timeline')) = [];
            end
        end
        % Load any specified misc individual data files.
        % Create full paths for files in `f`.
        fullPaths = cellfun(@(x) fullfile(eDir, x), f, 'UniformOutput', 0);
        % Try to load data from files.
        fdataE =...  % the loaded data from the datafiles for current `e`
            cellfun(@(x) loadMiscFile(x), fullPaths, 'UniformOutput', 0);
        if ~all(cellfun(@(x) isempty(x), fdataE))
            % Remove empty values for files data wasn't loaded from.
            nada = cellfun(@(x) isempty(x), fdataE);
            fdataE(nada) = [];
            loadedFiles = f(~nada);
            [~, fnames, ~] = cellfun(@(x) fileparts(x), loadedFiles,...
                                     'UniformOutput', 0);
            % Ensure the fieldname is struct compatible, and add the file's
            % data to `expInfo`.
            for i = 1:numel(fnames)
                % Remove expRef from fieldname.
                if contains(fnames{i}, expRef)
                    unders = strfind(fnames{i}, '_');
                    fnames{i} = fnames{i}(unders(end) + 1 : end);
                end
                % Replace `.` & '-' with `_`.
                fnames{i} = strrep(strrep(fnames{i}, '.', '_'), '-', '+');
                % Add to `expInfo`.
                expInfo(e).behavioralData.(fnames{i}) = fdataE{i};
                fdata(e).(fnames{i}) = fdataE{i};
            end
            % Remove loaded files from `f`.
            f(strcmpi(f, loadedFiles)) = [];
        end
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
        fprintf('\nDone.\n');
    catch ex
        fprintf('\nCould not load %s. Full error message: %s\n',...
             filepath, ex.message);
        return
    end
end

end