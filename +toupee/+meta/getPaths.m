function paths = getPaths()
% Returns the paths to subject data.
%
% Outputs:
% --------
% 'paths' : struct
%   Contains the cell array fields 'server' and 'local' for subject data
%   path(s) on remote server(s) and on the local computer.
%
% Examples:
% ---------
% 1) Get subject data paths:
%   paths = getPaths();

% Add as many data locations to these fields as you want. Put primary
% locations first in order to break the search-loop faster.

paths.local = {...
    'C:\Users\Jai\m2v1_expected_reward\data\subjects'}; 

paths.server = {...
    '\\znas.cortexlab.net\Subjects',...
    '\\zserver.cortexlab.net\Data\Subjects',... 
    '\\zubjects.cortexlab.net\Subjects'}; 

end