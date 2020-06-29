function paths = getPathsTemplate()
% A template file for returning the paths to subject data.
% A new file named `getPaths` should be created from this template,
% and the code should be edited to point to the relevant subject data
% locations.
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
%   [paths] = getPaths();
%

% Add as many data locations to these fields as desired. Put primary
% locations first in order to break the search-loop faster.

paths.local = {...
    {'C:\Users\username\local_data\subjects'}}; 

paths.server = {...
    {'\\znas.cortexlab.net\Subjects'},...
    {'\\zserver.cortexlab.net\Data\Subjects'},... 
    {'\\zubjects.cortexlab.net\Subjects'}}; 

end