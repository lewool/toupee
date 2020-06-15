function paths = getPathsTemplate()
% add as many data locations to this list as you want. put your primary one
% first to break the search-loop faster

paths.server = {...
    {'\\znas.cortexlab.net\Subjects'},...
    {'\\zserver.cortexlab.net\Data\Subjects'},... 
    {'\\zubjects.cortexlab.net\Subjects'}}; 

paths.local = {...
    {'C:\Users\username\local_data\subjects'}}; 

end