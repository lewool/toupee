function [block, Timeline] = loadData(varargin)
% 30 Jan 2018: This shit was written by LEW
% 14 Nov 2019: Updated to remove rigbox dat dependencies

%%

if nargin >= 3
	mouseName = varargin{1};
	expDate = varargin{2};
	expNum = varargin{3};
	
	if ischar(expNum)
		expNum = str2double(expNum);
	end 

    [expRef, expLog] = constructExpRef(mouseName,expDate,expNum);
else
    error('Check your inputs: mouseName, expDate, expNum')
end 

% add as many data locations to this list as you want. put your primary one
% first to break the search-loop faster
dataLocations = {...
    {'\\znas.cortexlab.net\Subjects'},...
    {'\\zserver.cortexlab.net\Data\Subjects'},... 
    {'\\zubjects.cortexlab.net\Subjects'}}; 

% load the block file
for i = 1:length(dataLocations)
    blockFilePath = makeExpFilePath(expRef, expLog, dataLocations{i}, 'block');
    try
        load(blockFilePath);
        if exist('block')
            break
        end
    catch
    end
end

% load the timeline file
for i = 1:length(dataLocations)
    timelineFilePath = makeExpFilePath(expRef, expLog, dataLocations{i}, 'timeline');
    try
        load(timelineFilePath);
        if exist('Timeline')
            break
        end
    catch
    end
end

end

function [expRef, expLog] = constructExpRef(mouseName,expDate,expNum)
    expRef = strcat(expDate,'_',num2str(expNum),'_',mouseName);
    expLog = {{mouseName};{expDate};{num2str(expNum)}};
end

function expFilePath = makeExpFilePath(expRef, expLog, dataFolder, fileType)
    fileType(1) = upper(fileType(1));
    suffix = strcat('_',fileType,'.mat');
    expFilePath = char(fullfile(dataFolder,expLog{1},expLog{2},expLog{3},strcat(expRef,suffix)));
end