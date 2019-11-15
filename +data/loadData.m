function [block, Timeline] = loadData(varargin)
% 30 Jan 2018: This shit was written by LEW
% 14 Nov 2019: Updated to remove rigbox dat dependencies

%check input arguments
if nargin >= 3
	mouseName = varargin{1};
	expDate = varargin{2};
	expNum = varargin{3};
	
	if ischar(expNum)
		expNum = str2double(expNum);
	end 

    [expRef, expLog] = data.constructExpRef(mouseName,expDate,expNum);
else
    error('Check your inputs: mouseName, expDate, expNum')
end 

% load all possible data locations
dataLocations = data.dataPaths();

% load the block file
for i = 1:length(dataLocations)
    blockFilePath = data.makeExpFilePath(expRef, expLog, dataLocations{i}, 'block');
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
    timelineFilePath = data.makeExpFilePath(expRef, expLog, dataLocations{i}, 'timeline');
    try
        load(timelineFilePath);
        if exist('Timeline')
            break
        end
    catch
    end
end

end