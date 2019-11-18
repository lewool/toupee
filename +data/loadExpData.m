function expInfo = loadExpData(expInfo)
% 30 Jan 2018: This shit was written by LEW
% 14 Nov 2019: Updated to remove rigbox dat dependencies

%check input arguments
if nargin > 1
    error('Check your inputs: mouseName, expDate, expNum')
else
	mouseName = expInfo.mouseName;
	expDate = expInfo.expDate;
	expNum = expInfo.expNum;
	
	if ischar(expNum)
		expNum = str2double(expNum);
        expInfo.expNum = expNum;
	end 

    [expRef, expLog] = data.constructExpRef(mouseName,expDate,expNum);
end 

% load all possible data locations
paths = data.dataPaths();
dataLocations = paths.server;

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

if ~exist('block')
    warning('No block file was found')
    block = [];
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

if ~exist('Timeline')
    warning('No timeline file was found')
    Timeline = [];
end

expInfo.block = block;
expInfo.Timeline = Timeline;

end

