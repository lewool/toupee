function [block, Timeline] = loadData(varargin)
% This shit was written by LEW on 30 Jan 2018

%%

if nargin >= 3
	mouseName = varargin{1};
	expDate = varargin{2};
	expNum = varargin{3};
	
	if ischar(expNum)
		expNum = str2double(expNum);
	end 

	expRef = dat.constructExpRef(mouseName,expDate,expNum);
    
elseif nargin == 2
	expRef = varargin{1};
    mouseName = varargin{2};

elseif nargin == 1
	expRef = varargin{1};
end 

oldDataFolder = '\\zserver.cortexlab.net\Data\Subjects'; 
newDataFolder = '\\zubjects.cortexlab.net\Subjects'; 

try
    blockFilePath = dat.expFilePath(expRef, 'block', 'master');
    load(blockFilePath);
    
catch
    try
        load(strrep(blockFilePath, oldDataFolder, newDataFolder));
    catch %if no TL file
        block = [];
    end

end

try 
    tlFilePath = dat.expFilePath(expRef, 'timeline', 'master');
    load(tlFilePath);
catch
    try
        load(strrep(tlFilePath, oldDataFolder, newDataFolder));
    catch %if no TL file
        Timeline = [];
    end
end

end