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

newDataFolder = '\\zserver.cortexlab.net\Data\Subjects'; 
oldDataFolder = '\\zubjects.cortexlab.net\Subjects'; 

blockFilePath = dat.expFilePath(expRef, 'block', 'master');
tlFilePath = dat.expFilePath(expRef, 'timeline', 'master');

try
    load(blockFilePath);    
    load(tlFilePath);
    
catch
    try
        load(strrep(blockFilePath, oldDataFolder, newDataFolder));
        load(strrep(tlFilePath, oldDataFolder, newDataFolder));
    catch %if no TL file
        block = block;
        Timeline = [];
    end

end

end